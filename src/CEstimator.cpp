#include <iostream>
#include <fstream>
#include <iomanip>
#include "CEstimator.h"
#include "ConfigParser.h"
#include "TF1.h"
#include "TFormula.h"

#define PI 3.14159265359

using namespace cana;

CEstimator::CEstimator(const std::string &path)
: formula(nullptr), step_range(30)
{
    if(!path.empty())
        LoadFormula(path);
}

CEstimator::~CEstimator()
{
    DeleteFormula();
}

void CEstimator::LoadFormula(const std::string &path)
{
    ConfigParser c_parser;

    if(!c_parser.OpenFile(path)) {
        std::cerr << "CEstimator Error: cannot open file "
                  << path << ", abort function loading."
                  << std::endl;
        return;
    }

    while(c_parser.ParseLine())
    {
        // first line will be the function
        if(c_parser.NbofElements() > 0)
            break;
    }

    std::string function;
    // glue all the parts together
    // since there probably is space in the function expression
    // which can be used as the splitter for parameters
    for(auto &ele : c_parser.TakeAll())
    {
        function += ele.String();
    }

    SetFormula(function.c_str());

    size_t idx = 0;
    // continue reading the parameters
    while(c_parser.ParseLine())
    {
        if(!c_parser.NbofElements())
            continue;

        if(idx < parameters.size())
            parameters[idx++] = Parameter(c_parser.TakeFirst().Double());
        else
            std::cout << "CEstimator Warning: Reading "
                      << idx << " parameters, but the function only accepts "
                      << parameters.size()
                      << std::endl;
    }

    UpdatePars();
    c_parser.CloseFile();
}

void CEstimator::SaveFormula(const std::string &path)
{
    UpdatePars();

    std::ofstream output(path);
    // write down function
    output << formula->GetExpFormula().Data() << std::endl;

    // write down parameters
    for(int i = 0; i < formula->GetNpar(); ++i)
        output << std::setw(12) << formula->GetParameter(i) << std::endl;

    output.close();
}

void CEstimator::DeleteFormula()
{
    if(formula)
        delete formula, formula = nullptr;
    parameters.clear();
    M_penalty_inv = CMatrix();
}

void CEstimator::SetFormula(const char *c)
{
    DeleteFormula();
    // name, expression, add to global list
    // the formula is taken care by this class
    // thus put false to not add it to root global list
    formula = new TFormula("myform", c, false);
    parameters.resize(formula->GetNpar());
}

void CEstimator::SetFormula(TF1 *tf)
{
    DeleteFormula();

    formula = new TFormula("myform", tf->GetExpFormula().Data(), false);
    for(int i = 0; i < tf->GetNpar(); ++i)
        parameters.emplace_back(tf->GetParameter(i));
}

void CEstimator::SetDataPoints(const std::vector<double> &x,
                               const std::vector<double> &y,
                               const std::vector<double> &err)
{
    if(x.size() != y.size() || x.size() != err.size())
    {
        std::cerr << "CEstimator Error: unmatched vector size from SetDataPoints!"
                  << std::endl;
        return;
    }

    data.clear();

    for(size_t i = 0; i < x.size(); ++i)
    {
        data.emplace_back(x[i], y[i], err[i]);
    }

    // automatically set the weight matrix as covariance
    CMatrix covariance_inv(err.size());
    for(size_t i = 0; i < err.size(); ++i)
    {
        covariance_inv(i, i) = 1/err.at(i)/err.at(i);
    }

    M_weight_inv = covariance_inv;
}

void CEstimator::SetDataPoints(const size_t &n, const double *x, const double *y, const double *err)
{
    data.clear();

    for(size_t i = 0; i < n; ++i)
    {
        data.emplace_back(x[i], y[i], err[i]);
    }
}

void CEstimator::SetParameter(const size_t &i, const double &p, const double &step, const double &fine_step)
{
    if(i >= parameters.size())
    {
        std::cerr << "CEstimator Error: Setting "
                  << i << " parameter, while it only has "
                  << parameters.size()
                  << " parameters."
                  << std::endl;
        return;
    }

    Parameter new_p(p, step, fine_step);
    parameters[i] = new_p;
}

void CEstimator::SetParameters(const std::vector<double> &p)
{
    for(size_t i = 0; i < p.size(); ++i)
        SetParameter(i, p.at(i));
}

void CEstimator::SetWeightMatrix(const CMatrix &m)
{
    // Set V, but only V^(-1) is used to evaluate, so do inverse here
    M_weight_inv = m.Inverse();
}

void CEstimator::SetPenaltyMatrix(const CMatrix &m)
{
    // Set Vbeta, but only Vbeta^(-1) is used to evaluate, so do inverse here
    M_penalty_inv = m.Inverse();
}

// fit the function to data
void CEstimator::Fit(int c_iter, int f_iter, bool verbose)
{
    // coarse tune
    int count = 1;
    do
    {
        if(verbose)
            std::cout << std::endl << "*Coarse* step optimization, iteration "
                      << count << ", current chi square is "
                      << GetReducedChiSquare()
                      << std::endl;
    } while(Optimize(step_range, false, verbose) && count++ < c_iter);

    // fine tune
    count = 1;
    do
    {
        if(verbose)
            std::cout << std::endl << "*Fine* step optimization, iteration "
                      << count << ", current chi square is "
                      << GetReducedChiSquare()
                      << std::endl;
        // set step range to be 100 since fine_step is about 1/100 of step
    } while(Optimize(100, true, verbose) && count++ < f_iter);

    if(verbose)
        std::cout << "Fit is done, final chi square is "
                  << GetReducedChiSquare()
                  << std::endl;

    // update the parameters to formula
    UpdatePars();
}


// optimize all the parameters independently
bool CEstimator::Optimize(int steps, bool fine, bool verbose)
{
    bool optimized = false;

    // lock all the parameters first
    for(size_t i = 0; i < parameters.size(); ++i)
        LockPar(i);

    // optimize one by one
    for(size_t i = 0; i < parameters.size(); ++i)
    {
        UnlockPar(i);

        int minimum = 0;
        double eval = Evaluate();
        do
        {
            // update parameters to the minimum
            if(minimum != 0) {
                NextStep(minimum, verbose);
                minimum = 0;
                optimized = true;
            }

            // calculate steps
            CalcStep(fine);

            // find the minimum within range
            for(int i = 1; i <= steps; ++i)
            {
                double this_val_m = Evaluate(-i);
                double this_val_p = Evaluate(i);
                if(this_val_m < eval) {
                    eval = this_val_m;
                    minimum = -i;
                }
                if(this_val_p < eval) {
                    eval = this_val_p;
                    minimum = i;
                }
            }

        } while(minimum != 0);

        LockPar(i);
    }

    return optimized;
}

double CEstimator::Evaluate(const double &factor)
{
    for(size_t i = 0; i < parameters.size(); ++i)
    {
        double value = parameters.at(i).value;

        if(!parameters.at(i).lock)
            value += parameters.at(i).step*factor;

        formula->SetParameter(i, value);
    }

    CMatrix p(1, data.size());
    for(size_t i = 0; i < data.size(); ++i)
        p(0, i) = data.at(i).val - formula->Eval(data.at(i).x);

    double result = p*M_weight_inv*transpose(p);

    // penalty matrix exists, calculate penalty term
    if(M_penalty_inv.DimN() == parameters.size() &&
       M_penalty_inv.DimM() == parameters.size())
    {
        CMatrix b(1, parameters.size());
        for(size_t i = 0; i < parameters.size(); ++i)
            b(0, i) = parameters.at(i).value - parameters.at(i).initial;

        result += b*M_penalty_inv*transpose(b);
    }

    return result;
}

CMatrix CEstimator::GetHessian()
{
    CMatrix J(data.size(), parameters.size());

    for(size_t j = 0; j < J.DimN(); ++j)
    {
        for(size_t i = 0; i < J.DimM(); ++i)
        {
            double gradient = formula->Eval(data.at(j).val);
            formula->SetParameter(i, parameters.at(i).value + parameters.at(i).fine_step);
            gradient -= formula->Eval(data.at(j).val);
            formula->SetParameter(i, parameters.at(i).value);
            gradient /= parameters.at(i).fine_step;
            if(gradient == 0.)
                J(j, i) = -1e-10; // put a small number
            else
                J(j, i) = -gradient;
        }
    }
    CMatrix J_w = M_weight_inv.Cholesky()*J;
    return J_w.Transpose()*J_w;
}

void CEstimator::CalcStep(bool fine)
{
    CMatrix rho_inv(1, parameters.size());
    for(size_t i = 0; i < rho_inv.DimN(); ++i)
    {
        if(fine)
            rho_inv(0, i) = 1/parameters.at(i).fine_step;
        else
            rho_inv(0, i) = 1/parameters.at(i).base_step;
    }

    CMatrix step = rho_inv*GetHessian();
    for(size_t i = 0; i < step.DimM(); ++i)
        parameters[i].step = 1/step(0, i);
}

void CEstimator::NextStep(const double &factor, bool verbose)
{
    // change the parameters according to its step size
    for(size_t i = 0; i < parameters.size(); ++i)
    {
        if(parameters[i].lock)
            continue;
        parameters[i].prev_value = parameters[i].value;
        parameters[i].value += parameters[i].step*factor;
        if(verbose) {
            std::cout << "Parameter " << i << " is updated, "
                      << "previous value: " << parameters[i].prev_value
                      << ", current value: " << parameters[i].value
                      << std::endl;
        }
    }
}

void CEstimator::UpdatePars()
{
    for(size_t i = 0; i < parameters.size(); ++i)
        formula->SetParameter(i, parameters.at(i).value);
}

TFormula *CEstimator::GetFormula()
{
    // make sure its the latest value, since Evaluate may modify it
    UpdatePars();
    return formula;
}

double CEstimator::GetFormulaVal(double x)
{
    return formula->Eval(x);
}

double CEstimator::GetReducedChiSquare()
{
    UpdatePars();
    double result = 0;
    size_t ndof = data.size() - parameters.size();
    for(size_t i = 0; i < data.size(); ++i)
    {
        double sig = data.at(i).error;
        double expect_val = formula->Eval(data.at(i).x);
        double data_val = data.at(i).val;
        if(sig != 0.)
            result += (data_val - expect_val)*(data_val - expect_val)/sig/sig;
        else
            ndof--;
    }

    return result/(double)ndof;
}

double CEstimator::GetPearsonChiSquare()
{
    UpdatePars();
    double result = 0;

    for(size_t i = 0; i < data.size(); ++i)
    {
        double expect_val = formula->Eval(data.at(i).x);
        double data_val = data.at(i).val;
        if(expect_val != 0.)
            result += (data_val - expect_val)*(data_val - expect_val)/expect_val;
    }

    return result;
}

double CEstimator::GetAbsoluteError()
{
    UpdatePars();
    double result = 0;
    for(size_t i = 0; i < data.size(); ++i)
    {
        result += abs(data.at(i).val - formula->Eval(data.at(i).x));
    }
    return result;
}

double CEstimator::GetRootMeanSquaredError()
{
    return sqrt(GetReducedChiSquare());
}

// get negative log likelihood for a Gaussian process
double CEstimator::GetNLL_Gaussian()
{
    UpdatePars();
    double result = 0;

    // get RSS, residual sum of squares
    for(size_t i = 0; i < data.size(); ++i)
    {
        double expect_val = formula->Eval(data.at(i).x);
        double data_val = data.at(i).val;
        result += (data_val - expect_val)*(data_val - expect_val);
    }

    result = log(result/(data.size() - parameters.size()) + 1.);
    return data.size()*(result + log(2.*PI) + 1.) - parameters.size();
}

// Akaike Information Criterion L* + 2m
double CEstimator::GetAkaikeCriterion()
{
    return GetNLL_Gaussian() + 2*parameters.size();
}

// Bayesian Information Criterion L* + mln(n)
double CEstimator::GetBayesianCriterion()
{
    return GetNLL_Gaussian() + parameters.size()*log((double)data.size());
}

// Hannan Criterion L* + cmln(ln(n))
double CEstimator::GetHannanCriterion(const double &c)
{
    if(c < 2) {
        std::cout << "CEstimator Warning: The constant c in Hannan Criterion should "
                  << "be larger or equal to 2."
                  << std::endl;
    }
    return GetNLL_Gaussian() + parameters.size()*log(log((double)data.size()))*c;
}

// Kashyap Criterion L* + mln(n/2pi) + ln|F_M|
// F_M is Fisher information matrix, to be implemented
double CEstimator::GetKashyapCriterion(const CMatrix &F_M)
{
    return GetNLL_Gaussian() + parameters.size()*log((double)data.size()/2./PI) + log(abs(det(F_M)));
}
