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
: formula(nullptr), fine_step_size(1.), coarse_step_size(100.)
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
    M_penalty = CMatrix();
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

    GenerateCovMatrix();
}

void CEstimator::SetDataPoints(const size_t &n, const double *x, const double *y, const double *err)
{
    data.clear();

    for(size_t i = 0; i < n; ++i)
    {
        data.emplace_back(x[i], y[i], err[i]);
    }

    GenerateCovMatrix();
}

// generate the covariance matrix from data
void CEstimator::GenerateCovMatrix()
{
    // automatically set the weight matrix as covariance
    CMatrix covariance_inv(data.size());
    for(size_t i = 0; i < data.size(); ++i)
    {
        if(data.at(i).error != 0.)
            covariance_inv(i, i) = 1/data.at(i).error/data.at(i).error;
        else
            covariance_inv(i, i) = 1e4/data.at(i).val/data.at(i).val; // set 1% level
    }

    SetWeightMatrix(covariance_inv);
}

void CEstimator::SetParameter(const size_t &i, const double &p, const double &step)
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

    Parameter new_p(p, step);
    parameters[i] = new_p;
}

void CEstimator::SetParameters(const std::vector<double> &p)
{
    for(size_t i = 0; i < p.size(); ++i)
        SetParameter(i, p.at(i));
}

void CEstimator::SetWeightMatrix(const CMatrix &m)
{
    M_weight = m;
    if(m.IsDiagonal())
        weight_diagonal = true;
    else
        weight_diagonal = false;
}

void CEstimator::SetPenaltyMatrix(const CMatrix &m)
{
    M_penalty = m;
    if(m.IsDiagonal())
        penalty_diagonal = true;
    else
        penalty_diagonal = false;
}

void CEstimator::SetStepFactor(const double &fine, const double &coarse)
{
    // set fine step
    if(fine == 0.) { // cannot be 0
        fine_step_size = 1.;
    } else {
        fine_step_size = fine;
    }

    // set coarse step
    if(coarse == 0.) { // cannot be 0
        coarse_step_size = 100.;
    } else if (coarse < fine) {
        std::cout << "CEstimator Warning: coarse step size cannot be smaller than "
                  << "fine step size, automatically set to " << fine*10
                  << std::endl;
        coarse_step_size = fine*10.;
    } else {
        coarse_step_size = coarse;
    }
}

// fit the function to data
void CEstimator::Fit(int c_iter, int f_iter, int range, bool verbose)
{
    // coarse tune parameters
    int iter = 0;
    while(iter++ < c_iter)
    {
        if(verbose) {
            std::cout << std::endl << "*Coarse* parameters optimization, iteration "
                      << iter << ", current chi square is "
                      << GetReducedChiSquare()
                      << std::endl;
        }

        // cannot optimize anymore
        if(!Optimize(range, coarse_step_size, verbose))
            break;

    }

    // fine tune parameters
    iter = 0;
    while(iter++ < f_iter)
    {
        if(verbose) {
            std::cout << std::endl << "*Fine* parameters optimization, iteration "
                      << iter << ", current chi square is "
                      << GetReducedChiSquare()
                      << std::endl;
        }

        if(!Optimize(range, fine_step_size, verbose))
            break;

    }

    if(verbose) {
        std::cout << "Fit is done, final chi square is "
                  << GetReducedChiSquare()
                  << std::endl;
    }

    // update the parameters to formula
    UpdatePars();
}


// optimize all the parameters independently
bool CEstimator::Optimize(int range, double factor, bool verbose)
{
    bool optimized = false;

    // lock all the parameters first
    for(size_t i = 0; i < parameters.size(); ++i)
        LockPar(i);

    // optimize one by one
    for(size_t i = 0; i < parameters.size(); ++i)
    {
        UnlockPar(i);

        double minimum = 0., step;
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
            CalcStep();
            // find the minimum within range
            for(int j = 1; j <= range; ++j)
            {
                step = j*factor;
                // negative direction
                double this_val_m = Evaluate(-step);
                if(this_val_m < eval) {
                    eval = this_val_m;
                    minimum = -step;
                }

                // positive direction
                double this_val_p = Evaluate(step);
                if(this_val_p < eval) {
                    eval = this_val_p;
                    minimum = step;
                }
            }

        } while(minimum != 0.);

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
    // calculate the difference between model and data
    for(size_t i = 0; i < data.size(); ++i)
        p(0, i) = data.at(i).val - GetFormulaVal(data.at(i).x);

    double result = 0;

    // calculate the diff term
    if(weight_diagonal) {
        for(size_t i = 0; i < data.size(); ++i)
            result += p(0, i)*p(0, i)*M_weight(i, i); // faster method for diagonal matrix
    } else {
        result = p*M_weight*transpose(p);
    }

    // penalty matrix exists, calculate penalty term
    if(M_penalty.DimN() == parameters.size() &&
       M_penalty.DimM() == parameters.size())
    {
        // calculate the parameter changes
        CMatrix b(1, parameters.size());
        for(size_t i = 0; i < parameters.size(); ++i)
        {
            b(0, i) = parameters.at(i).initial - parameters.at(i).value;

            // intended change
            if(!parameters.at(i).lock)
                b(0, i) -= parameters.at(i).step*factor;
        }

        // calculate penalty term
        if(penalty_diagonal) {
            for(size_t i = 0; i < parameters.size(); ++i)
                result += b(0, i)*b(0, i)*M_penalty(i, i);
        } else {
            result += b*M_penalty*transpose(b);
        }
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
            formula->SetParameter(i, parameters.at(i).value);
            double gradient = GetFormulaVal(data.at(j).x);
            formula->SetParameter(i, parameters.at(i).value + parameters.at(i).base_step*0.01);
            gradient -= GetFormulaVal(data.at(j).x);
            gradient /= parameters.at(i).base_step*0.01;

            if(gradient == 0.)
                J(j, i) = -1e-10; // put a small number
            else
                J(j, i) = -gradient;
        }
    }

    CMatrix J_w;

    if(weight_diagonal)
        J_w = M_weight.Cholesky_Diagonal()*J; // faster
    else
        J_w = M_weight.Cholesky()*J;

    return J_w.Transpose()*J_w;
}

void CEstimator::CalcStep()
{
    CMatrix rho(parameters.size(), 1);
    for(size_t i = 0; i < rho.DimN(); ++i)
    {
        rho(i, 0) = parameters.at(i).base_step;
    }

    CMatrix step = GetHessian()*rho;

    for(size_t i = 0; i < step.DimN(); ++i) {
        parameters[i].step = 1./step(i, 0);
    }
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

double CEstimator::GetFormulaVal(const double &x)
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
        double expect_val = GetFormulaVal(data.at(i).x);
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
        double expect_val = GetFormulaVal(data.at(i).x);
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
        result += abs(data.at(i).val - GetFormulaVal(data.at(i).x));
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
        double expect_val = GetFormulaVal(data.at(i).x);
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
