#include <iostream>
#include "CEstimator.h"
#include "TF1.h"
#include "TFormula.h"

CEstimator::CEstimator()
: formula(nullptr)
{
}

CEstimator::~CEstimator()
{
    DeleteFormula();
}

void CEstimator::DeleteFormula()
{
    if(formula)
        delete formula, formula = nullptr;
    parameters.clear();
}

void CEstimator::SetFormula(const char *c)
{
    DeleteFormula();
    // name, expression, add to global list
    // the formula is taken care by this class
    // thus put false to not add it to root global list
    formula = new TFormula("myform", c, false);
    parameters.resize(formula->GetNpar());
    std::cout << formula->GetExpFormula().Data() << std::endl;
}

void CEstimator::SetFormula(TF1 *tf)
{
    DeleteFormula();

    formula = new TFormula("myform", tf->GetExpFormula().Data(), false);
    for(int i = 0; i < tf->GetNpar(); ++i)
        parameters.emplace_back(tf->GetParameter(i));
}

void CEstimator::SetDataPoints(const std::vector<double> &x, const std::vector<double> &y)
{
    if(x.size() != y.size())
    {
        std::cerr << "CEstimator Error: unmatched vector size from SetDataPoints!"
                  << std::endl;
        return;
    }

    data.clear();

    for(size_t i = 0; i < x.size(); ++i)
    {
        data.emplace_back(x[i], y[i]);
    }
}

void CEstimator::SetDataPoints(const size_t &n, const double *x, const double *y)
{
    data.clear();

    for(size_t i = 0; i < n; ++i)
    {
        data.emplace_back(x[i], y[i]);
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

    CMatrix p_t = p.Transpose();

    CMatrix result = p*V_inv*p_t;

    return result.At(0, 0);
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
            J(j, i) = -gradient;
        }
    }
    CMatrix J_w = V_inv.Cholesky()*J;
    return J_w.Transpose()*J_w;
}

void CEstimator::CalcStep()
{
    CMatrix rho_inv(1, parameters.size());
    for(size_t i = 0; i < rho_inv.DimN(); ++i)
        rho_inv(0, i) = 1/parameters.at(i).base_step;

    CMatrix step = rho_inv*GetHessian();
    for(size_t i = 0; i < step.DimM(); ++i)
        parameters[i].step = 1/step(0, i);
}

void CEstimator::NextStep(const double &factor, bool verbose)
{
    for(size_t i = 0; i < parameters.size(); ++i)
    {
        if(parameters[i].lock)
            continue;
        parameters[i].prev_value = parameters[i].value;
        parameters[i].value += parameters[i].step*factor;
        if(verbose) {
            std::cout << "Parameter " << i << " is changed, "
                      << "previous value: " << parameters[i].prev_value
                      << ", current value: " << parameters[i].value
                      << std::endl;
        }
    }
}
