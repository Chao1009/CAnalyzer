#include <iostream>
#include "CEstimator.h"
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

void SetFormula(const char *c);
{
    DeleteFormula();
    // name, expression, add to global list
    // the formula is taken care by this class
    // thus put false to not add it to root global list
    formula = new TFormula("myform", c, false);
    for(int i = 0; i < formula->GetNPar(); ++i)
    {
        parameters.emplace_back((double)formula->GetParameter(i));
    }
}

void SetDataPoints(const std::vector<double> &x, const std::vector<double> &y
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

void SetDataPoints(size_t &n, const double *x, const double *y)
{
    data.clear();

    for(size_t i = 0; i < n; ++i)
    {
        data.emplace_back(x[i], y[i]);
    }
}

void SetParameter(const size_t &i, const double &p, const double &step, const double &fine_step)
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

    parameters.emplace_back(p, step, fine_step);
}
