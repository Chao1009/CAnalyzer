#ifndef C_ESTIMATOR_H
#define C_ESTIMATOR_H

#include <vector>
#include "CMatrix.h"
class TFormula

class CEstimater
{
public:
    struct DataPoint
    {
        double x;
        double val;

        DataPoint()
        : x(0.), val(0.)
        {};

        DataPoint(const double &xi, const double &vi)
        : x(xi), val(vi)
        {};
    };

    struct Parameter
    {
        double initial;
        double value;
        double step;
        double fine_step;
        bool lock;

        Parameter(const double &v, const double &s = 0., const double &fs = 0.)
        : initial(v), value(v), lock(false)
        {
            fine_step = (fs < s)? fs : s;

            if(s == 0.)
                step = value*0.01;
            if(fs == 0.)
                fine_step = value*0.0001;
        }
    };

public:
    CEstimator();
    virtual ~CEstimator();

    void DeleteFormula();
    void SetFormula(const char *c);
    void SetDataPoints(const std::vector<double> &x, const std::vector<double> &y);
    void SetDataPoints(const size_t &n, const double *x, const double *y);
    void SetParameter(const size_t &i,
                      const double &p,
                      const double &step = 0.,
                      const double &fine_step = 0.);
    void SetParameters(const std::vector<double> &p,
                       const std::vector<double> &step = std::vector<double>(),
                       const std::vector<double> &fine_step = std::vector<double>());
    void SetCovarianceMatrix(const CMatrix &m);
    virtual double Evaluate();

private:
    TFormula *formula;
    CMatrix Covariance;
    std::vector<DataPoint> data;
    std::vector<Parameter> parameters;
};

#endif
