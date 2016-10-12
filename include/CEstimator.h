#ifndef C_ESTIMATOR_H
#define C_ESTIMATOR_H

#include <vector>
#include "CMatrix.h"

class TF1;
class TFormula;

class CEstimator
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
        double prev_value;
        double value;
        double step;
        double base_step;
        double fine_step;
        bool lock;

        Parameter()
        : initial(0.), prev_value(0.), value(0.),
          step(0.), base_step(0.), fine_step(0.), lock(true)
        {};
        Parameter(const double &v, const double &s = 0., const double &fs = 0.)
        : initial(v), prev_value(v), value(v), step(0), lock(true)
        {
            fine_step = (fs < s)? fs : s;

            if(s == 0.)
                base_step = value*0.01;
            if(fs == 0.)
                fine_step = value*0.0001;
        }
    };

public:
    CEstimator();
    virtual ~CEstimator();

    void DeleteFormula();
    void SetFormula(TF1 *tf);
    void SetFormula(const char *c);
    void SetDataPoints(const std::vector<double> &x, const std::vector<double> &y);
    void SetDataPoints(const size_t &n, const double *x, const double *y);
    void SetParameter(const size_t &i,
                      const double &p,
                      const double &step = 0.,
                      const double &fine_step = 0.);
    void SetParameters(const std::vector<double> &p);
    void SetCovarianceMatrix(const CMatrix &m) {V_inv = m.Inverse();};
    void UnlockPar(const size_t &i) {parameters[i].lock = false;};
    void LockPar(const size_t &i) {parameters[i].lock = true;};
    virtual double Evaluate(const double &factor = 0.);
    virtual void NextStep(const double &factor, bool verbose = false);
    virtual CMatrix GetHessian();
    virtual void CalcStep();

    size_t GetNpar() {return parameters.size();};

    Parameter GetParameter(const size_t &i)
    {
        if(i >= parameters.size())
            return Parameter(0);
        return parameters.at(i);
    };

    TFormula *GetFormula() {return formula;};
    std::vector<Parameter> GetParameters() {return parameters;};

private:
    TFormula *formula;
    CMatrix V_inv;
    std::vector<DataPoint> data;
    std::vector<Parameter> parameters;
};

#endif
