#ifndef C_RAD_CORR_H
#define C_RAD_CORR_H

#include <cmath>
#include <vector>
#include <string>

class CRadCorr
{
public:
    CRadCorr();
    virtual ~CRadCorr();

    void ReadExpData(const std::string &path);

public:
    static double gamma(const double &z);
    static double spence(const double &z, const double &res = 1e-15);
    static double spence_tr(const double &z, const double &res, const int &nmax);

    static double simpson(double begin, double end, double (*f)(double), double step, int Nmin)
    {
        int Nsteps = (end - begin)/step;
        int Nbins = std::max(Nmin, Nsteps);
        double s = (end - begin)/(double)(2*Nbins);

        double result = (*f)(begin) + 4*(*f)(begin + s) + (*f)(end);
        double x = begin;
        for(int i = 0; i < Nbins; ++i)
        {
            result += 2*(*f)(x) + 4*(*f)(x + s);
            x += 2*s;
        }

        return result*s/3;
    }

    template<class T>
    static double simpson(double begin, double end, double (T::*f)(double), T *t, double step, int Nmin)
    {
        int Nsteps = (end - begin)/step;
        int Nbins = std::max(Nmin, Nsteps)/2;
        double s = (end - begin)/(double)(2*Nbins);

        double result = (t->*f)(begin) + 4*(t->*f)(begin + s) + (t->*f)(end);
        double x = begin;
        for(int i = 0; i < Nbins; ++i)
        {
            result += 2*(t->*f)(x) + 4.*(t->*f)(x + s);
            x += 2*s;
        }

        return result*s/3;
    }
};


#endif
