#ifndef C_RAD_CORR_H
#define C_RAD_CORR_H

#include <cmath>
#include <vector>
#include <string>
#include "ConfigObject.h"

class CRadCorr : public ConfigObject
{
public:
    struct DataPoint
    {
        // read from data file
        double nu;       // nu (MeV)
        double cxsn;     // cross section
        double stat;     // statistical error
        double syst;     // systematic error

        // calculated
        double Ep;       // final energy before coll. loss
        double PA;       // (Es - Ep)/Es
        double rad;      // Sigma_Radiated
        double born;     // Sigma_Born

        // constructors
        DataPoint() {};
        DataPoint(double n, double c, double st, double sy)
        : nu(n), cxsn(c), stat(st), syst(sy)
        {};
    };

    struct DataSet
    {
        // read from data file
        double energy;        // incident electron energy
        double radl_before;   // radiation length before target
        double radl_after;    // radiation length after target
        double coll_before;   // collision thickness before target
        double coll_after;    // collision thickness after target
        double error;         // relative error for RC
        double normalization; // normalization factor
        bool non_rad;         // non radiated means born cross section from file
        std::vector<DataPoint> data;

        // calculated
        double weight_mott;   // mott loss
        double eta;           // parameter related to target, energy and angle

        // constructors
        DataSet() {};
        DataSet(double e, double er, double no, bool m = true)
        : energy(e), radl_before(0), radl_after(0), coll_before(0), coll_after(0),
          error(er), normalization(no), non_rad(m)
        {};
        DataSet(double e, double rb, double ra, double cb, double ca, double er, double no, bool m = false)
        : energy(e), radl_before(rb), radl_after(ra), coll_before(cb), coll_after(ca),
          error(er), normalization(no), non_rad(m)
        {};
    };

public:
    CRadCorr();
    virtual ~CRadCorr();

    void Configure(const std::string &path);
    void ReadExpData(const char *path);
    void ReadExpData(const std::string &path);
    void ReadExpData(const std::vector<std::string> &filelist);
    bool SanityCheck();
    void RadiativeCorrection(int iters = 10);
    void Radiate();
    void SaveResult(const std::string &path);

private:
    void radcor(DataSet &set, bool radiate = false);
    void xyrad2d(DataSet &set, bool radiate = false);
    double fes(const double &Es);
    double fep(const double &Ep);
    double ftcs(const double &E0, const double &Eb);
    double terp(const DataSet &set, const double &w);
    double int_es(const double &Es);
    double int_ep(const double &Ep);
    double Iprob(const double &E0, const double &E2, const double &t);
    void calculateXI(DataSet &set);
    void readData(ConfigParser &p);

    // some lines
    double __E_max(double E);
    double __E_min(double E);
    double __phi(double x);
    double __Q2(double E, double E_);
    double __log_Q2m2(double E, double E_);
    double __F_bar(double E, double E_, double gamma_t);
    double __btr(double E, double E_);


    std::vector<DataSet> data_sets;
    bool internal_RC, external_RC, user_defined_XI, peak_approx;
    int n_sim;
    double sim_step;
    double target_Z, target_A, target_M;
    double angle, sin2, cos2;

    // parameters that will be shared between different functions
    double Es, Ep, R, delta;
    double F_mott, Schwinger;
    double BTB, BTA, BTR, XIB, XIA, GAMT;

public:
    static double gamma(const double &z);
    static double spence(const double &z, const double &res = 1e-15);
    static double spence_tr(const double &z, const double &res, const int &nmax);

    static double simpson(double begin, double end, double (*f)(const double&), double step, int Nmin)
    {
        int Nsteps = (end - begin)/step;
        int Nbins = std::max(Nmin, Nsteps)/2;
        double s = (end - begin)/(double)(2.*Nbins);

        // first bin
        double result = (*f)(begin) + 4.*(*f)(begin + s) + (*f)(end);
        double x = begin + 2.*s;
        int i = 1;
        while(i++ < Nbins)
        {
            result += 2.*(*f)(x) + 4.*(*f)(x + s);
            x += 2.*s;
        }

        return result*s/3.;
    }

    template<class T>
    static double simpson(double begin, double end, double (T::*f)(const double&), T *t, double step, int Nmin)
    {
        int Nsteps = (end - begin)/step;
        int Nbins = std::max(Nmin, Nsteps)/2;
        double s = (end - begin)/(double)(2.*Nbins);

        // first bin
        double result = (t->*f)(begin) + 4.*(t->*f)(begin + s) + (t->*f)(end);
        double x = begin + 2.*s;
        int i = 1;
        while(i++ < Nbins)
        {
            result += 2.*(t->*f)(x) + 4.*(t->*f)(x + s);
            x += 2.*s;
        }

        return result*s/3.;
    }
};

#endif
