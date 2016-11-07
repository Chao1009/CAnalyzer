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
        double stat;     // statistical error
        double syst;     // systematic error

        // calculated
        double Ep;       // final energy before coll. loss
        double v;        // (Es - Ep)/Es
        double rad;      // Radiated cross section
        double born;     // Born cross section
        double last;     // Save info from last iteration

        // constructors
        DataPoint() {};
        DataPoint(double n, double c, double st, double sy)
        : nu(n), stat(st), syst(sy), rad(c), born(c), last(c)
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
    void RadiativeCorrection(int iters = 0);
    void Radiate();
    void SaveResult(const std::string &path);

private:
    void iterByNumbers(int iters);
    void iterByPrecision();
    void radcor(DataSet &set, bool radiate = false);
    void xyrad2d(DataSet &set, bool radiate = false);
    double fes(const double &Es);
    double fep(const double &Ep);
    double int_es(const double &Es);
    double int_ep(const double &Ep, const double &Es);
    double int_esdp(const double &Es);
    double int_epds(const double &Ep);
    double get_cxsn(const double &E0, const double &Eb);
    double interp(const DataSet &set, const double &w);
    void calculateXI(DataSet &set);
    void readData(ConfigParser &p);

    // some lines
    void spectrum_init(DataSet &set);
    void point_init(DataPoint &point);
    double __Ep_max(double Es);
    double __Es_min(double Ep);
    double __phi(double x);
    double __eta(double Z);
    double __Q2(double E, double Epr);
    double __log_Q2m2(double E, double Epr);
    double __F_bar(double E, double Epr, double gamma_t);
    double __btr(double E, double Epr);
    double __I(double E0, double E, double delta, double bt);


    std::vector<DataSet> data_sets;
    bool internal_RC, external_RC, user_defined_XI, peak_approx;
    int n_sim, n_sim_2d;
    double iter_prec, sim_step, sim_step_2d;
    double target_Z, target_A, target_M;
    double angle, sin2, cos2;

    // parameters that will be shared between different functions
    double F_mott, Schwinger, delta, Bz;           // for whole data sets
    double Es, BTB, BTA, XIB, XIA, GAMT;           // for each spectrum
    double Ep, R, BTR, Epmin, Epmax, Esmin, Esmax; // for each data point

public:
    static double gamma(const double &z);
    static double spence(const double &z, const double &res = 1e-15);
    static double spence_tr(const double &z, const double &res, const int &nmax);

    // simpson integration
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

    template<class T, typename... Args>
    static double simpson(double begin, double end, double step, int Nmin,
                          double (T::*f)(const double&, const Args& ...), T *t, const Args&... args)
    {
        int Nsteps = (end - begin)/step;
        int Nbins = std::max(Nmin, Nsteps)/2;
        double s = (end - begin)/(double)(2.*Nbins);

        // first bin
        double result =  (t->*f)(begin, args...)
                       + 4.*(t->*f)(begin + s, args...)
                       + (t->*f)(end, args...);
        double x = begin + 2.*s;
        int i = 1;
        while(i++ < Nbins)
        {
            result += 2.*(t->*f)(x, args...) + 4.*(t->*f)(x + s, args...);
            x += 2.*s;
        }

        return result*s/3.;
    }
};

#endif
