#ifndef C_HE3_ELAS_H
#define C_HE3_ELAS_H

#include "ConfigObject.h"
#include "canalib.h"

class CHe3Elas : public ConfigObject
{
public:
    enum RadApproach
    {
        Mo_Tsai_Approx,
        Mo_Tsai_Exact,
        Bardin_Shumeiko,
    };

    CHe3Elas(const std::string &path = "");
    virtual ~CHe3Elas();

    void Configure(const std::string &path = "");
    void Initialize(bool uxi, bool pol, double pol_th);
    // unit: MeV, rad, radiation length, ub/MeV/sr
    double GetRadXS(double Es, double Ep, double theta,
                    double radl_in, double radl_out,
                    double xi_in = 0., double xi_out = 0.);
    static void GetEMFFs(double Q2, double &GE, double &GM);
    static double GetBornXS(double Es, double theta);

private:
    double xyradel(double Es, double Ep, double theta);
    double int_es(double Esx, double theta, double Es, double Ep);
    double __Ep_max(double _Es);
    double __Es_min(double _Ep);
    double __Q2(double _E, double _Epr);
    double __log_Q2m2(double _E, double _Epr);
    double __phi(double _x);
    double __eta(double _Z);
    double __F_bar(double _E, double _Epr, double _gamma_t);
    double __btr(double _E, double _Epr);
    double __I(double _E0, double _E, double _XI, double _bt);

private:
    double delta1, delta2, sin2, cos2, Schwinger;
    double BTB, BTA, Bz, GAMT, sim_step, xi_factor;
    double xi_before, xi_after;
    int n_sim, pol_angle;
    bool user_xi, polarized, xy_method;
    RadApproach approach;
};

#endif

