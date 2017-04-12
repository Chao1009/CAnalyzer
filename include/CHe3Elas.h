#ifndef C_HE3_ELAS_H
#define C_HE3_ELAS_H

#include "ConfigObject.h"
#include "canalib.h"

class CHe3Elas : public ConfigObject
{
public:
    CHe3Elas();
    virtual ~CHe3Elas();

    void Configure(const std::string &path = "");
    void Initialize(bool uxi, double xib, double xia, bool pol, double pol_th);
    // unit: MeV, rad, radiation length, ub/MeV/sr
    double GetRadXS(const double &Es, const double &Ep, const double &theta,
                    const double &radl_in, const double &radl_out);
    static void GetEMFFs(const double &Q2, double &GE, double &GM);
    static double GetBornXS(const double &Es, const double &angle);

private:
    double xyradel(const double &Es, const double &Ep, const double &theta,
                   const double &radl_before, const double &radl_after);
    double int_es(const double &Esx, const double &theta, const double &Es, const double &Ep);
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
    double BTB, BTA, Bz, GAMT, sim_step, xi_before, xi_after, xi_factor;
    int n_sim, pol_angle;
    bool user_xi, polarized, xy_method;
};

#endif

