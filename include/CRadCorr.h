#ifndef C_RAD_CORR_H
#define C_RAD_CORR_H

#include <cmath>
#include <vector>
#include <string>
#include "ConfigObject.h"
#include "CExpData.h"
#include "CModelWrapper.h"


class CRadCorr : public ConfigObject
{
public:
    CRadCorr();
    virtual ~CRadCorr();

    void Configure(const std::string &path);
    void RadiativeCorrection(CExpData &exp_data, int iters = 0);
    void Radiate(CExpData &exp_data);
    void Initialize(const CExpData &exp_data);
    bool SanityCheck(const CExpData &exp_data);

private:
    void radcor(CExpData::DataSet &dset, bool radiate = false);
    void xyrad2d(CExpData::DataSet &dset, bool radiate = false);
    double fes(const double &Es);
    double fep(const double &Ep);
    double int_es(const double &Es);
    double int_esep(const double &Ep, const double &Es);
    double int_ep(const double &Ep);
    double int_esdp(const double &Es);
    double get_cxsn(const double &E0, const double &Eb);
    void init_model(const CExpData &exp_data, bool radiate = false);
    void find_model_scale(const CExpData::DataSet &mset);
    double from_model(const double &E0, const double &Eb);
    template<typename T>
    void number_operation(const std::string &key, T &val);

    // some lines
    void spectrum_init(CExpData::DataSet &dset);
    void point_init(CExpData::DataPoint &point);
    double __Ep_max(double Es);
    double __Es_min(double Ep);
    double __phi(double x);
    double __eta(double Z);
    double __Q2(double E, double Epr);
    double __log_Q2m2(double E, double Epr);
    double __F_bar(double E, double Epr, double gamma_t);
    double __btr(double E, double Epr);
    double __I(double E0, double E, double xi, double bt);
    double __XI_Stein(double radl);
    double __ice_radl(double thickness);
    double __ice_coll(double thickness);

private:
    bool internal_RC, external_RC, user_defined_XI, peak_approx;
    int n_sim, n_sim_2d;
    double iter_prec, sim_step, sim_step_2d;
    double target_Z, target_A, target_M;
    double theta, sin2, cos2;

    // parameters that will be shared between different functions
    double F_mott, Schwinger, delta, delta1, delta2, Bz; // for whole data sets
    double Es, BTB, BTA, XIB, XIA, GAMT;                 // for each spectrum
    double Ep, R, BTR, Epmin, Epmax, Esmin, Esmax;       // for each data point

    // model related
    double model_scale;
    double model_shift;
};

#endif
