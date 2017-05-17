// He3 form factors
// Amroun et al. Nucl. Phys. A579, 596-626 (1994).
// Original Author: A. Deur (in fortran)
// Adapted to C++: Chao Peng

#include "CHe3Elas.h"
#include "canalib.h"
#include "rtails.h"
#include <cmath>
#include <iostream>



const double _target_Z = 2.0;
const double _target_A = 3.0149322473;
const double _target_M = _target_A*cana::amu;
const int _n_gaus = 12;
const double _gaus_rms = 0.65;
const double _gaus_mean[] = {0.1, 0.5, 0.9, 1.3,
                             1.6, 2.0, 2.4, 2.9,
                             3.4, 4.0, 4.6, 5.2};
const double _gaus_ampl_c[] = {0.027614, 0.170847, 0.219805, 0.170486,
                               0.134453, 0.100953, 0.074310, 0.053970,
                               0.023689, 0.017502, 0.002034, 0.004338};
const double _gaus_ampl_m[] = {0.059785, 0.138368, 0.281326, 0.000037,
                               0.289808, 0.019056, 0.114825, 0.042296,
                               0.028345, 0.018312, 0.007843, 0.};


// constructor
CHe3Elas::CHe3Elas(const std::string &path)
{
    if(path.size())
        Configure(path);
}

// destructor
CHe3Elas::~CHe3Elas()
{
    // place holder
}

// configuration
void CHe3Elas::Configure(const std::string &path)
{
    if(!path.empty())
        ConfigObject::Configure(path);

    delta1 = getDefConfig<double>("Delta1", 0.5);
    delta2 = getDefConfig<double>("Delta2", 0.5);
    n_sim = getDefConfig<int>("Min Number of Simpson Bins", 10000);
    sim_step = getDefConfig<double>("Simpson Step Size", 0.1);

    polarized = getDefConfig<bool>("Polarized", false);
    pol_angle = getDefConfig<int>("Polarization Angle", 180);

    xi_factor = (cana::pi*cana::ele_mass/2/cana::alpha)
              / (_target_Z + __eta(_target_Z))
              / log(183*std::pow(_target_Z, -1./3.));

    // B(z) Eq. A45 in STEIN
    Bz = 1./9.*(_target_Z + 1.)/(_target_Z + __eta(_target_Z));
    Bz /= log(183.*std::pow(_target_Z, -1./3.));
    Bz = 4./3.*(1. + Bz);


    // TODO need revisit this part, now keep the Stein formula for it
    /*
    user_xi = getDefConfig<bool>("User Defined XI", true);
    xi_before = getDefConfig<double>("XI Before", 0.007417);
    xi_after = getDefConfig<double>("XI After", 0.1069);
    */
    user_xi = false;
    rtails_init(user_xi, polarized, pol_angle);

    std::string app_str = GetConfig<std::string>("Radiative Approach");
    if(app_str.empty() ||
       ConfigParser::case_ins_equal(app_str, "MT Exact")) {
        approach = Mo_Tsai_Exact;
    } else if(ConfigParser::case_ins_equal(app_str, "MT Approx")) {
        approach = Mo_Tsai_Approx;
    } else if(ConfigParser::case_ins_equal(app_str, "BS")) {
        approach = Bardin_Shumeiko;
    }
}

// initialize
void CHe3Elas::Initialize(bool uxi, bool pol, double pol_th)
{
    user_xi = uxi;
    polarized = pol;
    pol_angle = pol_th;

    // pass inputs to initialize rtails
    rtails_init(user_xi, polarized, pol_angle);
}

// input: Four momentum transfer in fm^(-2)
// output: electromagnetic form factors GE and GM
void CHe3Elas::GetEMFFs(const double &Q2, double &GE, double &GM)
{
    double Q = sqrt(Q2);
    double G2 = _gaus_rms*_gaus_rms;

    GE = 0., GM = 0.;
    for(int i = 0; i < _n_gaus; ++i)
    {
        double B = 2.*_gaus_mean[i]*_gaus_mean[i]/G2;

        double QR = Q*_gaus_mean[i];
        if(QR > 0.) {
            double factor = (cos(QR) + B*sin(QR)/QR)/(1. + B);
            GE += factor*_gaus_ampl_c[i];
            GM += factor*_gaus_ampl_m[i];
        } else {
            GE += _gaus_ampl_c[i];
            GM += _gaus_ampl_m[i];
        }
    }

    // Z = 2.0, anomalous magnetic momentum k = -4.2
    GE *= exp(-Q2*G2/4.)*2.0;
    GM *= exp(-Q2*G2/4.)*2.0*(1.0 - 4.2);
}

// input: incident energy in MeV, scattering angle in rad
// output: elastic scattering cross section at Born level in ub/MeV/sr
double CHe3Elas::GetBornXS(const double &Es, const double &theta)
{
    // kinematics calculation
    double sinsq = sin(theta/2.)*sin(theta/2.);
    double recoil = 1./(1. + 2.*Es*sinsq/_target_M);
    double Ep = Es*recoil;
    double Qsq = 4.*Es*Ep*sinsq;
    double tau = Qsq/4./std::pow(_target_M, 2);
    double epsilon = 1./(1. + 2.*(1. + tau)*tan(theta/2.)*tan(theta/2.));
    double xs_mott = cana::alpha*cos(theta/2.)/2./Es/sinsq;
    xs_mott *= xs_mott;

    // get form factors
    double GE, GM;
    GetEMFFs(Qsq/cana::hbarc2, GE, GM);

    // convert unit from MeV to ubarn
    double unit_conv = cana::hbarc2*1e4;

    return xs_mott*recoil/(1. + tau)*(GE*GE + tau/epsilon*GM*GM)*unit_conv;
}

// output the radiated He3 elastic cross section
// input in MeV, rad, unit rad length
// output in ub/MeV/sr
double CHe3Elas::GetRadXS(const double &Es, const double &Ep, const double &theta,
                          const double &radl_in, const double &radl_out,
                          const double &xi_in, const double &xi_out)
{
    // update related variables
    BTB = radl_in*Bz, BTA = radl_out*Bz;
    if(user_xi) {
        xi_before = xi_in;
        xi_after = xi_out;
    } else {
        xi_before = radl_in*xi_factor;
        xi_after = radl_out*xi_factor;
    }
    rtails_set_radl(radl_in, radl_out, xi_before, xi_after);

    switch(approach)
    {
    case Bardin_Shumeiko:
        return rtails_rad_cxsn(Es, Ep, theta, false);
    case Mo_Tsai_Approx:
        return xyradel(Es, Ep, theta);
    case Mo_Tsai_Exact:
    default:
        return rtails_rad_cxsn(Es, Ep, theta, true);
    }
}

// output for xyradel, unit is ub/MeV/sr
// radiation length and collisional loss are determined in GetRadXS
double CHe3Elas::xyradel(const double &Es, const double &Ep, const double &theta)
{
    // update angle dependent terms
    sin2 = std::pow(sin(theta/2.), 2);
    cos2 = 1. - sin2;

    // check allowed kinematics
    double Epmax = __Ep_max(Es), Esmin = __Es_min(Ep);
    if(Ep > Epmax)
        return 0.;

    // Schwinger term will be used in internal radiation
    Schwinger = cana::pi*cana::pi/6 - cana::spence(cos2);

    // update radiation lengths and gamma_t, it is needed by Fbar
    GAMT = cana::gamma(1. + BTB) * cana::gamma(1. + BTA);

    // singularity parts
    double Fbar = __F_bar(Es, Ep, GAMT);
    double BTR = __btr(Es, Ep);

    // Es singular, Es - delta2 to Es
    double sgl_Es = std::pow(delta2/Es, BTB + BTR)/cana::gamma(1. + BTB + BTR)
                  * (1. - xi_before/(1 - BTB - BTR)/delta1)
                  * __I(Epmax, Ep, xi_after, BTA + BTR) * Fbar * GetBornXS(Es, theta);

    double delta1_p = __Ep_max(Esmin + delta1) - __Ep_max(Esmin);
    // Ep singular, Ep to Ep + delta1', transformed from Es_min to Es_min + delta1
    double sgl_Ep = std::pow(delta1_p/Ep, BTA + BTR)/cana::gamma(1. + BTA + BTR)
                  * (1. - xi_after/(1 - BTA - BTR)/delta1_p)
                  * __I(Es, Esmin, xi_before, BTB + BTR)* Fbar * GetBornXS(Esmin, theta)
                  * pow(1./(1. - 2.*Epmax*sin2/_target_M), 2);

    // get integration range
    double Es_beg = Esmin + delta1, Es_end = Es - delta2;
    double int_Es;
    if((Es_beg <= 0) || (Es_end <= 0) || (Es_beg >= Es_end))
    {
        /*
        std::cout << "Warning: skipped point Es = " << Es << ", Ep = " << Ep
                  << ", angle = " << theta * cana::rad2deg
                  << ", since the integration range is not reasonable "
                  << "(Es_min = " << Es_beg << ", Es_max = " << Es_end << ")"
                  << std::endl;
        */
        int_Es = 0.;
    } else {
        int_Es = cana::simpson(Es_beg, Es_end, sim_step, n_sim,
                               &CHe3Elas::int_es, this, theta, Es, Ep);
    }

    return sgl_Es + sgl_Ep + int_Es;
}

// integral over es for elastic process
double CHe3Elas::int_es(const double &Esx, const double &theta,
                          const double &Es, const double &Ep)
{
    // elastic
    double Epx = __Ep_max(Esx);
    // not allowed
    if(Epx <= Ep)
        return 0.;

    double FBAR = __F_bar(Esx, Epx, GAMT);
    double TRx = __btr(Esx, Epx);

    double lost = __I(Es, Esx, xi_before, BTB + TRx)*__I(Epx, Ep, xi_after, BTA + TRx);
    return lost*FBAR*GetBornXS(Esx, theta);
}

// some inlines
// Emax for integration, sin2 and target_M is pre-calculated inside the class
inline double CHe3Elas::__Ep_max(double _Es)
{
    return _Es/(1. + 2.*_Es*sin2/_target_M);
}

// Emin for integration, sin2 and target_M is pre-calculated inside the class
inline double CHe3Elas::__Es_min(double _Ep)
{
    return _Ep/(1. - 2.*_Ep*sin2/_target_M);
}

// Q2 for E and E', sin2 is pre-calculated inside the class
inline double CHe3Elas::__Q2(double _E, double _Epr)
{
    return 4.*_E*_Epr*sin2;
}

// log(Q2/m2) for E and E', used inline __Q2
inline double CHe3Elas::__log_Q2m2(double _E, double _Epr)
{
    return log(__Q2(_E, _Epr)/cana::ele_mass/cana::ele_mass);
}

// phi
inline double CHe3Elas::__phi(double _x)
{
    return 1. - _x + 3.*_x*_x/4.;
}

// eta(Z)
inline double CHe3Elas::__eta(double _Z)
{
    return log(1440.*std::pow(_Z, -2./3.))/log(183.*std::pow(_Z, -1./3.));
}

// Get Fbar(Q2), used inline __log_Q2m2, Schwinger term is pre-calculated
// improvement from J. Singh, higher order terms DHO is exponentiated,
// 1+0.5772*bt term is replaced by two gamma normalization, see GAMT for details
// exp(DHO)/gamma(1+bt) = (1 + 0.5772*bt + ...)*(1 + DHO + ...)
//                      = 1 + 0.5772*bt + DHO + ... [Eq. (2.8) in TSAI71]
inline double CHe3Elas::__F_bar(double _E, double _Epr, double _gamma_t)
{
    double LogQ2m2 = __log_Q2m2(_E, _Epr);
    double Log2EsEp = std::pow(log(_E/_Epr), 2);

    double DHO = 2.*(3./4.*LogQ2m2 - 1.);   // vertex correction
    DHO += 2.*(LogQ2m2/3. - 5./9.);         // vacuum correction
    DHO += Schwinger;                       // Schwinger term, angle dependent
    DHO += -0.5*Log2EsEp;                   // Correction to angle peaking approx.

    DHO *= cana::alpha/cana::pi;                        // common factor

    return exp(DHO)/_gamma_t;
}

// Get tr(Q2), effective radiator thickness before and after the scattering from
// external Bremsstrahlung process, used inline __log_Q2m2
inline double CHe3Elas::__btr(double _E, double _Epr)
{
    return cana::alpha/cana::pi*(__log_Q2m2(_E, _Epr) - 1.);
}

// Probability function I(E0, E, t)
// __XI is accounted for collisional loss 
// NOTICE here we are using b(z)t instead of t
inline double CHe3Elas::__I(double _E0, double _E, double _XI, double _bt)
{
    double _dE = _E0 - _E;
    return std::pow(_dE/_E0, _bt) / cana::gamma(1. + _bt) * (__phi(_dE/_E0)*_bt + _XI/_dE)/_dE;
}
