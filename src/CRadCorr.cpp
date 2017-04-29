#include "CRadCorr.h"
#include "canalib.h"
#include "ConfigParser.h"
#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>

// a hack pointer to make passing parameters for simpson integration easier
const CExpData *interp_source = nullptr;



// constructor
CRadCorr::CRadCorr()
: model_scale(1.0), model_shift(0.0)
{
    // place holder
}

// destructor
CRadCorr::~CRadCorr()
{
    // place holder
}

// configuration
void CRadCorr::Configure(const std::string &path)
{
    // the configuration value is supposed to be provided in config file
    // if the term is not included in configuration file, it will use the
    // default value
    if(!path.empty())
        ConfigObject::Configure(path);

    // get configuration value from the file
    // functions
    internal_RC = getDefConfig<bool>("Internal RC", true);
    external_RC = getDefConfig<bool>("External RC", true);
    user_defined_XI = getDefConfig<bool>("User Defined XI", true);
    peak_approx = getDefConfig<bool>("Energy Peaking Approximation", false);
    use_model = getDefConfig<bool>("Extrapolation By Model", false);

    // calculation related
    iter_prec = getDefConfig<double>("Iteration Precision", 0.005);
    n_sim = getDefConfig<int>("Min Number of Simpson Bins", 10000);
    sim_step = getDefConfig<double>("Simpson Step Size", 0.1);
    // for non-peaking-approximation only
    // 2d integral is much slower than 1d, thus configure the step carefully
    n_sim_2d = getDefConfig<int>("Min 2D Simpson Bins", 100);
    sim_step_2d = getDefConfig<double>("2D Simpson Step Size", 1.);

    // data related
    delta = getDefConfig<double>("IR DIV Delta", 10);
    // for non-peaking-approximation only
    delta1 = getDefConfig<double>("Delta1", 5);
    delta2 = getDefConfig<double>("Delta2", 5);
}


// pre-calculate some terms that won't change during radiation calculation
void CRadCorr::Initialize(const CExpData &exp_data, bool radiate)
{
    // target
    target_Z = exp_data.TargetZ();
    target_A = exp_data.TargetA();
    // scattering angle
    theta = exp_data.Angle()*cana::deg2rad;

    // calculate values based on the input configuration
    // common value for all spectrums
    target_M = target_A * cana::amu;
    sin2 = std::pow(sin(theta/2.), 2);
    cos2 = std::pow(cos(theta/2.), 2);

    // Schwinger term in internal radiation
    Schwinger = cana::pi*cana::pi/6 - cana::spence(cos2);

    // B(z) Eq. A45 in STEIN
    Bz = 1./9.*(target_Z + 1.)/(target_Z + __eta(target_Z));
    Bz /= log(183.*std::pow(target_Z, -1./3.));
    Bz = 4./3.*(1. + Bz);

    if(use_model) {
        // force using model for initializing
        interp_source = nullptr;
        init_model(exp_data, radiate);
    }

    interp_source = &exp_data;
}

// test some configuration values, warn some simple mistake
bool CRadCorr::SanityCheck(const CExpData &exp_data)
{
    if(!internal_RC && !external_RC) {
        std::cout << "Both internal and external RC are OFF, no need to run."
                  << std::endl;
        return false;
    }

    if(exp_data.Empty()) {
        std::cout << "There are no data."
                  << std::endl;
        return false;
    }

    if(sim_step <= 0 || sim_step_2d <= 0) {
        std::cout << "Simpson integration step size must be > 0"
                  << std::endl;
        return false;
    }

    if(delta <= 0 || delta1 <= 0 || delta2 <= 0) {
        std::cout << "Delta must be > 0 MeV"
                  << std::endl;
        return false;
    }

    return true;
}

// do radiative correction on data sets
// if iteration number is provided and > 0, it will do iterations by the input
// number
// if else, it will do iterations until the result diff. reached the iter_prec
void CRadCorr::RadiativeCorrection(CExpData &exp_data, int iters)
{
    Initialize(exp_data, false);

    if(!SanityCheck(exp_data))
        return;

    bool end_by_prec = (iters <= 0)? true : false;

    for(int iter = 1; (iter <= iters) || end_by_prec; ++iter)
    {
        // do radiative correction for all data sets
        for(auto &dset : exp_data.GetSets())
        {
            // skip Born Level data/model
            if(dset.non_rad)
                continue;

            std::cout << "Iteration " << iter << ", spectrum E = " << dset.energy
                      << std::endl;

            // energy peaking or not
            if(peak_approx)
                radcor(dset);
            else
                xyrad2d(dset);
        }

        // after correction, check if this iteration reaches the precision
        if(end_by_prec) {
            // get maximum difference between current value and last iteration value
            double max_rel_diff = 0.;
            for(auto &dset : exp_data.GetSets())
            {
                for(auto &point : dset.data)
                {
                    double rel_diff = fabs(1. - point.last/point.born);
                    if(rel_diff > max_rel_diff)
                        max_rel_diff = rel_diff;
                }
            }

            // end iteration message
            std::cout << "Iteration " << iter << " done."
                      << "Maximum change = " << max_rel_diff*100. << "%, "
                      << "required precision = " << iter_prec*100. << "%. "
                      << std::endl;

            // need continue
            if(max_rel_diff > iter_prec) {
                std::cout << "Not converging within the required precision, "
                          << "continue iterations..."
                          << std::endl;
            // Done the radiative correction
            } else {
                std::cout << "Converged within the required precision, done!"
                          << std::endl;
                return;
            }
        }

    }
}

// radiate data sets
void CRadCorr::Radiate(CExpData &exp_data)
{
    Initialize(exp_data, true);

    if(!SanityCheck(exp_data))
        return;

    for(auto &s : exp_data.GetSets())
    {
        // only radiate for Born Level
        if(!s.non_rad)
            continue;

        std::cout << "Radiate, spectrum energy: " << s.energy << std::endl;

        if(peak_approx)
            radcor(s, true);
        else
            xyrad2d(s, true);
    }
}

//==============================================================================
// radiative correction for one spectrum                                        
//========================ORIGINAL AUTHORS======================================
// RADCOR FORTRAN CODE                                                          
// ***by Randy Roy Whitney, Phys. Rev. C 9, 2230 - 2235 (1974)                  
//------------------------------------------------------------------------------
// MODIFIED RADCOR - 01/29/2003                                                 
// ***by K. Slifer, Temple University                                           
// downloadable at  http://www.jlab.org/~slifer/codes.html                      
//------------------------------------------------------------------------------
// MODIFIED RADCOR - 04/01/2007                                                 
// ***by Jaideep Singh                                                          
//========================LATEST CHANGES========================================
// 1. Adapted the code in C++                                                   
// 2. Corrected the calculation in cross section interpolation                  
// ***by Chao Peng, Duke University, 11/8/2016                                  
//========================REFERENCES============================================
// MOTS69    Radiative Corrections to Elastic and Inelastic ep and mu-p         
//           Scattering                                                         
//           by: L.W. Mo and Y.S. Tsai                                          
//           Reviews of Modern Physics, 41, p205-235, Jan 1969                  
//           http://link.aps.org/abstract/RMP/v41/p205                          
//------------------------------------------------------------------------------
// TSAI71    Radiative Corrections to Electron Scattering                       
//           by: Yung-Su Tsai, SLAC PUB 848, Jan 1971                           
//           http://www.slac.stanford.edu/pubs/slacpubs/0750/slac-pub-0848.pdf  
//------------------------------------------------------------------------------
// STEIN     Electron scattering at 4° with energies of 4.5-20 GeV              
//           by: S. Stein, W. B. Atwood, E. D. Bloom, R. L. A. Cottrell,        
//               H. DeStaebler, C. L. Jordan§, H. G. Piel, C. Y. Prescott,      
//               R. Siemann, and R. E. Taylor                                   
//           Phys. Rev. D 12, 1884 - 1919 (October 1975)                        
//           http://link.aps.org/abstract/PRD/v12/p1884                         
//------------------------------------------------------------------------------
// Miller72  Inelastic Electron-Proton Scattering at Large Momentum Transfers   
//           and the Inelastic Structure Functions of the Proton                
//           by: G. Miller, E. D. Bloom, G. Buschhorn,  D. H. Coward,           
//               H. DeStaebler, J. Drees, C. L. Jordan, L. W. Mo, R. E. Taylor, 
//               J. I. Friedman, G. C. Hartmanna, H. W. Kendall, and R. Verdier 
//           Phys. Rev. D 5, 528 - 544 (Feb. 1972)                              
//           http://link.aps.org/abstract/PRD/v5/p528                           
//==============================================================================
void CRadCorr::radcor(CExpData::DataSet &s, bool radiate, bool verbose)
{
    spectrum_init(s);

    int count = 0;
    // iteration on data points
    for(auto &point : s.data)
    {
        if(verbose) {
            std::cout << "Data points: " << count++ << "/" << s.data.size() << "\r"
                      << std::flush;
        }

        point_init(point);

        // components of SIGRAD
        double SIGLOW = 0, SIGBEF = 0, SIGAFT = 0;

        double FBAR = __F_bar(Es, Ep, GAMT);
        SIGLOW = FBAR * std::pow(R*delta/Es, BTB + BTR)
                      * std::pow(delta/Ep, BTA + BTR)
                      * (1. - (XIB+XIA)/delta/(1. - BTB - BTA - 2.*BTR));

        // calculate integral along dEs for fixed Ep SIGBEF
        if((Esmin <= 0) || (Esmax <= 0) || (Esmin >= Esmax)) {
            if(verbose) {
                std::cout << "Skip point at Ep = " << Ep
                          << ", in spectrum E = " << Es
                          << ", ESMIN = " << Esmin
                          << ", ESMAX = " << Esmax
                          << std::endl;
            }
        } else {
            SIGBEF = cana::simpson(Esmin, Esmax, &CRadCorr::fes, this, sim_step, n_sim);
        }

        // calculate integral along dEp for fixed Es SIGAFT
        if((Epmin <= 0) || (Epmax <= 0) || (Epmin >= Epmax)) {
            if(verbose) {
                std::cout << "Skip point at Ep = " << Ep
                          << ", in spectrum E = " << Es
                          << ", EPMIN = " << Esmin
                          << ", EPMAX = " << Esmax
                          << std::endl;
            }
        } else {
            SIGAFT = cana::simpson(Epmin, Epmax, &CRadCorr::fep, this, sim_step, n_sim);
        }

        if(radiate) {
            // radiate, update radiated cross section
            point.rad = SIGLOW*point.born + (SIGBEF + SIGAFT);
        } else {
            // save last iteration
            point.last = point.born;
            // radiative correction, update the born cross section
            point.born = (point.rad - (SIGBEF+SIGAFT))/SIGLOW;
        }
    }

    if(verbose) {
        std::cout << "Data points: " << count << "/" << s.data.size() << ", finished!"
                  << std::endl;
    }
}

// for integral along dEs
double CRadCorr::fes(const double &Esx)
{
    double FBAR = __F_bar(Esx, Ep, GAMT);
    double TRx = __btr(Esx, Ep);

    // calculate effect of multiple soft photon emission
    double FSOFT = std::pow(1. - Esx/Es, BTB+TRx)*std::pow((Es-Esx)/Ep/R, BTA+TRx);

    double FES = (BTB + TRx)/(Es - Esx)*__phi(1. - Esx/Es);
    FES += XIB/(Es - Esx)/(Es - Esx);
    FES *= FSOFT*FBAR*get_cxsn(Esx, Ep);
    return FES;
}

// for integral along dEp
double CRadCorr::fep(const double &Epx)
{
    double FBAR = __F_bar(Es, Epx, GAMT);
    double TRx = __btr(Es, Epx);

    // calculate effect of multiple soft photon emission
    double FSOFT = std::pow(1. - Ep/Epx, BTA+TRx)*std::pow((Epx-Ep)*R/Es, BTB+TRx);

    double FEP = (BTA + TRx)/(Epx - Ep)*__phi(1. - Ep/Epx);
    FEP += XIA/(Epx - Ep)/(Epx - Ep);
    FEP *= FSOFT*FBAR*get_cxsn(Es, Epx);
    return FEP;
}

//==============================================================================
// radiative correction for one spectrum without peaking approximation          
//========================ORIGINAL AUTHORS======================================
// XYRAD2D FORTRAN CODE                                                         
// ***by Xuefei Yan, Duke University,  10/11/2016                               
//========================LATEST CHANGES========================================
// 1. Adapted the code in C++                                                   
// 2. Now get cross section from interpolation of input data instead of model   
// 3. Add external radiative correction part                                    
// 4. Using exponentialted higher order term and gamma term in Fbar calculation 
// ***by Chao Peng, Duke University, 11/8/2016                                  
//========================REFERENCES============================================
// MOTS69    Radiative Corrections to Elastic and Inelastic ep and mu-p         
//           Scattering                                                         
//           by: L.W. Mo and Y.S. Tsai                                          
//           Reviews of Modern Physics, 41, p205-235, Jan 1969                  
//           http://link.aps.org/abstract/RMP/v41/p205                          
//------------------------------------------------------------------------------
// TSAI71    Radiative Corrections to Electron Scattering                       
//           by: Yung-Su Tsai, SLAC PUB 848, Jan 1971                           
//           http://www.slac.stanford.edu/pubs/slacpubs/0750/slac-pub-0848.pdf  
//==============================================================================
void CRadCorr::xyrad2d(CExpData::DataSet &s, bool radiate, bool verbose)
{
    spectrum_init(s);

    int count = 0;
    // iteration on data points
    for(auto &point : s.data)
    {
        if(verbose) {
            std::cout << "Data points: " << count++ << "/" << s.data.size() << "\r"
                      << std::flush;
        }

        point_init(point);

        // components of SIGRAD
        double int_2d = 0, sgl_Es = 0, sgl_Ep = 0, sgl_both = 0;
        double FBAR = __F_bar(Es, Ep, GAMT);

                        // Es singularity
        sgl_both = FBAR*std::pow(delta1/Es, BTB + BTR)/cana::gamma(1. + BTB + BTR)
                        // coll. loss term
                       *(1. - XIB/(1. - BTB - BTR)/delta1)
                        // Ep singularity
                       *std::pow(delta2/Ep, BTA + BTR)/cana::gamma(1. + BTA + BTR)
                        // coll. loss term
                       *(1. - XIA/(1. - BTA - BTR)/delta2);

        if((Esmin <= 0) || (Esmax <= 0) || (Esmin >= Esmax)) {
            if(verbose) {
                std::cout << "Skip point at Ep = " << Ep
                          << ", in spectrum E = " << Es
                          << ", ESMIN = " << Esmin
                          << ", ESMAX = " << Esmax
                          << std::endl;
            }
        } else {
            int_2d = cana::simpson(Esmin, Esmax, &CRadCorr::int_es, this, sim_step_2d, n_sim_2d);
            sgl_Ep = cana::simpson(Esmin, Esmax, &CRadCorr::int_esdp, this, sim_step, n_sim);
        }

        if((Epmin <= 0) || (Epmax <= 0) || (Epmin >= Epmax)) {
            if(verbose) {
                std::cout << "Skip point at Ep = " << Ep
                          << ", in spectrum E = " << Es
                          << ", EPMIN = " << Esmin
                          << ", EPMAX = " << Esmax
                          << std::endl;
            }
        } else {
            sgl_Es = cana::simpson(Epmin, Epmax, &CRadCorr::int_ep, this, sim_step, n_sim);
            // Es singularity
            sgl_Es *= std::pow(delta1/Es, BTB + BTR)/cana::gamma(1. + BTB + BTR);
            // coll. loss term
            sgl_Es *= 1. - XIB/(1 - BTB - BTR)/delta1;
        }

        if(radiate) {
            // radiate, update radiated cross section
            point.rad = sgl_both*point.born + (sgl_Es + sgl_Ep + int_2d);
        } else {
            // save last iteration
            point.last = point.born;
            // radiative correction, update the born cross section
            point.born = (point.rad - (sgl_Es + sgl_Ep + int_2d))/sgl_both;
        }
    }

    if(verbose) {
        std::cout << "Data points: " << count << "/" << s.data.size() << ", finished!"
                  << std::endl;
    }

}

// it is a 2d integral int_es(int_ep)
double CRadCorr::int_es(const double &Esx)
{
    double Ep_max = __Ep_max(Esx);
    double Ep_min = delta2 + Ep; // delta2
    if(Ep_max < Ep_min)
        return 0.;

    return cana::simpson(Ep_min, Ep_max, sim_step_2d, n_sim_2d, &CRadCorr::int_esep, this, Esx);
}

// should be inside int_es
double CRadCorr::int_esep(const double &Epx, const double &Esx)
{
    double FBAR = __F_bar(Esx, Epx, GAMT);
    double TRx = __btr(Esx, Epx);
    return FBAR*get_cxsn(Esx, Epx)*__I(Epx, Ep, XIA, BTA + TRx)*__I(Es, Esx, XIB, XIB + TRx);
}

// integral over ep and es singularity
double CRadCorr::int_ep(const double &Epx)
{
    double FBAR = __F_bar(Es, Epx, GAMT);
    double TRx = __btr(Es, Epx);

    double lost = __I(Epx, Ep, XIA, BTA + TRx);
    return lost*FBAR*get_cxsn(Es, Epx);
}

// integral over es and ep singularity
double CRadCorr::int_esdp(const double &Esx)
{
    double FBAR = __F_bar(Esx, Ep, GAMT);
    double TRx = __btr(Esx, Ep);

    double lost = __I(Es, Esx, XIB, BTB + TRx);
    double int_dp = std::pow(delta2/Esx, BTA + TRx)/cana::gamma(1. + BTA + TRx);
    // coll. loss term
    int_dp *= 1. - XIA/(1 - BTA - TRx)/delta2;

    return lost*int_dp*FBAR*get_cxsn(Esx, Ep);
}

// interpolates or extrapolates
inline double CRadCorr::get_cxsn(const double &E0, const double &Eb)
{
    // not allowed
    if(Eb > __Ep_max(E0))
        return 0;

    // use model for extrapolation
    if(!interp_source || (use_model && !interp_source->InRange(E0))) {
        return from_model(E0, Eb);
    // interpolation from model
    } else {
        return interp_source->GetCrossSection(E0, Eb);
    }
}

// a helper function to find the peak point
typedef std::vector<CExpData::DataPoint>::const_iterator CIter;
CIter find_peak(const std::vector<CExpData::DataPoint> &data, bool born_level)
{
    // find the peak point of the first data set
    CIter res;
    double peak_xs = 0.;

    for(auto it = data.cbegin(); it != data.cend(); ++it)
    {
        if(!born_level && it->rad > peak_xs) {
            peak_xs = it->rad;
            res = it;
        }

        if(born_level && it->born > peak_xs) {
            peak_xs = it->born;
            res = it;
        }
    }

    return res;
}

// initialize the model
// determine the scale and shift from the data at lowest energy
void CRadCorr::init_model(const CExpData &exp_data, bool born_level)
{
    model_scale = 1.0, model_shift = 0.0;

    // copy the first data set
    CExpData::DataSet model_data;
    for(auto &dset : exp_data.GetSets())
    {
        if(!born_level && !dset.non_rad) {
            model_data = dset;
            break;
        }

        if(born_level && dset.non_rad) {
            model_data = dset;
            break;
        }
    }

    // cannot find the data set or there are no data points
    if(model_data.data.empty())
        return;

    std::cout << "Initializing model scaling and shifting factors from data "
              << "E = " << model_data.energy
              << std::endl;

    auto max_it = find_peak(model_data.data, born_level);
    size_t max_idx = max_it - model_data.data.begin();
    CExpData::DataPoint max_point(*max_it);

    // keep the i +- 20 points around the data set
    std::vector<CExpData::DataPoint> model_points;
    size_t beg = (max_idx > 20) ? (max_idx - 20) : 0;
    size_t end = max_idx + 20;
    for(size_t i = beg; i < end && i < model_data.data.size(); ++i)
    {
        CExpData::DataPoint point = model_data.data.at(i);
        point.born = from_model(model_data.energy, point.Ep);
        model_points.push_back(point);
    }
    model_data.data = model_points;

    // if non-radiated, we can get the model scale factor now
    if(born_level) {
        CExpData::DataPoint max_model(*find_peak(model_data.data, born_level));
        model_scale = max_point.born/max_model.born;

    // radiated will be a little bit more complicated
    } else {
        // iteratively determine the model scale factor
        int model_iter = 0;
        while(model_iter++ < 20)
        {
            if(peak_approx)
                radcor(model_data, true, false);
            else
                xyrad2d(model_data, true, false);

            CExpData::DataPoint max_model(*find_peak(model_data.data, born_level));

            model_scale *= max_point.rad/max_model.rad;

            // good agreement, no need to continue
            if(std::abs(1.0 - max_point.rad/max_model.rad) < 0.001)
                break;

            // update for next iteration
            for(auto &point : model_data.data)
            {
                point.born = from_model(model_data.energy, point.Ep);
            }
        }
    }

    std::cout << "Initialized, model scale = " << model_scale
              << ", model shift = " << model_shift
              << std::endl;
}

// get cross sections from model
inline double CRadCorr::from_model(const double &E0, const double &Eb)
{
    // bosted model
    double cxsn;
    QFS_xs(target_Z, target_A, E0, Eb - model_shift, theta, &cxsn);
    //Bosted_xs(target_Z, target_A, E0, Eb - model_shift, angle, &cxsn);
    return cxsn*model_scale;
}

// spectrum based kinematics intialization
inline void CRadCorr::spectrum_init(CExpData::DataSet &s)
{
    // update these parameters when external RC is ON
    if(external_RC) {
        // LATEST CHANGE: replaced b(z) = 4./3. with Eq. A45 from STEIN
        // originally it was an approximated value regardless of Z dependence,
        // ignoring it introducing an error on a few percent level.
        BTB = s.radl_before*Bz;
        BTA = s.radl_after*Bz;
        if(user_defined_XI) {
            XIB = s.coll_before;
            XIA = s.coll_after;
        } else {
            XIB = __XI_Stein(s.radl_before);
            XIA = __XI_Stein(s.radl_after);
        }
    } else {
        BTB = 0;
        BTA = 0;
        XIB = 0;
        XIA = 0;
    }

    Es = s.energy;
    // improvements by J. Singh
    // two cana::gamma function normalization so
    // 1/cana::gamma(1 + b*tb)/gamma(1 + b*ta) = 1 + 0.5772b*(tb + ta) + ...
    GAMT = cana::gamma(1. + BTB) * cana::gamma(1. + BTA);
}

// data point based kinematics initialization
inline void CRadCorr::point_init(CExpData::DataPoint &point)
{
    // kinematics
    Ep = point.Ep;

    Esmin = __Es_min(Ep);
    Epmax = __Ep_max(Es);

    if(peak_approx) {
        R = Esmin/Epmax * Es/Ep;
        Esmax = Es - R*delta;
        Epmin = Ep + delta;
    } else {
        Esmax = Es - delta1;
        Epmin = Ep + delta2;
    }

    // equivalent radiator from Bremsstrahlung
    BTR = __btr(Es, Ep);
}

// some inline functions
//
// Emax for integration, sin2 and target_M is pre-calculated inside the class
inline double CRadCorr::__Ep_max(double _Es)
{
    return _Es/(1. + 2.*_Es*sin2/target_M);
}

// Emin for integration, sin2 and target_M is pre-calculated inside the class
inline double CRadCorr::__Es_min(double _Ep)
{
    return _Ep/(1. - 2.*_Ep*sin2/target_M);
}

// Q2 for E and E', sin2 is pre-calculated inside the class
inline double CRadCorr::__Q2(double _E, double _Epr)
{
    return 4.*_E*_Epr*sin2;
}

// log(Q2/m2) for E and E', used inline __Q2
inline double CRadCorr::__log_Q2m2(double _E, double _Epr)
{
    return log(__Q2(_E, _Epr)/cana::ele_mass/cana::ele_mass);
}

// phi
inline double CRadCorr::__phi(double _x)
{
    return 1. - _x + 3.*_x*_x/4.;
}

// eta(Z)
inline double CRadCorr::__eta(double _Z)
{
    return log(1440.*std::pow(_Z, -2./3.))/log(183.*std::pow(_Z, -1./3.));
}

// Get Fbar(Q2), used inline __log_Q2m2, Schwinger term is pre-calculated
// improvement from J. Singh, higher order terms DHO is exponentiated,
// 1+0.5772*bt term is replaced by two gamma normalization, see GAMT for details
// exp(DHO)/gamma(1+bt) = (1 + 0.5772*bt + ...)*(1 + DHO + ...)
//                      = 1 + 0.5772*bt + DHO + ... [Eq. (2.8) in TSAI71]
inline double CRadCorr::__F_bar(double _E, double _Epr, double _gamma_t)
{
    if(!internal_RC)
        return 1./_gamma_t;

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
inline double CRadCorr::__btr(double _E, double _Epr)
{
    if(!internal_RC)
        return 0.;

    return cana::alpha/cana::pi*(__log_Q2m2(_E, _Epr) - 1.);
}

// Probability function I(E0, E, t)
// __XI is accounted for collisional loss 
// NOTICE here we are using b(z)t instead of t
inline double CRadCorr::__I(double _E0, double _E, double _XI, double _bt)
{
    double _dE = _E0 - _E;
    return std::pow(_dE/_E0, _bt) / cana::gamma(1. + _bt) * (__phi(_dE/_E0)*_bt + _XI/_dE)/_dE;
}

// calculate XI based on Stein's formula, require radiation length as input
inline double CRadCorr::__XI_Stein(double radl)
{
    // formula from Stein, but old and probably wrong estimation
    double xi = (cana::pi*cana::ele_mass/2/cana::alpha);
    xi /= (target_Z + __eta(target_Z));
    xi /= log(183*std::pow(target_Z, -1./3.));

    return xi*radl;
}

template<typename T>
void CRadCorr::number_operation(const std::string &key, T &val)
{
    auto operations = ConfigParser::split(GetConfig<std::string>(key), ",");

    while(operations.size())
    {
        std::string op = ConfigParser::trim(operations.front(), " \t");
        operations.pop_front();
        try {
            double number = stod(op.substr(1));
            if(op.at(0) == '+') {
                val += number;
            } else if(op.at(0) == '-') {
                val -= number;
            } else if(op.at(0) == '*') {
                val *= number;
            } else if(op.at(0) == '/') {
                val /= number;
            } else {
                std::cout << "Does not support operator type " << op.at(0)
                          << ", skip operation of " << op
                          << std::endl;
            }
        } catch(std::exception &e) {
            std::cout << "Error: " << e.what()
                      << ", skip operation of " << op
                      << std::endl;
        }
    }
}
