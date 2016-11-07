#include "CRadCorr.h"
#include "ConfigParser.h"
#include <iostream>
#include <iomanip>
#include <algorithm>


#define ALPHA 7.297352568E-3    // 1./137.03599911
#define PI 3.1415926535897932   // pi
#define RADDEG 57.2957795131    // rad to degree
#define ELECM 0.510998918       // MeV
#define HBARC 197.326968        // hbar*c (MeV*fm)
#define AMUMEV 931.494043       // MeV per amu


// constructor
CRadCorr::CRadCorr()
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
    readConfigFile(path);

    internal_RC = getConfig<bool>("Internal RC", true);
    external_RC = getConfig<bool>("External RC", true);
    user_defined_XI = getConfig<bool>("User Defined XI", false);
    peak_approx = getConfig<bool>("Peaking Approximation", true);

    iter_prec = getConfig<double>("Iteration Precision", 0.005);
    n_sim = getConfig<int>("Min Number of Simpson Bins", 10000);
    sim_step = getConfig<double>("Simpson Step Size", 0.1);
    n_sim_2d = getConfig<int>("Min 2D Simpson Bins", 100);
    sim_step_2d = getConfig<double>("2D Simpson Step Size", 1.);

    delta = getConfig<double>("IR DIV Delta", 10);
    delta1 = getConfig<double>("Delta1", 5);
    delta2 = getConfig<double>("Delta2", 5);
    target_Z = getConfig<double>("Target Z", 2);
    target_A = getConfig<double>("Target A", 3.0149322473);
    angle = getConfig<double>("Scattering Angle", 6.1);

    // common value for all spectrums
    angle /= RADDEG;
    target_M = target_A * AMUMEV;
    sin2 = sin(angle/2)*sin(angle/2);
    cos2 = cos(angle/2)*cos(angle/2);

    // mott cross section
    F_mott = HBARC*ALPHA*cos(angle/2)/2/sin2;
    F_mott = F_mott*F_mott*1e7; // MeV^2*nb/sr

    // Schwinger term in internal radiation
    Schwinger = PI*PI/6 - spence(cos2);

    // B(z) Eq. A45 in STEIN
    Bz = 1./9.*(target_Z + 1.)/(target_Z + __eta(target_Z));
    Bz /= log(183.*std::pow(target_Z, -1./3.));
    Bz = 4./3.*(1. + Bz);

}

// read one data file, see comments in readData() for format details
void CRadCorr::ReadExpData(const char *path)
{
    std::vector<std::string> list;
    list.emplace_back(path);
    ReadExpData(list);
}

void CRadCorr::ReadExpData(const std::string &path)
{
    std::vector<std::string> list;
    list.push_back(path);
    ReadExpData(list);
}

// read several data files
void CRadCorr::ReadExpData(const std::vector<std::string> &file_list)
{
    ConfigParser c_parser;

    // clear data sets first
    data_sets.clear();

    for(auto &path: file_list)
    {
        if(!c_parser.OpenFile(path)) {
            std::cout << "Cannot Open Data File "
                      << "\"" << path << "\""
                      << std::endl;
            continue;
        }

        readData(c_parser);

        c_parser.CloseFile();
    }

    // sort data sets by energy
    std::sort(data_sets.begin(), data_sets.end(),
              [] (const DataSet &set1, const DataSet &set2)
              {
                  return set1.energy < set2.energy;
              });

    // calculate collision thickness if not using users' input
    if(!user_defined_XI) {
        for(auto &set : data_sets)
            calculateXI(set);
    }

    std::cout << "Read " << data_sets.size() << " data sets."
              << std::endl;

    for(auto &set : data_sets)
    {
        std::cout << "Energy: " << std::setw(8) << set.energy
                  << ",    DataPoints:" << std::setw(8) << set.data.size()
                  << std::endl;
    }

}

bool CRadCorr::SanityCheck()
{
    if(!internal_RC && !external_RC) {
        std::cout << "Both internal and external RC are OFF, no need to run."
                  << std::endl;
        return false;
    }

    if(data_sets.empty()) {
        std::cout << "There are no data."
                  << std::endl;
        return false;
    }

    bool radiated_data = false;
    for(auto &set : data_sets)
    {
        if(!set.non_rad) {
            radiated_data = true;
            break;
        }
    }

    if(!radiated_data) {
        std::cout << "There are no radiated data."
                  << std::endl;
        return false;
    }

    if(sim_step < 0) {
        std::cout << "Simpson integration step size must be > 0"
                  << std::endl;
        return false;
    }

    if(delta < 0) {
        std::cout << "Delta must be > 0 MeV"
                  << std::endl;
        return false;
    }

    return true;
}

void CRadCorr::RadiativeCorrection(int iters)
{
    if(!SanityCheck())
        return;

    if(iters > 0)
        iterByNumbers(iters);
    else
        iterByPrecision();
}

void CRadCorr::Radiate()
{
    if(!SanityCheck())
        return;

    for(auto &set : data_sets)
    {
        if(set.non_rad)
            continue;

        std::cout << "Radiate, spectrum energy: " << set.energy << std::endl;

        if(peak_approx)
            radcor(set);
        else
            xyrad2d(set);
    }
}

void CRadCorr::SaveResult(const std::string &path)
{
    std::ofstream output(path);
    if(!output.is_open()) {
        std::cerr << "Cannot open output file "
                  << "\"" << path << "\""
                  << std::endl;
        return;
    }

    output << "#"
           << std::setw(7) << "ENERGY"
           << std::setw(8) << "NU"
           << std::setw(15) << "SIGRAD"
           << std::setw(15) << "SIGBORN"
           << std::setw(15) << "STAT. ERR."
           << std::setw(15) << "SYST. ERR."
//           << std::setw(15) << "ITER. CHANGE"
           << std::endl;
    for(auto &set : data_sets)
    {
        if(set.non_rad)
            continue;

        for(auto &point : set.data)
        {
            double xsr = point.rad*set.weight_mott;
            double xsb = point.born*set.weight_mott;
//            double xsl = point.last*set.weight_mott;
            double scale = xsb/xsr;
            double stat_err = point.stat*scale;
            double syst_err = point.syst*scale;
            double rc_err = set.error*std::abs(xsb - xsr);
            syst_err = sqrt(syst_err*syst_err + rc_err*rc_err);
            output << std::setw(8) << set.energy
                   << std::setw(8) << point.nu
                   << std::setw(15) << xsr
                   << std::setw(15) << xsb
                   << std::setw(15) << stat_err
                   << std::setw(15) << syst_err
//                   << std::setw(15) << xsb - xsl
                   << std::endl;
        }
    }
}

void CRadCorr::iterByNumbers(int iters)
{
    for(int i = 1; i <= iters; ++i)
    {
        for(auto &set : data_sets)
        {
            // do nothing for non-radiated data
            if(set.non_rad)
                continue;

            std::cout << "RC iteration: " << i << " of " << iters
                      << ", spectrum E = " << set.energy
                      << std::endl;
            if(peak_approx)
                radcor(set);
            else
                xyrad2d(set);
        }
    }
}

void CRadCorr::iterByPrecision()
{
    // do first iteration anyway
    bool iter = true;
    int i = 1;

    while(iter)
    {
        for(auto &set : data_sets)
        {
            // do nothing for non-radiated data
            if(set.non_rad)
                continue;

            std::cout << "RC iteration: " << i
                      << ", spectrum E = " << set.energy
                      << std::endl;
            if(peak_approx)
                radcor(set);
            else
                xyrad2d(set);
        }

        // after correction, check if this iteration reaches the precision
        iter = false;
        for(auto &set : data_sets)
        {
            for(auto point : set.data)
            {
                double rel_diff = (point.born - point.last)/point.born;
                // the rel. diff. comes from this iteration is larger than wanted
                // precision
                if(std::abs(rel_diff) > iter_prec)
                {
                    iter = true; // still need more iterations
                    break; // no need to check other points
                }
            }

            if(iter) {
                std::cout << "Iteration is not converging within the precision = "
                          << iter_prec << ", conitnue..."
                          << std::endl;
                break; // no need to check other sets
            }
        }

        ++i;
    }
}


//==============================================================================
// radiative correction for one spectrum                                        
//========================ORIGINAL AUTHORS======================================
// RADCOR FORTRAN CODE                                                          
// ***by Randy Roy Whitney, Phys. Rev. C 9, 2230 - 2235 (1974)                  
//------------------------------------------------------------------------------
// MODIFIED RADCOR - 01/29/2003                                                 
// ***by K. Slifer, Temple University, 01/29/03                                 
// downloadable at  http://www.jlab.org/~slifer/codes.html                      
//------------------------------------------------------------------------------
// MODIFIED RADCOR - 04/01/2007                                                 
// ***by Jaideep Singh                                                          
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
void CRadCorr::radcor(DataSet &set, bool radiate)
{
    spectrum_init(set);

    // iteration on data points
    for(auto &point : set.data)
    {
        point_init(point);

        // components of SIGRAD
        double SIGLOW = 0, SIGBEF = 0, SIGAFT = 0;

        double FBAR = __F_bar(Es, Ep, GAMT);
        SIGLOW = FBAR * std::pow(R*delta/Es, BTB + BTR)
                      * std::pow(delta/Ep, BTA + BTR)
                      * (1. - (XIB+XIA)/delta/(1. - BTB - BTA - 2.*BTR));

        // calculate integral along dEs for fixed Ep SIGBEF
        if((Esmin <= 0) || (Esmax <= 0) || (Esmin >= Esmax)) {
            std::cout << "Skip point at Ep = " << Ep
                      << ", in spectrum E = " << Es
                      << ", ESMIN = " << Esmin
                      << ", ESMAX = " << Esmax
                      << std::endl;
        } else {
            SIGBEF = simpson(Esmin, Esmax, &CRadCorr::fes, this, sim_step, n_sim);
        }

        // calculate integral along dEp for fixed Es SIGAFT
        if((Epmin <= 0) || (Epmax <= 0) || (Epmin >= Epmax)) {
            std::cout << "Skip point at Ep = " << Ep
                      << ", in spectrum E = " << Es
                      << ", EPMIN = " << Esmin
                      << ", EPMAX = " << Esmax
                      << std::endl;
        } else {
            SIGAFT = simpson(Epmin, Epmax, &CRadCorr::fep, this, sim_step, n_sim);
        }

        if(radiate) {
            // radiate, update radiated cross section
            point.rad = SIGLOW*set.weight_mott*point.born + SIGBEF + SIGAFT;
        } else {
            // save last iteration
            point.last = point.born;
            // radiative correction, update the born cross section
            point.born = (point.rad - (SIGBEF+SIGAFT)/set.weight_mott)/SIGLOW;
        }
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
// ***by. Xuefei Yan, Duke University,  10/11/2016                              
//                                                                              
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
void CRadCorr::xyrad2d(DataSet &set, bool radiate)
{
    spectrum_init(set);

    // iteration on data points
    for(auto &point : set.data)
    {
        point_init(point);

        // components of SIGRAD
        double int_2d = 0, sgl_Es = 0, sgl_Ep = 0, sgl_both = 0;
        double FBAR = __F_bar(Es, Ep, GAMT);

        sgl_both = FBAR*std::pow((delta1 + XIB)/Es, BTB + BTR)/gamma(1. + BTB + BTR)
                       *std::pow((delta2 + XIA)/Ep, BTA + BTR)/gamma(1. + BTA + BTR);

        if((Esmin <= 0) || (Esmax <= 0) || (Esmin >= Esmax)) {
            std::cout << "Skip point at Ep = " << Ep
                      << ", in spectrum E = " << Es
                      << ", ESMIN = " << Esmin
                      << ", ESMAX = " << Esmax
                      << std::endl;
        } else {
            int_2d = simpson(Esmin, Esmax, &CRadCorr::int_es, this, sim_step_2d, n_sim_2d);
            sgl_Ep = simpson(Esmin, Esmax, &CRadCorr::int_esdp, this, sim_step, n_sim);
        }

        if((Epmin <= 0) || (Epmax <= 0) || (Epmin >= Epmax)) {
            std::cout << "Skip point at Ep = " << Ep
                      << ", in spectrum E = " << Es
                      << ", EPMIN = " << Esmin
                      << ", EPMAX = " << Esmax
                      << std::endl;
        } else {
            sgl_Es = simpson(Epmin, Epmax, &CRadCorr::int_epds, this, sim_step, n_sim);
        }

        if(radiate) {
            // radiate, update radiated cross section
            point.rad = sgl_both*set.weight_mott*point.born + sgl_Es + sgl_Ep + int_2d;
        } else {
            // save last iteration
            point.last = point.born;
            // radiative correction, update the born cross section
            point.born = (point.rad - (sgl_Es + sgl_Ep + int_2d)/set.weight_mott)/sgl_both;
        }
    }

}

double CRadCorr::int_es(const double &Esx)
{
    double Ep_max = __Ep_max(Esx);
    double Ep_min = delta2 + Ep; // delta2
    if(Ep_max < Ep_min)
        return 0.;

    double lost = __I(Es, Esx, XIB, BTB + BTR);
    double inner_int = simpson(Ep_min, Ep_max, sim_step_2d, n_sim_2d, &CRadCorr::int_ep, this, Esx);
    return lost*inner_int;
}

double CRadCorr::int_ep(const double &Epx, const double &Esx)
{
    double FBAR = __F_bar(Esx, Epx, GAMT);
    double TRx = __btr(Esx, Epx);
    return FBAR*get_cxsn(Esx, Epx)*__I(Epx, Ep, XIA, BTA + TRx);
}

double CRadCorr::int_epds(const double &Epx)
{
    double FBAR = __F_bar(Es, Epx, GAMT);
    double TRx = __btr(Es, Epx);

    double lost = __I(Epx, Ep, XIA, BTA + TRx);
    double int_ds = std::pow((delta1 + XIB)/Epx, BTB + TRx)/gamma(1. + BTB + TRx);
    return lost*int_ds*FBAR*get_cxsn(Es, Epx);
}

double CRadCorr::int_esdp(const double &Esx)
{
    double FBAR = __F_bar(Esx, Ep, GAMT);
    double TRx = __btr(Esx, Ep);

    double lost = __I(Es, Esx, XIB, BTB + TRx);
    double int_dp = std::pow((delta2 + XIA)/Esx, BTA + TRx)/gamma(1. + BTA + TRx);

    return lost*int_dp*FBAR*get_cxsn(Esx, Ep);
}

// interpolates or extrapolates
double CRadCorr::get_cxsn(const double &E0, const double &Eb)
{
    if(Eb >= __Ep_max(E0))
        return 0;

    double WEXC = 1 - Eb/E0;

    // less than the lowest energy we have
    if(E0 <= data_sets.at(0).energy) {
        return interp(data_sets.at(0), WEXC)*F_mott/E0/E0;
    // within the energy range in spectrum
    } else if (E0 <= data_sets.back().energy) {
        size_t i = 1;
        for(; i < data_sets.size(); ++i)
        {
            if(data_sets.at(i).energy >= E0)
                break;
        }
        double E1 = data_sets.at(i-1).energy;
        double E2 = data_sets.at(i).energy;
        double FTCS = (interp(data_sets.at(i-1), WEXC)*(E2-E0)
                    + interp(data_sets.at(i), WEXC)*(E0-E1))/(E2-E1);
        return F_mott/E0/E0*FTCS;
    } else {
        std::cerr << "Required energy E0 = " << E0
                  << " exceeds the highest energy in data Emax = "
                  << data_sets.back().energy
                  << std::endl;
        exit(-1);
    }
    return 0.;
}

// recursive binary_search for terp
int __rc_binary_search(const CRadCorr::DataPoint *array, const double &key, const int &first, const int &last)
{
    if(first >= last)
        return -1;

    int mid = (first + last - 1)/2;

    if(key < array[mid].v)
        return __rc_binary_search(array, key, first, mid);
    else if(key > array[mid+1].v)
        return __rc_binary_search(array, key, mid+1, last);
    else
        return mid;
}

double CRadCorr::interp(const DataSet &set, const double &w)
{
    // exceeds the boundary
    if(w < set.data.front().v)
        return 0.;
    if(w >= set.data.back().v)
        return set.data.back().born;

    int j = __rc_binary_search(&set.data[0], w, 0, set.data.size());
    if(j < 0)
        return 0;

    // do 2 points interpolation
    if(j <= 1) {
        const DataPoint &p1 = set.data.at(j);
        const DataPoint &p2 = set.data.at(j+1);

        double _interp = p1.born*(p2.v - w) + p2.born*(w - p1.v);
        // unknown constant 0.00001 here, probably some protection for two same points
        return _interp/(p2.v - p1.v);// + 0.00001);
    }

    // 3 point parabolic fit
    const DataPoint &p1 = set.data.at(j-1);
    const DataPoint &p2 = set.data.at(j);
    const DataPoint &p3 = set.data.at(j+1);
    double x, xp, xm, a, b, c;
    x = w - p2.v;
    xp = p3.v - p2.v;
    xm = p1.v - p2.v;
    a = p2.born;
    c = (p3.born - p1.born)*(xp + xm)/(xp - xm) - (p3.born + p1.born - 2*p2.born);
    c /= xp*xm*2;
    b = (p3.born - p1.born)/(xp - xm) - c*(xp + xm);
    return a + b*x + c*x*x;
}

// calculate XI based on STEIN's formula
void CRadCorr::calculateXI(DataSet &set)
{
    // formula from STEIN, but old and probably wrong estimation
    double xi = (PI*ELECM/2/ALPHA)*(set.radl_before + set.radl_after);
    xi /= (target_Z + __eta(target_Z));
    xi /= log(183*std::pow(target_Z, -1./3.));

    set.coll_before = xi/2;
    set.coll_after = xi/2;
}

// spectrum based kinematics intialization
inline void CRadCorr::spectrum_init(DataSet &set)
{
    // update these parameters when external RC is ON
    if(external_RC) {
        // originally b(z) = 4./3., which is an approximated value regardless of
        // Z dependence, it has about 2% difference with the value from Eq. A45
        // in STEIN, thus b(z) is updated according to that equation.
        BTB = set.radl_before*Bz;
        BTA = set.radl_after*Bz;
        XIB = set.coll_before;
        XIA = set.coll_after;
    } else {
        BTB = 0;
        BTA = 0;
        XIB = 0;
        XIA = 0;
    }

    Es = set.energy;
    // improvements by J. Singh
    // two gamma function normalization so
    // 1/gamma(1 + b*tb)/gamma(1 + b*ta) = 1 + 0.5772b*(tb + ta) + ...
    GAMT = gamma(1. + BTB) * gamma(1. + BTA);
}

// data point based kinematics initialization
inline void CRadCorr::point_init(DataPoint &point)
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
    return log(__Q2(_E, _Epr)/ELECM/ELECM);
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
    if(peak_approx)
        DHO += -0.5*Log2EsEp;               // Correction to peaking approx. ?

    DHO *= ALPHA/PI;                        // common factor

    return exp(DHO)/_gamma_t;
}

// Get tr(Q2), effective radiator thickness before and after the scattering from
// external Bremsstrahlung process, used inline __log_Q2m2
inline double CRadCorr::__btr(double _E, double _Epr)
{
    if(!internal_RC)
        return 0.;

    return ALPHA/PI*(__log_Q2m2(_E, _Epr) - 1.);
}

// Probability function I(E0, E, t)
// E0 is replaced by E0 - delta to include the ffect from ionization
// NOTICE here we are using b(z)t instead of t
inline double CRadCorr::__I(double _E0, double _E, double _D, double _bt)
{
    double _v = (_E0 - _D - _E)/_E0;
    return _bt/gamma(1. + _bt)*std::pow(_v, _bt)/(_E0 - _D - _E)*__phi(_v);
}

// read experimental data in the format line by line
// comment marks are # and //
// splitter can be tab, space and comma
// additional white spaces (space and tab) will be trimmed
// *data header*
// *For non-radiated data*
// energy, rel. error, norm
// *For radiated data*
// energy, rel. error, norm, radl bef, radl aft, coll thick bef, coll thick aft
// *data*
// energy
// nu, cxsn, stat. error, syst. error
// nu, cxsn, stat. error, syst. error
// ......
void CRadCorr::readData(ConfigParser &c_parser)
{
    int set_idx = -1;
    double energy, radl_bef, radl_aft, coll_bef, coll_aft, error, norm;
    double nu, cxsn, stat, syst;

    while(c_parser.ParseLine())
    {
        if (c_parser.NbofElements() == 4) { // new data points

            if(set_idx == -1) {
                std::cout << "Cannot determine which data set the data belong to "
                          << "at line " << c_parser.LineNumber()
                          << std::endl
                          << "\"" << c_parser.CurrentLine() << "\""
                          << std::endl;
                continue;
            }

            DataSet &set = data_sets.at(set_idx);

            c_parser >> nu >> cxsn >> stat >> syst;
            // apply normalization
            cxsn *= set.normalization/set.weight_mott;
            stat *= set.normalization;
            DataPoint new_point(nu, cxsn, stat, syst);
            // calculate ep
            new_point.Ep = set.energy - nu;
            new_point.v = nu/set.energy;

            set.data.push_back(std::move(new_point));

        } else if(c_parser.NbofElements() == 1) { // data set indication

            c_parser >> energy;
            size_t i = 0;
            // find the proper data set for inserting data points
            for(; i < data_sets.size(); ++i)
            {
                if(data_sets.at(i).energy == energy) {
                    set_idx = i;
                    break;
                }
            }
            // did not find the energy in data_sets
            if(i >= data_sets.size()) {
                std::cout << "Cannot find data set definition for energy "
                          << energy
                          << ", please make sure you add this energy in the data header"
                          << std::endl;
                set_idx = -1;
            }

        } else if(c_parser.NbofElements() == 3) { // new data set (non-radiated)

            c_parser >> energy >> error >> norm;

            DataSet new_set(energy, error, norm);
            new_set.weight_mott = F_mott/energy/energy;

            data_sets.push_back(std::move(new_set));

        } else if(c_parser.NbofElements() == 7) { // new data set

            c_parser >> energy >> error >> norm
                     >> radl_bef >> radl_aft >> coll_bef >> coll_aft;

            DataSet new_set(energy, radl_bef, radl_aft, coll_bef, coll_aft, error, norm);
            new_set.weight_mott = F_mott/energy/energy;

            data_sets.push_back(std::move(new_set));

        } else {
            std::cout << "Skipped line " << c_parser.LineNumber()
                      << ", format is unrecognized." << std::endl
                      << "\"" << c_parser.CurrentLine() << "\""
                      << std::endl;
        }
    }

}

#define __GAMMA_NUM 9
#define __GAMMA_G 7.0
static double __gamma_c[] = {0.9999999999998099,
                             6.7652036812188510E2,
                            -1.2591392167224028E3,
                             7.7132342877765313E2,
                            -1.7661502916214059E2,
                             1.2507343278686905E1,
                            -1.3857109526572012E-1,
                             9.9843695780195716E-6,
                             1.5056327351493116E-7};

// gamma function
double CRadCorr::gamma(const double &z)
{

    if(z < 1) {
        return PI*(1 - z)/sin(PI*(1 - z))/gamma(2 - z);
    } else if(z == 1) {
        return 1;
    } else {
        double ag = __gamma_c[0];
        for(int k = 1; k < __GAMMA_NUM; ++k)
            ag += __gamma_c[k]/(z - 1. + k);

        double output = 0.5*log(2*PI)
                        + (z - 0.5)*log(z - 0.5 + __GAMMA_G)
                        - (z - 0.5 + __GAMMA_G)
                        + log(ag);
        return exp(output);
    }
}

#define __SPENCE_NUM 8
#define __SPENCE_NMAX 50
#define __SPENCE_NMAX_PREC 1000
static double __spence_c[] = {-1.1741940560772957946600E-1,
                              -2.7618966846029390643791E-2,
                              -8.0493987190845793511240E-3,
                              -2.7095568666150792944136E-3,
                              -7.1455906877666711465857E-4,
                               4.1757495974272487715417E-2,
                              -4.9028996486663818655E-2,
                              -9.08073640732783360};

// spence function
double CRadCorr::spence(const double &z, const double &res)
{
    if(z > 1) {
        return 2*PI*PI/6 - log(z)*log(z)/2 - spence(1/z, res);
    } else if (z == 1) {
        return PI*PI/6;
    } else if (z > 0.5) {
        return PI*PI/6 - log(1-z)*log(z) - spence(1-z, res);
    } else if (z == 0.5) {
        return PI*PI/6/2 - log(0.5)*log(0.5)/2;
    } else if (z > 0) {
        return spence_tr(z, res, __SPENCE_NMAX); // do nothing, fall into the bottom session
    } else if (z == 0) {
        return 0;
    } else if (z > -0.95) {
        return spence_tr(z, res, __SPENCE_NMAX_PREC);
    } else if (z > -1.05) {
        // poly fit
        double output = 0;
        double dz = z + 1;
        for(int i = 0; i < __SPENCE_NUM; ++i)
            output += __spence_c[i]*std::pow(dz, i);
        return -(1 + output*dz*dz)*PI*PI/6/2 + dz*log(2);
    } else {
        return -PI*PI/6 - log(-z)*log(-z)/2 - spence(1/z, res);
    }
}

// truncation of spence function
double CRadCorr::spence_tr(const double &z, const double &res, const int &nmax)
{
    // calculate spence until res is reached
    double output = 0.;
    int n = 0;
    while(++n <= nmax)
    {
        double nth = std::pow(z, n)/n/n;
        output += nth;

        if(std::abs(nth) < res)
            break;
    }

    return output;
}

