#include "CRadCorr.h"
#include "canalib.h"
#include "ConfigParser.h"
#include <iostream>
#include <iomanip>
#include <algorithm>

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
    // the configuration value is supposed to be provided in config file
    // if the term is not included in configuration file, it will use the
    // default value
    readConfigFile(path);

    // get configuration value from the file
    // functions
    internal_RC = getDefConfig<bool>("Internal RC", true);
    external_RC = getDefConfig<bool>("External RC", true);
    user_defined_XI = getDefConfig<bool>("User Defined XI", true);
    peak_approx = getDefConfig<bool>("Peaking Approximation", false);

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
    // target
    target_Z = getDefConfig<double>("Target Z", 2);
    target_A = getDefConfig<double>("Target A", 3.0149322473);
    // scattering angle
    angle = getDefConfig<double>("Scattering Angle", 6.1);

    // calculate values based on the input configuration
    // common value for all spectrums
    angle /= RADDEG;
    target_M = target_A * AMUMEV;
    sin2 = sin(angle/2.)*sin(angle/2.);
    cos2 = cos(angle/2.)*cos(angle/2.);

    // mott cross section
    F_mott = HBARC*ALPHA*cos(angle/2)/2/sin2;
    F_mott = F_mott*F_mott*1e7; // MeV^2*nb/sr

    // Schwinger term in internal radiation
    Schwinger = PI*PI/6 - cana::spence(cos2);

    // B(z) Eq. A45 in STEIN
    Bz = 1./9.*(target_Z + 1.)/(target_Z + __eta(target_Z));
    Bz /= log(183.*std::pow(target_Z, -1./3.));
    Bz = 4./3.*(1. + Bz);

    // read data files
    std::string file_str = GetConfig<std::string>("Data File");
    auto files = ConfigParser::split(file_str, ",");

    std::vector<std::string> flist;
    while(files.size())
    {
        flist.push_back(ConfigParser::trim(files.front(), " \t"));
        files.pop();
    }
    if(!flist.empty())
        ReadExpData(flist);
}

// read one data file, see comments in readData() for format details
void CRadCorr::ReadExpData(const char *path)
{
    std::vector<std::string> flist;
    flist.emplace_back(path);
    ReadExpData(flist);
}

void CRadCorr::ReadExpData(const std::string &path)
{
    std::vector<std::string> flist;
    flist.push_back(path);
    ReadExpData(flist);
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

    // sort data points by nu
    for(auto &s : data_sets)
    {
        std::sort(s.data.begin(), s.data.end(),
                  [] (const DataPoint &p1, const DataPoint &p2)
                  {
                    return p1.nu < p2.nu;
                  });
    }

    // calculate collision thickness if not using users' input
    if(!user_defined_XI) {
        for(auto &s : data_sets)
            calculateXI(s);
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

// test some configuration values, warn some simple mistake
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
void CRadCorr::RadiativeCorrection(int iters)
{
    if(!SanityCheck())
        return;

    if(iters > 0)
        iterByNumbers(iters);
    else
        iterByPrecision();
}

// radiate data sets
void CRadCorr::Radiate()
{
    if(!SanityCheck())
        return;

    for(auto &s : data_sets)
    {
        if(s.non_rad)
            continue;

        std::cout << "Radiate, spectrum energy: " << s.energy << std::endl;

        if(peak_approx)
            radcor(s, true);
        else
            xyrad2d(s, true);
    }
}

// save result to the path
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
        for(auto &s : data_sets)
        {
            // do nothing for non-radiated data
            if(s.non_rad)
                continue;

            std::cout << "RC iteration: " << i << " of " << iters
                      << ", spectrum E = " << s.energy
                      << std::endl;
            if(peak_approx)
                radcor(s);
            else
                xyrad2d(s);
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
        for(auto &s : data_sets)
        {
            // do nothing for non-radiated data
            if(s.non_rad)
                continue;

            std::cout << "RC iteration: " << i
                      << ", spectrum E = " << s.energy
                      << std::endl;
            if(peak_approx)
                radcor(s);
            else
                xyrad2d(s);
        }

        // after correction, check if this iteration reaches the precision
        iter = false;
        for(auto &s : data_sets)
        {
            for(auto point : s.data)
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
//========================LATEST CHANGES========================================
// 1. Adapted the code in C++                                                   
// 2. Replaced B(z) with analytical form                                        
// 3. Corrected the calculation in cross section interpolation                  
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
void CRadCorr::radcor(DataSet &s, bool radiate)
{
    spectrum_init(s);

    // iteration on data points
    for(auto &point : s.data)
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
            SIGBEF = cana::simpson(Esmin, Esmax, &CRadCorr::fes, this, sim_step, n_sim);
        }

        // calculate integral along dEp for fixed Es SIGAFT
        if((Epmin <= 0) || (Epmax <= 0) || (Epmin >= Epmax)) {
            std::cout << "Skip point at Ep = " << Ep
                      << ", in spectrum E = " << Es
                      << ", EPMIN = " << Esmin
                      << ", EPMAX = " << Esmax
                      << std::endl;
        } else {
            SIGAFT = cana::simpson(Epmin, Epmax, &CRadCorr::fep, this, sim_step, n_sim);
        }

        if(radiate) {
            // radiate, update radiated cross section
            point.rad = SIGLOW*point.born + (SIGBEF + SIGAFT)/s.weight_mott;
        } else {
            // save last iteration
            point.last = point.born;
            // radiative correction, update the born cross section
            point.born = (point.rad - (SIGBEF+SIGAFT)/s.weight_mott)/SIGLOW;
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
// ***by Xuefei Yan, Duke University,  10/11/2016                               
//========================LATEST CHANGES========================================
// 1. Adapted the code in C++                                                   
// 2. Now get cross section from interpolation of input data instead of model   
// 3. Add external radiative correction part                                    
// 4. Replaced B(z) with analytical form                                        
// 5. Using exponentialted higher order term and gamma term in Fbar calculation 
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
void CRadCorr::xyrad2d(DataSet &s, bool radiate)
{
    spectrum_init(s);

    // iteration on data points
    for(auto &point : s.data)
    {
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
            std::cout << "Skip point at Ep = " << Ep
                      << ", in spectrum E = " << Es
                      << ", ESMIN = " << Esmin
                      << ", ESMAX = " << Esmax
                      << std::endl;
        } else {
            int_2d = cana::simpson(Esmin, Esmax, &CRadCorr::int_es, this, sim_step_2d, n_sim_2d);
            sgl_Ep = cana::simpson(Esmin, Esmax, &CRadCorr::int_esdp, this, sim_step, n_sim);
        }

        if((Epmin <= 0) || (Epmax <= 0) || (Epmin >= Epmax)) {
            std::cout << "Skip point at Ep = " << Ep
                      << ", in spectrum E = " << Es
                      << ", EPMIN = " << Esmin
                      << ", EPMAX = " << Esmax
                      << std::endl;
        } else {
            sgl_Es = cana::simpson(Epmin, Epmax, &CRadCorr::int_ep, this, sim_step, n_sim);
            // Es singularity
            sgl_Es *= std::pow(delta1/Es, BTB + BTR)/cana::gamma(1. + BTB + BTR);
            // coll. loss term
            sgl_Es *= 1. - XIB/(1 - BTB - BTR)/delta1;
        }

        if(radiate) {
            // radiate, update radiated cross section
            point.rad = sgl_both*point.born + (sgl_Es + sgl_Ep + int_2d)/s.weight_mott;
        } else {
            // save last iteration
            point.last = point.born;
            // radiative correction, update the born cross section
            point.born = (point.rad - (sgl_Es + sgl_Ep + int_2d)/s.weight_mott)/sgl_both;
        }
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
double CRadCorr::get_cxsn(const double &E0, const double &Eb)
{
    if(Eb >= __Ep_max(E0))
        return 0;

    double weight = 1. -  Eb/E0;

    // search the position of E0 in data sets
    auto it_pair = cana::binary_search_interval(data_sets.begin(), data_sets.end(), E0);
    // out of spectrum
    if(it_pair.second == data_sets.end() || it_pair.first == data_sets.end()) {
        return from_model(E0, Eb);
    // exact matched
    } else if(it_pair.first == it_pair.second) {
        return interp(*(it_pair.first), weight)*F_mott/E0/E0;
    // find in between
    } else {
        double E1 = it_pair.first->energy;
        double E2 = it_pair.second->energy;
        double FTCS = (interp(*(it_pair.first), weight)*(E2-E0)
                       + interp(*(it_pair.second), weight)*(E0-E1))/(E2-E1);
        return F_mott/E0/E0*FTCS;
    }
}

double CRadCorr::interp(const DataSet &s, const double &w)
{
    // search the position of w
    auto it_pair = cana::binary_search_interval(s.data.begin(), s.data.end(), w);

    // not found
    if(it_pair.second == s.data.end() || it_pair.first == s.data.end()) {
        // return cross section normalized by Mott, to be consistent
        //double cxsn = from_model(s.energy, (1. - w)*s.energy);
        //return cxsn/(F_mott/s.energy/s.energy);
        return 0.;
    }

    // exact matched
    if(it_pair.first == it_pair.second)
        return (it_pair.first->born);

    // only have 2 points, do a straight line interpolation
    if(it_pair.first == s.data.begin()) {
        const DataPoint &p1 = *it_pair.first;
        const DataPoint &p2 = *it_pair.second;

        double _interp = p1.born*(p2.v - w) + p2.born*(w - p1.v);

        // LATEST CHANGE: removed unknown constant 0.00001 here, 
        // probably some protection for two same points
        // return _interp/(p2.v - p1.v + 0.00001);
        return _interp/(p2.v - p1.v);
    }

    // 3 points parabolic fit
    const DataPoint &p1 = *(it_pair.first - 1);
    const DataPoint &p2 = *it_pair.first;
    const DataPoint &p3 = *it_pair.second;
    double x, xp, xm, a, b, c;
    x = w - p2.v;
    xp = p3.v - p2.v;
    xm = p1.v - p2.v;
    a = p2.born;
    c = (p3.born - p1.born)*(xp + xm)/(xp - xm) - (p3.born + p1.born - 2*p2.born);
    c /= xp*xm*2;
    b = (p3.born - p1.born)/(xp - xm) - c*(xp + xm);
    return (a + b*x + c*x*x);
}

inline double CRadCorr::from_model(const double &E0, const double &Eb)
{
    // bosted model
    double cxsn;
    Bosted_xs(target_Z, target_A, E0/1000., Eb/1000., angle, &cxsn);
    return 1000.*cxsn;
}

// calculate XI based on STEIN's formula
void CRadCorr::calculateXI(DataSet &s)
{
    // formula from STEIN, but old and probably wrong estimation
    double xi = (PI*ELECM/2/ALPHA)*(s.radl_before + s.radl_after);
    xi /= (target_Z + __eta(target_Z));
    xi /= log(183*std::pow(target_Z, -1./3.));

    s.coll_before = xi/2;
    s.coll_after = xi/2;
}

// spectrum based kinematics intialization
inline void CRadCorr::spectrum_init(DataSet &s)
{
    // update these parameters when external RC is ON
    if(external_RC) {
        // LATEST CHANGE: replaced b(z) = 4./3. with Eq. A45 from STEIN
        // originally it was an approximated value regardless of Z dependence,
        // ignoring it introducing an error on a few percent level.
        BTB = s.radl_before*Bz;
        BTA = s.radl_after*Bz;
        XIB = s.coll_before;
        XIA = s.coll_after;
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
// __XI is accounted for collisional loss 
// NOTICE here we are using b(z)t instead of t
inline double CRadCorr::__I(double _E0, double _E, double _XI, double _bt)
{
    double _dE = _E0 - _E;
    double _v = _dE/_E0;
    return _bt/cana::gamma(1. + _bt)*std::pow(_v, _bt)/_dE*(__phi(_v) + _XI/2./_dE);
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

            DataSet &s = data_sets.at(set_idx);

            c_parser >> nu >> cxsn >> stat >> syst;
            // apply normalization
            cxsn *= s.normalization/s.weight_mott;
            stat *= s.normalization;
            DataPoint new_point(nu, cxsn, stat, syst);
            // calculate ep
            new_point.Ep = s.energy - nu;
            new_point.v = nu/s.energy;

            s.data.push_back(std::move(new_point));

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

            number_operation("Radiation Length Before", radl_bef);
            number_operation("Radiation Length After", radl_aft);

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

template<typename T>
void CRadCorr::number_operation(const std::string &key, T &val)
{
    auto operations = ConfigParser::split(GetConfig<std::string>(key), ",");

    while(operations.size())
    {
        std::string op = ConfigParser::trim(operations.front(), " \t");
        operations.pop();
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
