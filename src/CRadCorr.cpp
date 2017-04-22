#include "CRadCorr.h"
#include "canalib.h"
#include "ConfigParser.h"
#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>



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
    angle *= cana::deg2rad;
    target_M = target_A * cana::amu;
    sin2 = sin(angle/2.)*sin(angle/2.);
    cos2 = cos(angle/2.)*cos(angle/2.);

    // mott cross section
    F_mott = cana::hbarc*cana::alpha*cos(angle/2)/2/sin2;
    F_mott = F_mott*F_mott*1e7; // MeV^2*nb/sr

    // Schwinger term in internal radiation
    Schwinger = cana::pi*cana::pi/6 - cana::spence(cos2);

    // B(z) Eq. A45 in STEIN
    Bz = 1./9.*(target_Z + 1.)/(target_Z + __eta(target_Z));
    Bz /= log(183.*std::pow(target_Z, -1./3.));
    Bz = 4./3.*(1. + Bz);

    // read data files
    std::string file_str = GetConfig<std::string>("Data Config File");
    auto files = ConfigParser::split(file_str, ",");

    for(auto &f : files)
    {
        readDataConf(f);
    }

    // sort data sets by energy
    std::sort(data_sets.begin(), data_sets.end(),
             [] (const DataSet &set1, const DataSet &set2)
             {
                 return set1.energy < set2.energy;
             });

    // show how many data sets are read-in
    std::cout << "Read " << data_sets.size() << " data sets."
              << std::endl;

    for(auto &dset : data_sets)
    {
        std::cout << "Energy: " << std::setw(8) << dset.energy
                  << ",    DataPoints:" << std::setw(8) << dset.data.size()
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
    for(auto &dset : data_sets)
    {
        if(!dset.non_rad) {
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

    bool end_by_prec = (iters <= 0)? true : false;

    for(int iter = 1; (iter <= iters) || end_by_prec; ++iter)
    {
        //init_model();

        // do radiative correction for all data sets
        for(auto &dset : data_sets)
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
            for(auto &dset : data_sets)
            {
                for(auto &point : dset.data)
                {
                    double rel_diff = fabs(1. - point.last/point.born);
                    if(rel_diff > max_rel_diff)
                        max_rel_diff = rel_diff;
                }
            }

            // check if it should continue or not
            if(max_rel_diff > iter_prec) {
                std::cout << "Iteration " << iter << " is not converging within "
                          << "the required precision = " << iter_prec*100. <<"%, "
                          << "the maximum difference = " << max_rel_diff*100. << "%, "
                          << "continue iterations..."
                          << std::endl;
            } else {
                // Done the radiative correction
                return;
            }
        }

    }
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
    for(auto &dset : data_sets)
    {
        if(dset.non_rad)
            continue;

        for(auto &point : dset.data)
        {
            double xsr = point.rad*dset.weight_mott;
            double xsb = point.born*dset.weight_mott;
//            double xsl = point.last*set.weight_mott;
            double scale = xsb/xsr;
            double stat_err = point.stat*scale;
            double syst_err = point.syst*scale;
            double rc_err = dset.error*std::abs(xsb - xsr);
            syst_err = sqrt(syst_err*syst_err + rc_err*rc_err);
            output << std::setw(8) << dset.energy
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
void CRadCorr::radcor(DataSet &s, bool radiate)
{
    spectrum_init(s);

    int count = 0;
    // iteration on data points
    for(auto &point : s.data)
    {
        std::cout << "Data points: " << count++ << "/" << s.data.size() << "\r"
                  << std::flush;
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
    std::cout << "Data points: " << count << "/" << s.data.size() << std::endl;
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
void CRadCorr::xyrad2d(DataSet &s, bool radiate)
{
    spectrum_init(s);

    int count = 0;
    // iteration on data points
    for(auto &point : s.data)
    {
        std::cout << "Data points: " << count++ << "/" << s.data.size() << "\r"
                  << std::flush;
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
    std::cout << "Data points: " << count << "/" << s.data.size() << std::endl;

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

    double weight = (E0 - Eb)/E0;

    // less than the lowest energy we have
    if(E0 < data_sets.at(0).energy) {
        // extrapolate
        return interp(data_sets.at(0), weight)*F_mott/E0/E0;
        //return from_model(E0, Eb);
    // within the energy range in spectrum
    } else if (E0 <= data_sets.back().energy) {
        size_t i = 0;
        for(; i < data_sets.size(); ++i)
        {
            if(data_sets.at(i).energy >= E0)
                break;
        }
        if(data_sets.at(i).energy == E0)
            return interp(data_sets.at(i), weight)*F_mott/E0/E0;

        double E1 = data_sets.at(i-1).energy;
        double E2 = data_sets.at(i).energy;
        double FTCS = (interp(data_sets.at(i-1), weight)*(E2-E0)
                       + interp(data_sets.at(i), weight)*(E0-E1))/(E2-E1);
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

double CRadCorr::interp(const DataSet &s, const double &w)
{
    // exceeds the boundary
    if(w < s.data.front().v)
        return 0.;
    if(w >= s.data.back().v)
        return s.data.back().born;

    // search the position of w
    auto it_pair = cana::binary_search_interval(s.data.begin(), s.data.end(), w);

    // not found
    if(it_pair.second == s.data.end() || it_pair.first == s.data.end())
        return 0;

    // exact matched
    if(it_pair.first == it_pair.second)
        return it_pair.first->born;

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
    return a + b*x + c*x*x;
}

// initialize the model
// determine the scale and shift from the data at lowest energy
void CRadCorr::init_model()
{
    // data sets should be in a energy transcendent order
    for(unsigned int i = 0; i < data_sets.size(); ++i)
    {
        if(!data_sets.at(i).non_rad) {
            find_model_scale(data_sets.at(i));
            break;
        }
    }
}

// Determines the scale and shift from the peak value, so the data set should has
// QE peak or delta resonance peak, otherwise it will be problematic
void CRadCorr::find_model_scale(const DataSet &mset)
{
    // find the peak point
    unsigned int ip = 0;
    double max_xs = 0.;
    for(unsigned int i = 0; i < mset.data.size(); ++i)
    {
        const DataPoint &p = mset.data.at(i);
        if(p.born > max_xs) {
            ip = i;
            max_xs = p.born;
        }
    }

    // the maximum point
    const DataPoint &mpoint = mset.data.at(ip);

    max_xs = 0.;
    double cxsn = 0; //, max_nu = 0.;
    // we search the QE peak in the model for a +-20 MeV range
    for(double nu = std::max(5., mpoint.nu - 20.); nu < mpoint.nu + 20.; nu += 1.0)
    {
        QFS_xs(target_Z, target_A, mset.energy, (mset.energy - nu), angle, &cxsn);
        if(max_xs < cxsn)
        {
            max_xs = cxsn;
            //max_nu = nu;
        }
    }

    model_scale = mpoint.born*mset.weight_mott/max_xs;
    // do not apply shift now, need more study
//    model_shift = mpoint.nu - max_nu;

    std::cout << "Determined model scale and shift from spectrum at "
              << mset.energy << " MeV. "
              << std::endl
              << "Scale = " << model_scale << ", "
              << "shift = " << model_shift
              << std::endl;
}

// get cross sections from model
inline double CRadCorr::from_model(const double &E0, const double &Eb)
{
    // bosted model
    double cxsn;
    QFS_xs(target_Z, target_A, E0, Eb - model_shift, angle, &cxsn);
    //Bosted_xs(target_Z, target_A, E0, Eb - model_shift, angle, &cxsn);
    return cxsn*model_scale;
}

// spectrum based kinematics intialization
inline void CRadCorr::spectrum_init(DataSet &s)
{
    // update these parameters when external RC is ON
    if(external_RC) {
        // LATEST CHANGE: replaced b(z) = 4./3. with Eq. A45 from STEIN
        // originally it was an approximated value regardless of Z dependence,
        // ignoring it introducing an error on a few percent level.
        BTB = (s.radl_before + __ice_radl(s.ice_before))*Bz;
        BTA = (s.radl_after + __ice_radl(s.ice_after))*Bz;
        if(user_defined_XI) {
            XIB = s.coll_before + __ice_coll(s.ice_before);
            XIA = s.coll_after + __ice_coll(s.ice_after);
        } else {
            XIB = __XI_Stein(s.radl_before + __ice_radl(s.ice_before));
            XIA = __XI_Stein(s.radl_after + __ice_radl(s.ice_after));
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

// calculate ice radiation length, input is thickness in mm
inline double CRadCorr::__ice_radl(double thickness)
{
    // hard coded for ice ONLY
    // http://pdg.lbl.gov/2008/AtomicNuclearProperties/HTML_PAGES/325.html
    return thickness/393.1;
}

// calculate ice collisional loss, input is thickness in mm
// it assumes beta is close to 1, so for electron beam, the beam energy should
// be far greater than electron's mass
inline double CRadCorr::__ice_coll(double thickness)
{
    // hard coded for ice ONLY
    // See J. Singh's technical note for details
    // http://hallaweb.jlab.org/experiment/E97-110/tech/radlength_sagdhv130.pdf
    // Z/A = 0.55509 mol/g
    // a = 0.15353747 MeV*(cm^2/mol)
    // rho = 0.918 g/cm^3
    return 0.55509*0.15353747*0.918*thickness/10.;

}

void CRadCorr::readDataConf(const std::string &path)
{
    ConfigParser c_parser;
    auto contents = c_parser.ReadFileInBlocks(path, "{", "}", "DATASET", false);

    for(auto &config : contents)
    {
        createDataSet(config, path);
    }
}

template<typename T>
inline void update_config(const ConfigObject &conf, const std::string &var, T &val)
{
    ConfigValue term = conf.GetConfigValue(var);
    if(!term.IsEmpty())
        val = term.Convert<T>();
}

// create a new data set based on configuration text
void CRadCorr::createDataSet(const std::string &config, const std::string &path)
{
    ConfigObject conf_obj;

    conf_obj.ReadConfigString(config);

    ConfigValue Econf = conf_obj.GetConfigValue("Energy");
    if(Econf.IsEmpty() || cana::is_in(Econf.Double(), data_sets.begin(), data_sets.end())) {
        std::cerr << "Error: Missing energy information or there exists a data"
                  << "set with the same energy in configuration: " << config
                  << std::endl;
        return;
    }

    // create a data set with default values
    DataSet new_set(Econf.Double());

    new_set.weight_mott = F_mott/std::pow(new_set.energy, 2);

    // connect data set variables to configurations
    update_config(conf_obj, "Radiation Length Before", new_set.radl_before);
    update_config(conf_obj, "Radiation Length After", new_set.radl_after);
    update_config(conf_obj, "Collisional Loss Before", new_set.coll_before);
    update_config(conf_obj, "Collisional Loss After", new_set.coll_after);
    update_config(conf_obj, "Ice Before", new_set.ice_before);
    update_config(conf_obj, "Ice After", new_set.ice_after);
    update_config(conf_obj, "RC Error", new_set.error);
    update_config(conf_obj, "Born Level", new_set.non_rad);
    update_config(conf_obj, "Normalization", new_set.normalization);

    // read-in the data points in this set
    std::string data_file, data_label;
    update_config(conf_obj, "Data File", data_file);
    update_config(conf_obj, "Data Label", data_label);

    std::string dir = ConfigParser::decompose_path(path).dir;
    readData(new_set, ConfigParser::form_path(dir, data_file), data_label);

    // sort data points by nu
    std::sort(new_set.data.begin(), new_set.data.end(),
              [] (const DataPoint &p1, const DataPoint &p2)
              {
                  return p1.nu < p2.nu;
              });

    data_sets.emplace_back(std::move(new_set));
}

// read experimental data in the format line by line
// comment marks are # and //
// splitter can be tab, space and comma
// additional white spaces (space and tab) will be trimmed
// expected 4 ~ 5 columns (label can be neglected if input data_label is empty)
// *label*, nu, cxsn, stat. error, syst. error
void CRadCorr::readData(DataSet &dset, const std::string &path, const std::string &label)
{
    ConfigParser c_parser;

    if(!c_parser.ReadFile(path)) {
        std::cerr << "Cannot open file \"" << path << "\", no data points read "
                  << "for data set at energy = " << dset.energy
                  << std::endl;
        return;
    }

    while(c_parser.ParseLine())
    {
        // check if label agrees
        if(c_parser.NbofElements() == 5) {
            std::string plabel = c_parser.TakeFirst().String();
            if(plabel != label)
                continue;
        }

        // check format
        if(!c_parser.CheckElements(4))
            continue;

        // new data points
        double nu, cxsn, stat, syst;
        c_parser >> nu >> cxsn >> stat >> syst;
        // apply normalization
        cxsn *= dset.normalization/dset.weight_mott;
        stat *= dset.normalization;

        DataPoint new_point(nu, cxsn, stat, syst);
        // calculate ep
        new_point.Ep = dset.energy - nu;
        new_point.v = nu/dset.energy;

        dset.data.push_back(std::move(new_point));
    }
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
