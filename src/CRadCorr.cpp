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
    n_simpson_bins = getConfig<int>("Min Number of Simpson Bins", 10000);
    simpson_bin_size = getConfig<double>("Simpson Step Size", 0.1);
    delta = getConfig<double>("IR DIV Delta", 10);
    target_Z = getConfig<double>("Target Z", 2);
    target_A = getConfig<double>("Target A", 3.0149322473);
    angle = getConfig<double>("Scattering Angle", 6.1);

    // common value for all spectrums
    angle /= RADDEG;
    target_M = target_A * AMUMEV;
    sin2 = sin(angle/2)*sin(angle/2);
    F_mott = HBARC*ALPHA*cos(angle/2)/2/sin2;
    F_mott = F_mott*F_mott*1e7; // MeV^2*nb/sr
    cos2 = cos(angle/2)*cos(angle/2);
    Schwinger = PI*PI/6 - spence(cos2);
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

    if(simpson_bin_size < 0) {
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

    for(int i = 1; i <= iters; ++i)
    {
        for(auto &set : data_sets)
        {
            // do nothing for non-radiated data
            if(set.non_rad)
                continue;

            std::cout << "RC iteration: " << i
                      << ", spectrum energy: " << set.energy
                      << std::endl;

            radcor(set);
        }
    }
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

        radcor(set, true);
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

    for(auto &set : data_sets)
    {
        if(set.non_rad)
            continue;

        for(auto &point : set.data)
        {
            double xsi = point.rad*set.weight_mott;
            double xsf = point.born*set.weight_mott;
            double scale = xsf/xsi;
            double stat_err = point.stat*scale;
            double syst_err = point.syst*scale;
            double rc_err = set.error*std::abs(xsf - xsi);
            syst_err = sqrt(syst_err*syst_err + rc_err*rc_err);
            output << std::setw(8) << set.energy
                   << std::setw(8) << point.nu
                   << std::setw(15) << xsi
                   << std::setw(15) << xsf
                   << std::setw(15) << stat_err
                   << std::setw(15) << syst_err
                   << std::endl;
        }
    }
}

//==============================================================================
// radiative correction for one spectrum                                        
//========================ORIGINAL AUTHORS======================================
// RADCOR FORTRAN CODE                                                          
// ***by Randy Roy Whitney, Phys. Rev. C 9, 2230 - 2235 (1974)                  
//                                                                              
// ***modified by K. Slifer, Temple University, 01/29/03                        
// downloadable at  http://www.jlab.org/~slifer/codes.html                      
//                                                                              
// ***modified by Jaideep Singh, 04/01/07                                       
//                                                                              
//========================REFERENCES============================================
// MOTS69    Radiative Corrections to Elastic and Inelastic ep and mu-p         
//           Scattering                                                         
//           by: L.W. Mo and Y.S. Tsai                                          
//           Reviews of Modern Physics, 41, p205-235, Jan 1969                  
//           http://link.aps.org/abstract/RMP/v41/p205                          
//==============================================================================
// TSAI71    Radiative Corrections to Electron Scattering                       
//           by: Yung-Su Tsai, SLAC PUB 848, Jan 1971                           
//           http://www.slac.stanford.edu/pubs/slacpubs/0750/slac-pub-0848.pdf  
//==============================================================================
// STEIN     Electron scattering at 4° with energies of 4.5-20 GeV              
//           by: S. Stein, W. B. Atwood, E. D. Bloom, R. L. A. Cottrell,        
//               H. DeStaebler, C. L. Jordan§, H. G. Piel, C. Y. Prescott,      
//               R. Siemann, and R. E. Taylor                                   
//           Phys. Rev. D 12, 1884 - 1919 (October 1975)                        
//           http://link.aps.org/abstract/PRD/v12/p1884                         
//==============================================================================
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
    // update these parameters when external RC is ON
    if(external_RC) {
        BTB = set.radl_before*4./3.;
        BTA = set.radl_after*4./3.;
        XIB = set.coll_before;
        XIA = set.coll_after;
    } else {
        BTB = 0;
        BTA = 0;
        XIB = 0;
        XIA = 0;
    }

    Es = set.energy;
    GAMT = gamma(1 + BTB) * gamma(1 + BTA);

    // iteration on data points
    for(auto &point : set.data)
    {

        // kinematics
        Ep = point.Ep;
        double ETAP = (1 - 2*Ep*sin2/target_M);
        R = set.eta/ETAP;

        // calculate low energy corner SIGLOW
        double DHO, BTR;
        internalRC(Es, Ep, DHO, BTR);
        double FBAR = exp(DHO)/GAMT;
        double SIGLOW = FBAR * std::pow(R*delta/Es, BTB+BTR)
                             * std::pow(delta/Ep, BTA+BTR)
                             * (1 - (XIB+XIA)/delta/(1-BTB-BTA-2*BTR));

        // calculate integral along dEs for fixed Ep SIGBEF
        double Esmax = Es - R*delta;
        double Esmin = Ep/ETAP;
        double SIGBEF = 0.;
        if((Esmin <= 0) || (Esmax <= 0) || (Esmin >= Esmax)) {
            std::cout << "Skip point at Ep = " << Ep
                      << ", in spectrum E = " << Es
                      << ", ESMIN = " << Esmin
                      << ", ESMAX = " << Esmax
                      << std::endl;
        } else {
            SIGBEF = simpson(Esmin, Esmax, &CRadCorr::fes, this, simpson_bin_size, n_simpson_bins);
        }

        // calculate integral along dEp for fixed Es SIGAFT
        double Epmax = Es/set.eta;
        double Epmin = Ep+delta;
        double SIGAFT = 0.;
        if((Epmin <= 0) || (Epmax <= 0) || (Epmin >= Epmax)) {
            std::cout << "Skip point at Ep = " << Ep
                      << ", in spectrum E = " << Es
                      << ", EPMIN = " << Esmin
                      << ", EPMAX = " << Esmax
                      << std::endl;
        } else {
            SIGAFT = simpson(Epmin, Epmax, &CRadCorr::fep, this, simpson_bin_size, n_simpson_bins);
        }

        if(radiate) {
            // radiate, update radiated cross section
            point.rad = SIGLOW*set.weight_mott*point.born + SIGBEF + SIGAFT;
        } else {
            // radiative correction, update the born cross section
            point.born = (point.rad - (SIGBEF+SIGAFT)/set.weight_mott)/SIGLOW;
        }
    }
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
//==============================================================================
// TSAI71    Radiative Corrections to Electron Scattering                       
//           by: Yung-Su Tsai, SLAC PUB 848, Jan 1971                           
//           http://www.slac.stanford.edu/pubs/slacpubs/0750/slac-pub-0848.pdf  
//==============================================================================
void CRadCorr::xyrad2d(DataSet &set, bool radiate)
{

}

inline double __rc_phi(double x)
{
    return 1 - x + 3./4.*std::pow(x, 2.);
}

// for integral along dEs
double CRadCorr::fes(const double &Esx)
{
    double DHO, BTR;
    internalRC(Esx, Ep, DHO, BTR);

    double FBAR = exp(DHO)/GAMT;

    // calculate effect of multiple soft photon emission
    double FSOFT = std::pow(1 - Esx/Es, BTB+BTR)*std::pow((Es-Esx)/Ep/R, BTA+BTR);

    double FES = (BTB + BTR)/(Es - Esx)*__rc_phi(1 - Esx/Es);
    FES += XIB/(Es - Esx)/(Es - Esx);
    FES *= FSOFT*FBAR*ftcs(Esx, Ep);
    return FES;
}

// for integral along dEp
double CRadCorr::fep(const double &Epx)
{
    double DHO, BTR;
    internalRC(Es, Epx, DHO, BTR);

    double FBAR = exp(DHO)/GAMT;

    // calculate effect of multiple soft photon emission
    double FSOFT = std::pow(1 - Ep/Epx, BTA+BTR)*std::pow((Epx-Ep)*R/Es, BTB+BTR);

    double FEP = (BTA + BTR)/(Epx - Ep)*__rc_phi(1 - Ep/Epx);
    FEP += XIA/(Epx - Ep)/(Epx - Ep);
    FEP *= FSOFT*FBAR*ftcs(Es, Epx);
    return FEP;
}

// interpolates or extrapolates
double CRadCorr::ftcs(const double &E0, const double &Eb)
{
    if(Eb >= E0/(1+2*E0*sin2/target_M))
        return 0;

    double WEXC = 1 - Eb/E0;

    // less than the lowest energy we have
    if(E0 <= data_sets.at(0).energy) {
        return terp(data_sets.at(0), WEXC)*F_mott/E0/E0;
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
        double FTCS = (terp(data_sets.at(i-1), WEXC)*(E2-E0) + terp(data_sets.at(i), WEXC)*(E0-E1))/(E2-E1);
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

    if(key < array[mid].PA)
        return __rc_binary_search(array, key, first, mid);
    else if(key > array[mid+1].PA)
        return __rc_binary_search(array, key, mid+1, last);
    else
        return mid;
}

double CRadCorr::terp(const DataSet &set, const double &w)
{
    // exceeds the boundary
    if(w < set.data.front().PA)
        return 0.;
    if(w >= set.data.back().PA)
        return set.data.back().born;

    int j = __rc_binary_search(&set.data[0], w, 0, set.data.size());
    if(j < 0)
        return 0;

    // do 2 points interpolation
    if(j <= 1) {
        const DataPoint &p1 = set.data.at(j);
        const DataPoint &p2 = set.data.at(j+1);

        double interp = p1.born*(p2.PA - w) + p2.born*(w - p1.PA);
        // unknown constant 0.00001 here, probably some protection for two same points
        return interp/(p2.PA - p1.PA);// + 0.00001);
    }

    // 3 point parabolic fit
    const DataPoint &p1 = set.data.at(j-1);
    const DataPoint &p2 = set.data.at(j);
    const DataPoint &p3 = set.data.at(j+1);
    double x, xp, xm, a, b, c;
    x = w - p2.PA;
    xp = p3.PA - p2.PA;
    xm = p1.PA - p2.PA;
    a = p2.born;
    c = (p3.born - p1.born)*(xp + xm)/(xp - xm) - (p3.born + p1.born - 2*p2.born);
    c /= xp*xm*2;
    b = (p3.born - p1.born)/(xp - xm) - c*(xp + xm);
    return a + b*x + c*x*x;
}


// do internal RC
// input Es Ep
// output DHO BTR
void CRadCorr::internalRC(const double &Esx, const double &Epx, double &DHO, double &BTR)
{
    double XLQM = log(4.*Esx*Epx*sin2/ELECM/ELECM);
    if(internal_RC) {
        double vertex = 2.*(3./4.*XLQM - 1.); // vertex correction
        double vacuum = 2.*(XLQM/3. - 5./9.); // vacuum correction
        double z0 = -0.5*log(Esx/Epx)*log(Esx/Epx);
        DHO = ALPHA/PI*(vertex + vacuum + Schwinger + z0);
        BTR = ALPHA/PI*(XLQM - 1);
    } else {
        DHO = 0.;
        BTR = 0.;
    }
}

// calculate XI based on STEIN's formula
void CRadCorr::calculateXI(DataSet &set)
{
    // formula from STEIN, but old and probably wrong estimation
    double rad_logp = log(1440*std::pow(target_Z, -2./3.));
    double rad_log = log(183*std::pow(target_Z, -1./3.));
    double xeta = rad_logp/rad_log;
    double xi = (PI*ELECM/2/ALPHA)*(set.radl_before + set.radl_after);
    xi /= (target_Z+xeta)*rad_log;
    set.coll_before = xi/2;
    set.coll_after = xi/2;
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
            DataPoint new_point(nu, cxsn, stat, syst);
            // apply normalization
            new_point.cxsn *= set.normalization;
            new_point.stat *= set.normalization;
            // calculate ep
            new_point.Ep = set.energy - nu;
            new_point.PA = nu/set.energy;
            // initial values for sigrad and sigborn
            new_point.rad = new_point.cxsn/set.weight_mott;
            new_point.born = new_point.rad;

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
            new_set.eta = 1 + 2*energy*sin2/target_M;

            data_sets.push_back(std::move(new_set));

        } else if(c_parser.NbofElements() == 7) { // new data set

            c_parser >> energy >> error >> norm
                     >> radl_bef >> radl_aft >> coll_bef >> coll_aft;

            DataSet new_set(energy, radl_bef, radl_aft, coll_bef, coll_aft, error, norm);
            new_set.weight_mott = F_mott/energy/energy;
            new_set.eta = 1 + 2*energy*sin2/target_M;

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
