//============================================================================//
// A C++ version of the punchcol_url.f rtails.f and mulsc_mod.f               //
// Generate the elastic tail for Helium-3 with radiation effect and punch     //
// through effect from the target collimators                                 //
//                                                                            //
// Original Authors: A. Deur, K. Sliffer, J. Singh, V. Sulkosky               //
// Adapted to C++: C. Peng, Feb. 13, 2017                                     //
//============================================================================//

#include "CElasTails.h"
#include "canalib.h"



const static double __gauss[] = {0.1306, 0.1169, 0.938, 0.674, 0.434, 0.250,
                                 0.129, 0.0598, 0.0248, 0.00921, 0.00306, 0.000912};

CElasTails::CElasTails(const std::string &path)
{
    if(!path.empty())
        Configure(path);
}

CElasTails::~CElasTails()
{
    // place holder
}

void CElasTails::Configure(const std::string &path)
{
    ConfigObject::Configure(path);

    ang_range = getDefConfig<double>("Angle Range", 6.0);
    ang_step = getDefConfig<double>("Angle Step", 0.06);
    zt_min = getDefConfig<double>("Target Z Min", -20);
    zt_max = getDefConfig<double>("Target Z Max", 20);
    zt_step = getDefConfig<double>("Target Z Step", 0.2);
    step = getDefConfig<double>("Step", 10.0);
    fine_step = getDefConfig<double>("Fine Step", 1.0);
    finer_step = getDefConfig<double>("Finer Step", 0.1);
    fine_range = getDefConfig<double>("Fine Range", 50.0);
    finer_range = getDefConfig<double>("Finer Range", 15.0);

    // pass this file to the model
    he3_model.Configure(path);
}

void CElasTails::setupColl(const std::string &path)
{
    ConfigParser c_parser;
    if(!c_parser.ReadFile(path)) {
        std::cerr << "Cannot open file "
                  << "\"" << path << "\""
                  << std::endl;
        return;
    }

    std::cout << "Setting up collimators from "
              << "\"" << path << "\"."
              << std::endl;

    c_parser.ParseAll();

    c_parser >> alpha_d >> alpha_u
             >> ad_x >> ad_y >> ld
             >> au_x >> au_y >> lu;

    alpha_d /= cana::rad2deg;
    alpha_u /= cana::rad2deg;
    return;
}

void CElasTails::Initialize(const CExpData &data, int i)
{
    // hard coded value for wall thickness
    radl_wall = 0.06;
    scat_angle = data.Angle();
    in_energy = data.GetSet(i).energy;

    double sin_scat = sin(scat_angle*cana::deg2rad);
    double sin2 = std::pow(sin(scat_angle*cana::deg2rad/2.), 2);
    radl_in = data.GetSet(i).radl_before;
    radl_out = data.GetSet(i).radl_after - radl_wall/10.61/sin_scat;

    setupColl(data.GetSet(i).coll_file);
    acpt.Read(data.GetSet(i).accpt_file);

    nu_max = data.GetSet(i).data.back().nu + 10.;
    double target_M = data.TargetA()*cana::amu;
    nu_elas = int(in_energy - in_energy/(1. + 2.*in_energy*sin2/target_M)) + 1.;
    nu_min = nu_elas;
}

void CElasTails::Generate(double nu_beg, double nu_end)
{
    if(nu_beg > 0.) nu_min = nu_beg;
    if(nu_end > 0.) nu_max = nu_end;

    // angle range
    double ang_min = scat_angle - ang_range/2.0, ang_max = scat_angle + ang_range/2.0;

    int zn = int((zt_max - zt_min)/zt_step + 0.5) + 1;
    int an = int((ang_max - ang_min)/ang_step + 0.5) + 1;

    int flag;
    double lc;
    int count = 0;
    for(int i = 0; i < zn; ++i)
    {
        double z = zt_min + i*zt_step;
        for(int j = 0; j < an; ++j)
        {
            std::cout << "Sampling... "
                      << ++count << "/" << zn*an
                      << "\r" << std::flush;

            double ang = ang_min + j*ang_step;
            flag = calcCollLength(z, ang/cana::rad2deg, lc);

            // blocked by collimators
            if(flag == 4 || flag == 7)
            {
                continue;
            }

            simElasTails(flag, ang, lc);
        }
    }
    std::cout << std::endl;
}

int CElasTails::calcCollLength(const double &z, const double &phi, double &lc)
{
    double bux, buy, bdx, bdy, cux, cuy, cdx, cdy;
    // where the track intersect the collimator lines
    bux = (z*tan(phi) - au_x*tan(alpha_u) + au_y)/(tan(phi) - tan(alpha_u));
    buy = (bux - au_x)*tan(alpha_u) + au_y;
    bdx = (z*tan(phi) - ad_x*tan(alpha_d) + ad_y)/(tan(phi) - tan(alpha_d));
    bdy = (bdx - ad_x)*tan(alpha_d) + ad_y;

    cux = (z*tan(phi) - au_x*tan(alpha_u + cana::pi/2.) + au_y)
          / (tan(phi) - tan(alpha_u + cana::pi/2.));
    cuy = (cux - au_x)*tan(alpha_u + cana::pi/2.) + au_y;
    cdx = (z*tan(phi) - ad_x*tan(alpha_d + cana::pi/2.) + ad_y)
          / (tan(phi) - tan(alpha_d + cana::pi/2.));
    cdy = (cdx - ad_x)*tan(alpha_d + cana::pi/2.) + ad_y;

    // go through upper coll.
    if((bux >= au_x) && (bux <= au_x + lu*cos(alpha_u)) &&
       (buy >= au_y) && (buy <= au_y + lu*sin(alpha_u)))
    {
        if(phi >= alpha_u)
        {
            // length of collimator crossed
            lc = (lu - std::fabs(au_x - bux)/cos(alpha_u))*cos(phi - alpha_u);
            return 1;
        }
        else
        {
            lc = std::fabs(au_x - bux)/cos(alpha_u)*cos(phi - alpha_u);
            return 2;
        }
    }
    // stopped by upper coll.
    else if((cuy >= au_y) && (cux <= au_x))
    {
        return 4;
    }
    // go through down coll.
    else if((bdx >= ad_x) && (bdx <= ad_x + ld*cos(alpha_d)) &&
            (bdy >= ad_y) && (bdy <= ad_y + ld*sin(alpha_d)))
    {
        if(phi >= alpha_d)
        {
            lc = std::fabs(ad_x - bdx)/cos(alpha_d)*cos(phi - alpha_d);
            return 5;
        }
        else
        {
            lc = (ld - std::fabs(ad_x - bdx)/cos(alpha_d))*cos(phi - alpha_d);
            return 6;
        }
    }
    // stopped by the down coll.
    else if((cdy <= ad_y) && (cdx >= ad_x))
    {
        return 7;
    }
    // did not touch any coll.
    else
    {
        lc = 0.;
        return 3;
    }
}

void CElasTails::simElasTails(int flag, double angle, double lc)
{
    // kinematics
    double theta = angle*cana::deg2rad;

    // compute radiation length with collimator
    double rloutp = radl_out + radl_wall/10.61/sin(theta) + lc/0.35;

    double nu_in = nu_min;
    while(nu_in <= nu_max)
    {
        // convert ub/MeV/sr to nb/MeV/sr
        double sigrad = 1000.*he3_model.GetRadXS(in_energy, in_energy - nu_in, theta, radl_in, rloutp);
        if(sigrad > 0.)
            fillData(flag, nu_in, sigrad, lc/0.35, angle);

        // set nu_in for the next step
        if(nu_in + 0.9*finer_step <= finer_range + nu_elas)
        {
            nu_in += finer_step;
        }
        else if(nu_in + 0.9*fine_step <= fine_range + nu_elas)
        {
            nu_in += fine_step;
        }
        else
        {
            nu_in += step;
        }
    }
}

// it is a 2d integral int_es(int_ep)
inline CElasTails::TailPoint CElasTails::punchThrough(const double &nu,
                                                      const double &tail,
                                                      const double &phi,
                                                      const double &rl)
{
    TailPoint point(phi, 0., 0.);
    double phi2, wgt;
    double ramoli = 13.6/(in_energy - nu/2.)*74.*sqrt(rl)*(1. + 0.038*log(rl));

    for(int j = 0; j < 12; ++j)
    {
        // first half
        phi2 = phi + ramoli*(j+1)/3.;
        wgt = acpt.Eval(phi2);
        point.tail += tail*wgt*__gauss[j];
        point.weight += wgt;

        // second half
        phi2 = phi - ramoli*(j+1)/3.;
        wgt = acpt.Eval(phi2);
        point.tail += tail*wgt*__gauss[j];
        point.weight += wgt;

    }

    point.tail /= 24.;
    point.weight /= 24.;

    return point;
}

void CElasTails::fillData(const int &flag,
                          const double &nu,
                          const double &tail,
                          const double &rlcoll,
                          const double &phi)
{
    TailPoint point;

    // no punch through effect
    if(flag == 3)
    {
        point.phi = phi;
        point.weight = acpt.Eval(phi);
        point.tail = point.weight*tail;
    }
    else
    {
        point = punchThrough(nu, tail, phi, rlcoll);
    }

    // no need to save this point
    if(point.weight == 0.)
        return;

    auto it = tset.find(nu);
    if(it == tset.end())
    {
        std::vector<TailPoint> points;
        points.reserve(5000);
        points.push_back(point);
        tset[nu] = points;
    }
    else
    {
        it->second.push_back(point);
    }
}

double __average_tail(const std::vector<CElasTails::TailPoint> &data)
{
    double total_tail = 0., total_weight = 0.;
    for(auto &val : data)
    {
        total_tail += val.tail;
        total_weight += val.weight;
    }

    return total_tail/total_weight;
}

void CElasTails::Output(const std::string &path)
{
    std::cout << "Save result to " << path << std::endl;
    std::ofstream outfile(path);

    for(auto &it : tset)
    {
        outfile << std::setprecision(12);
        outfile << std::setw(8) << in_energy
                << std::setw(16) << it.first
                << std::setw(20) << __average_tail(it.second)
                << std::endl;
    }

}


// initialized constructor
CElasTails::Acceptance::Acceptance()
: p0_rise(0.), p1_rise(0.), p2_rise(0.), p3_rise(0.), p4_rise(0.),
  p0_main(0.), p1_main(0.),
  p0_fall(0.), p1_fall(0.), p2_fall(0.), p3_fall(0.), p4_fall(0.),
  begin_rise(0.), end_rise(0.), begin_fall(0.), end_fall(0.)
{
    // place holder
};

// get acceptance
double CElasTails::Acceptance::Eval(const double &pt)
const
{
    // rising edge
    if ((pt > begin_rise) && (pt < end_rise))
    {
        double pt2 = pt*pt;
        return p0_rise + p1_rise*pt + p2_rise*pt2 + p3_rise*pt*pt2 +
               p4_rise*pt2*pt2;
    }
    // flat range
    else if ((pt >= end_rise) && (pt <= begin_fall))
    {
        return p0_main + p1_main*pt;
    }
    // falling edge
    else if ((pt > begin_fall) && (pt < end_fall))
    {
        double pt2 = pt*pt;
        return p0_fall + p1_fall*pt + p2_fall*pt2 + p3_fall*pt*pt2 +
               p4_fall*pt2*pt2;
    }
    // not in acceptance
    else
    {
        return 0.;
    }
}

void CElasTails::Acceptance::Read(const std::string &path, bool verbose)
{
    ConfigParser c_parser;
    c_parser.ReadFile(path);

    // first line
    c_parser.ParseAll();

    c_parser >> p0_rise >> p1_rise >> p2_rise >> p3_rise >> p4_rise
             >> p0_main >> p1_main
             >> p0_fall >> p1_fall >> p2_fall >> p3_fall >> p4_fall
             >> begin_rise >> end_rise >> begin_fall >> end_fall;

    if(verbose)
    {
        std::cout << "Read acceptance from " << "\"" << path << "\"." << std::endl;
        std::cout << GetDescription() << std::endl;
    }
}

#include <initializer_list>
std::string __polynomial(const std::initializer_list<double> &para)
{
    std::string res;
    int order = 0;

    for(double ele : para)
    {
        if (ele == 0.)
        {
            ++order;
            continue;
        }

        if(order > 0)
        {
            if(ele < 0)
                res += " - ";
            else
                res += " + ";
            res += std::to_string(fabs(ele)) + "*x^" + std::to_string(order);
        }
        else
        {
            res += std::to_string(ele);
        }
        ++order;
    }

    return res;
}

std::string CElasTails::Acceptance::GetDescription()
const
{
    std::string res;
    res =   "Rising edge from " + std::to_string(begin_rise)
          + " to " + std::to_string(end_rise)
          + ":\n\t"
          + __polynomial({p0_rise, p1_rise, p2_rise, p3_rise, p4_rise}) + "\n"
          + "Flat region from " + std::to_string(end_rise)
          + " to " + std::to_string(begin_fall)
          + ":\n\t"
          + __polynomial({p0_main, p1_main}) + "\n"
          + "Falling edge from " + std::to_string(begin_fall)
          + " to " + std::to_string(end_fall)
          + ":\n\t"
          + __polynomial({p0_fall, p1_fall, p2_fall, p3_fall, p4_fall});

    return res;
}
