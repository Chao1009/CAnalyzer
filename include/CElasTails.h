#ifndef C_ELAS_TAILS_H
#define C_ELAS_TAILS_H

#include "ConfigObject.h"
#include "CExpData.h"
#include "CHe3Elas.h"
#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <iterator>
#include <string>

class CElasTails : public ConfigObject
{
public:
    struct NuPoint
    {
        double nu, tail, weight;

        NuPoint(double n, double t) : nu(n), tail(t), weight(0.)
        {}
    };

    class Acceptance
    {
    public:
        Acceptance();

        void Read(const std::string &path, bool verbose = true);
        double Eval(double pt) const;
        std::string GetDescription() const;

    private:
        // 5 parameters to define rising edge
        double p0_rise, p1_rise, p2_rise, p3_rise, p4_rise;
        // 2 parameters to define flat acceptance
        double p0_main, p1_main;
        // 5 parameters to define falling edge
        double p0_fall, p1_fall, p2_fall, p3_fall, p4_fall;
        // define acceptance range
        double begin_rise, end_rise, begin_fall, end_fall;
    };

public:
    CElasTails(const std::string &path = "");
    virtual ~CElasTails();

    void Configure(const std::string &path);
    CHe3Elas &GetModel() {return he3_model;};
    void Initialize(const CExpData::DataSet &dset);
    void Generate(double nu_beg, double nu_end, double prec = 5e-3);
    void Output(const std::string &path);

private:
    void initGrids(double nu_beg, double nu_end, double prec);
    void setupColl(const std::string &path);
    int calcCollLength(double z, double phi, double &lc);
    void fillData(NuPoint &np, int flag, double xs, double rlcoll, double phi);
    void simElasTails(int flag, double angle, double rloutp);

private:
    CHe3Elas he3_model;
    Acceptance acpt;
    std::vector<NuPoint> points;

    double Es, scat_angle;
    double radl_wall, radl_in, radl_out;
    // target collimator geometry
    // Z-position for upstream and downstream target window (cm)
    double alpha_d, alpha_u, ad_x, ad_y, au_x, au_y, ld, lu;

    // sampling setup
    double zt_min, zt_max, zt_step, ang_range, ang_step, nu_step;
};

std::ostream &operator <<(std::ostream &os, const CElasTails::Acceptance &acpt);

#endif
