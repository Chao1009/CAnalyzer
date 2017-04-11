#ifndef C_ELAS_TAILS_H
#define C_ELAS_TAILS_H

#include "ConfigObject.h"
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
    struct DataPoint
    {
        double phi, tail, weight;

        DataPoint() {};
        DataPoint(const double &p) : phi(p) {};
        DataPoint(const double &p, const double &t, const double &w)
        : phi(p), tail(t), weight(w)
        {};
    };
    typedef std::map<double, std::vector<DataPoint>> DataSet;

    class Acceptance
    {
    public:
        Acceptance();

        void Read(const std::string &path, bool verbose = true);
        double Eval(const double &pt) const;
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
    CElasTails();
    virtual ~CElasTails();

    void Configure(const std::string &path = "");
    CHe3Elas &GetModel() {return he3_model;};
    void Generate();
    void Output(const std::string &path);

private:
    void setupColl(const std::string &path);
    int calcCollLength(const double &z, const double &phi, double &lc);
    void fillData(const int &flag,
                  const double &nu,
                  const double &tail,
                  const double &rlcoll,
                  const double &phi);
    void rtails(int flag, double angle, double rloutp);
    DataPoint punchThrough(const double &nu,
                           const double &tail,
                           const double &phi,
                           const double &rl);
private:
    CHe3Elas he3_model;
    Acceptance acpt;
    DataSet dset;

    double in_energy, scat_angle, nu_min, nu_max, deltas, deltap, Mtarg;
    double radl_wall, radl_in, radl_out;
    bool col_loss;
    // target collimator geometry
    // Z-position for upstream and downstream target window (cm)
    double alpha_d, alpha_u, ad_x, ad_y, au_x, au_y, ld, lu;

    // sampling setup
    double zt_min, zt_max, zt_step, ang_range, ang_step;
    double step, fine_step, finer_step, fine_range, finer_range;
};

std::ostream &operator <<(std::ostream &os, const CElasTails::Acceptance &acpt);

#endif