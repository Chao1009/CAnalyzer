#include "CElasTails.h"
#include "CExpData.h"
#include "canalib.h"
#include <string>
#include <iostream>

using namespace std;

void gen_tail(const string &data_conf = "configs/data_sets_9deg.conf")
{
    CExpData data(data_conf);
    CElasTails eltail("configs/elas_tail.conf");

    double theta = data.Angle()*cana::deg2rad;
    double target_M = data.TargetA()*cana::amu;
    for(auto &dset : data.GetSets())
    {
        double sin_scat = sin(theta);
        double sin2 = std::pow(sin(theta/2.), 2);
        double Es = dset.energy;
        double nu_min = int(Es - Es/(1. + 2.*Es*sin2/target_M)) + 1.;
        double nu_max = int(Es) - 1.;

        double nu_beg = cana::clamp(dset.data.front().nu - 10., nu_min, nu_max);
        double nu_end = cana::clamp(dset.data.back().nu + 10., nu_min, nu_max);

        eltail.Initialize(dset);
        eltail.Generate(nu_beg, nu_end, 0.005);
        string output = "output/tails/tail_" + to_string((int)Es) + "_para.out";
        eltail.Output(output);
    }
}

