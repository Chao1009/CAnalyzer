#include "CExpData.h"
#include "CRadCorr.h"

using namespace std;

void rad_corr(const string &data_conf = "configs/data_sets_9deg.conf")
{
    CRadCorr rad_cor;
    rad_cor.Configure("configs/rad_corr.conf");

    cout << "Internal RC is ";
    if(rad_cor.GetConfig<bool>("Internal RC"))
        cout << "ON" << endl;
    else
        cout << "OFF" << endl;

    cout << "External RC is ";
    if(rad_cor.GetConfig<bool>("External RC"))
        cout << "ON" << endl;
    else
        cout << "OFF" << endl;

    cout << "User Defined XI is ";
    if(rad_cor.GetConfig<bool>("User Defined XI"))
        cout << "ON" << endl;
    else
        cout << "OFF" << endl;

    cout << "Peaking Approximation is ";
    if(rad_cor.GetConfig<bool>("Energy Peaking Approximation"))
        cout << "ON" << endl;
    else
        cout << "OFF" << endl;

    CExpData data;
    data.ReadConfigFile(data_conf);

//    rad_cor.Radiate(data);
    rad_cor.RadiativeCorrection(data, 1);

    data.SaveResult("output/radcor_out.dat");
}


