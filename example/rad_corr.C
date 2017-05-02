#include "CExpData.h"
#include "CRadCorr.h"

using namespace std;

void rad_corr()
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
    data.ReadConfigFile("configs/data_sets_6deg.conf");

//    rad_cor.Radiate();
    rad_cor.RadiativeCorrection(data);

    data.SaveResult("output/radcor_out.dat");
}


