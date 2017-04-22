#include "CElasTails.h"
#include "canalib.h"
#include <string>
#include <iostream>

using namespace std;

void tail_test(int energy = 1147, double angle = 9.03, double rli = 0.002071, double rlo = 0.06293)
{
    string conf_folder = to_string(int(angle + 0.5)) + "degs/";
    string conf_file = "configs/" + conf_folder + to_string(energy) + ".conf";

    CHe3Elas model;
    model.Configure(conf_file);

    TGraph *g1 = new TGraph();
    TGraph *g2 = new TGraph();
    TGraph *g3 = new TGraph();

    double tarm = 3.015*cana::amu;
    double sin2 = pow(sin(angle*cana::deg2rad/2.), 2);
    double Ep_el = energy/(1. + 2.*energy*sin2/tarm);
    double Ep_step = 1.0;
    double Ep_max = int(Ep_el - 0.01*energy);

    for(double Ep = Ep_max; Ep > energy*0.3; Ep -= Ep_step)
    {
        g1->SetPoint(g1->GetN(), energy - Ep, 1000.*model.GetRadXS(energy, Ep, angle*cana::deg2rad, rli, rlo));
    }

    model.SetConfigValue("Use X. Yan's Formula", false);
    model.Configure();

    for(double Ep = Ep_max; Ep > energy*0.3; Ep -= Ep_step)
    {
        g2->SetPoint(g2->GetN(), energy - Ep, 1000.*model.GetRadXS(energy, Ep, angle*cana::deg2rad, rli, rlo));
    }

    TCanvas *c1 = new TCanvas("elastic tail","elastic tail",200,10,700,800);
    c1->SetGrid();
    c1->Divide(1, 2);
    c1->cd(1);

    // plot data
    g1->SetLineColor(1);
    g1->SetLineWidth(1);
    g1->SetMarkerStyle(8);
    g1->SetMarkerColor(1);
    g1->SetMarkerSize(0.7);
    g1->Draw("AL");
    g2->SetMarkerStyle(8);
    g2->SetMarkerColor(2);
    g2->SetMarkerSize(0.7);
    g2->SetLineColor(2);
    g2->SetLineWidth(1);
    g2->Draw("L");

    Double_t *y1 = g1->GetY();
    Double_t *y2 = g2->GetY();
    for(int i = 0; i < g1->GetN() && i < g2->GetN(); ++i)
    {
        if(y1[i] != 0.)
            g3->SetPoint(i, energy - Ep_max + i*Ep_step, (y1[i] - y2[i])/y1[i]*100.);
    }

    c1->cd(2);
    g3->Draw("ALP");

}

void gen_tail()
{
    vector<int> energies = {1147, 2234, 3319, 3775, 4404,
                            2135, 2845, 4209};
    vector<string> conf_folder = {"9degs/", "9degs/", "9degs/", "9degs/", "9degs/",
                                  "6degs/", "6degs/", "6degs/"};
    for(size_t i = 0; i < 1 /*energies.size()*/; ++i)
    {
        int energy = energies.at(i);
        string folder = conf_folder.at(i);
        string type;

        CElasTails eltail;
        // read configuration
        eltail.Configure("configs/" + folder + to_string(energy) + ".conf");

        // print out information
        cout << "Generating unpolarized tail for energy = " << energy << endl;
        bool pol = eltail.GetConfig<bool>("Polarized");
        if(pol) {
            int angle = eltail.GetConfig<int>("Polarization Angle");
            cout << "Polarized, polarization angle is "
                 << angle
                 << endl;
            if(angle == 180)
                type = "para";
            else if(angle == 270)
                type = "perp";
            else
                type = to_string(angle);
        } else {
            cout << "Unpolarized." << endl;
            type = "unpl";
        }

        eltail.Generate();

        string out_path = "output/" + to_string(energy) + "_complete_" + type;
        string cut_file = eltail.GetConfig<string>("Acceptance File");
        auto numbers = ConfigParser::find_integers(cut_file);
        if(numbers.size())
            out_path += "_" + to_string(numbers.back());

        eltail.Output(out_path + ".out");
    }
}

