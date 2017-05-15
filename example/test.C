#include "CAnalyzer.h"
#include "CMatrix.h"
#include "CEstimator.h"
#include "CExpData.h"
#include "CRadCorr.h"
#include "CElasTails.h"
#include "CNeuralNetwork.h"

void fit_test()
{
    // read formula
    CEstimator estimator("fit_par_1147.dat");

    // read data
    CAnalyzer c_ana;
    c_ana.ReadData("1147.dat");
    auto nu = c_ana.GetColumn(0);
    auto cxsn = c_ana.GetColumn(1);
    auto stat = c_ana.GetColumn(2);

    // set data to the estimator
    estimator.SetDataPoints(nu, cxsn, stat);

    // plots to show data
    TGraphErrors *g1 = new TGraphErrors();
    TGraph *g2 = new TGraph();
    TGraph *g3 = new TGraph();

    // fill the function value before fit
    for(size_t i = 0; i < cxsn.size(); ++i)
    {
        g2->SetPoint(g2->GetN(), nu.at(i), estimator.GetFormulaVal(nu.at(i)));
    }

    // fill data points
    for(size_t i = 0; i < cxsn.size(); ++i)
    {
        auto pN = g1->GetN();
        g1->SetPoint(pN, nu.at(i), cxsn.at(i));
        g1->SetPointError(pN, 0, stat.at(i));
    }

    // fit data
    estimator.Fit();

    // save the optimized parameters
    estimator.SaveFormula("fit_par_1147.dat");

    // fill the function value after fit
    for(size_t i = 0; i < cxsn.size(); ++i)
        g3->SetPoint(g3->GetN(), nu.at(i), estimator.GetFormulaVal(nu.at(i)));

    // create canvas
    TCanvas *c1 = new TCanvas("saGDH fit","fitting band",200,10,700,500);

    // set frame range
    c1->SetGrid();
    c1->DrawFrame(20, 0, 700, 1600);

    // plot data
    g1->SetMarkerStyle(1);
    g1->SetMarkerColor(5);
    g1->SetMarkerSize(0.7);
    g1->Draw("P");
    // plot function before fitting
    g2->SetLineColor(2);
    g2->SetLineWidth(1);
    g2->Draw("C");
    /*
    // plot function after fitting
    g3->SetLineColor(4);
    g3->SetLineWidth(1);
    g3->Draw("C");
    */

}

void matrix_test()
{
    CMatrix m(6);
    m = {2, 8, -1, 1, 2, 2,
         4, 2, 1, 6, 4, 5,
         9, 1, 2, 7, 9, 7,
         1, 2, 5, 2, 8, 6,
         1, 2, 3, 4, 5, 6,
         6, 5, 4, 3, 2, 1};
    CMatrix l = m;

//    cout << m << endl;
//    cout << l << endl;
    cout << m.Det_Leibniz() << endl;
    cout << m.Det_Laplace() << endl;
/*
    cout << m.UpperLeft(1) << endl;
    cout << m.Det(1) << endl;
    cout << m.UpperLeft(2) << endl;
    cout << m.Det(2) << endl;
    cout << m.UpperLeft(3) << endl;
    cout << m.Det(3) << endl;
    cout << m.UpperLeft(4) << endl;
    cout << m.Det(4) << endl;
    auto n = cholesky(m);
    cout << n << endl;
    cout << n * transpose(n) << endl;
    cout << tril(m) << endl;
    cout << triu(m) << endl;
    cout << power(m, 0) << endl;
    cout << power(m, 1) << endl;
    cout << power(m, 10) << endl;
    cout << m*m*m << endl;
    cout << transpose(m) << endl;
*/
}

void cana_test()
{
    cout << cana::double2int(-123.6)
         << ", "
         << cana::double2int(-12345678.9)
         << endl;
}

double quad(double x)
{
    return x*x;
}

class MyQuad
{
    double c;
public:
    MyQuad(double cc = 1) : c(cc) {};
    double eval(const double &x) {return c*x*x;};
};

void function_test()
{
    vector<int> test = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    vector<double> search = {0.5, 1.5, 8.5, 0, 9, 9.5, -0.5};
    for(auto &val : search)
    {
        auto it_pair = cana::binary_search_interval(test.begin(), test.end(), val);
        if(it_pair.first == test.end() || it_pair.second == test.end())
            cout << val << " not found." << endl;
        else
            cout << val << " is in " << *it_pair.first << ", " << *it_pair.second << endl;
    }

    vector<double> gn = {0.9, 1.0, 1.1};
    for(auto &g : gn)
        cout << setw(10) << "gamma: "
             << setw(8) << g << ", "
             << setw(15) << cana::gamma(g) << endl;

    vector<double> sn = {1.1, 1.0, 0.8, 0.5, 0.3, 0.0, -0.9, -1.01, -1.1};
    for(auto &s : sn)
        cout << setw(10) << "spence: "
             << setw(8) << s << ", "
             << setw(15) << cana::spence(s) << endl;

    MyQuad myq(2.0);

    cout << cana::simpson(0, 10, &quad, 0.01, 1000) << endl;
    cout << cana::simpson(0, 10, &MyQuad::eval, &myq, 0.01, 1000) << endl;
    auto lmd_quad = [](double x) {return x*x;};
    cout << cana::simpson(0, 10, 1000, lmd_quad) << endl;
}

void model_wrapper_test(double energy, double angle, double nu_min, double nu_max)
{
    double Z = 2.0, A = 3.0, Q2 = 1.0, W2 = 0.5;

    CModelWrapper wrapper;
    wrapper.SetRange(6000, 50000, 500, 900, 1900, 1000);
    ofstream outf("test_model.dat");
    double xs1, xs2, xs3;
    for(double nu = nu_min; nu <= nu_max; nu += 1.0)
    {
        Bosted_xs(Z, A, energy, energy - nu, angle*cana::deg2rad, &xs1);
        QFS_xs(Z, A, energy, energy - nu,angle*cana::deg2rad, &xs2);
        xs3 = wrapper.GetCrossSection(energy, energy - nu, angle);
        outf << setw(8) << nu
             << setw(20) << xs1
             << setw(20) << xs2
             << setw(20) << xs3
             << endl;
    }
}

void init_model_test()
{
    CRadCorr rad_cor;
    rad_cor.Configure("configs/rad_corr.conf");

    CExpData data;
    data.ReadConfigFile("configs/data_sets_6deg.conf");

    rad_cor.Initialize(data);

    for(int nu = 10; nu < 200; ++nu)
    {
        cout << nu << ", "
             << rad_cor.GetModel().GetCrossSection(2135, 2135 - nu, 6.10*cana::deg2rad)
             << endl;
    }
}

void interp_test(double energy = 1000, double nu_beg = 10, double nu_end = 100)
{
    CExpData data;
    data.ReadConfigFile("data_sets_9deg.conf");
    for(int nu = 10; nu < 100; ++nu)
    {
        cout << nu << ", " << data.GetCrossSection(energy, energy - nu) << endl;
    }
}

#include "../../sys_study/scripts/utils.C"

void fill_graph_radcor(TGraph *gr,
                       const string &rad_file,
                       int nu_col1,
                       int xs_col1,
                       int born_col,
                       const string &data_file,
                       int nu_col2,
                       int xs_col2)
{
    CAnalyzer cana;
    cana.ReadData(rad_file);
    auto nu1 = cana.GetColumn(nu_col1);
    auto xs1 = cana.GetColumn(xs_col1);
    auto born = cana.GetColumn(born_col);

    cana.ReadData(data_file);
    auto nu2 = cana.GetColumn(nu_col2);
    auto xs2 = cana.GetColumn(xs_col2);

    for(size_t i = 0; i < nu2.size(); ++i)
    {
        for(size_t j = 0; j < nu1.size(); ++j)
        {
            if(nu2.at(i) == nu1.at(j)) {
                double cxsn = born.at(j)/xs1.at(j)*xs2.at(i);
                gr->SetPoint(gr->GetN(), nu2.at(i), cxsn);
            }
        }
    }
}

void qe_compare(int energy = 2135, double scale = 1., double shift = 0.)
{
    // plots to show data
    TGraph *g1a = new TGraph();
    TGraph *g1b = new TGraph();
    TGraph *g1c = new TGraph();
    TGraph *g2 = new TGraph();
    TGraph *g3 = new TGraph();

    string rad_file1 = "output/radcor_" + to_string(energy) + "_0mm.dat";
    string rad_file2 = "output/radcor_" + to_string(energy) + "_2mm.dat";
    string rad_file3 = "output/radcor_" + to_string(energy) + "_4mm.dat";
    string data_file1 = "data/" + to_string(energy) + "_tailsub_0mm.dat";
    string data_file2 = "data/" + to_string(energy) + "_tailsub_2mm.dat";
    string data_file3 = "data/" + to_string(energy) + "_tailsub_4mm.dat";

    string calc_file1 = "../../qe_calc/Hannover/s36.dat";
    string calc_file2 = "../../qe_calc/Golak/Golak.dat";


    fill_graph_radcor(g1a, rad_file1, 1, 2, 3, data_file1, 0, 1);
    fill_graph_radcor(g1b, rad_file2, 1, 2, 3, data_file2, 0, 1);
    fill_graph_radcor(g1c, rad_file3, 1, 2, 3, data_file3, 0, 1);

    // hannover is calculated at 6 and 9 degree, need to be corrected
    vector<int> sets_9deg = {1147, 2234, 3319, 3775, 4404};
    double h_scale = 1.;
    if(cana::is_in(energy, sets_9deg.begin(), sets_9deg.end()))
    {
        h_scale = scale_mott(9.0, 9.03);
    }
    else
    {
        h_scale = scale_mott(6.0, 6.10);
    }

    fill_graph(g2, calc_file1, 1, 2, energy, h_scale*scale, shift);
    fill_graph(g3, calc_file2, 1, 2, energy, scale, shift);

    TCanvas *c1 = new TCanvas("unpol tail","unpol tail",200,10,700,500);
    c1->SetGrid();

    double x1, x2, y1, y2;
    get_frame_range(x1, y1, x2, y2, {g2, g3});
    c1->DrawFrame(x1, y1, x2, y2);

    // plot data
    g1a->SetMarkerStyle(20);
    g1a->SetMarkerSize(0.7);
    g1a->SetMarkerColor(4);
    g1b->SetMarkerStyle(21);
    g1b->SetMarkerSize(0.7);
    g1b->SetMarkerColor(1);
    g1c->SetMarkerStyle(22);
    g1c->SetMarkerSize(0.7);
    g1c->SetMarkerColor(2);
    g1a->Draw("P");
    g1b->Draw("P");
    g1c->Draw("P");
    g2->SetLineColor(8);
    g2->SetLineWidth(1);
    g2->Draw("CP");
    g3->SetLineColor(9);
    g3->SetLineWidth(1);
    g3->Draw("CP");
}

void tail_compare(int energy = 1147)
{
    // plots to show data
    TGraphErrors *g1 = new TGraphErrors();
    TGraph *g2a = new TGraph();
    TGraph *g2b = new TGraph();
    TGraph *g2c = new TGraph();

    string data_file = "../../data/cxsn_raw/He3_bc_cxsn_" + to_string(energy) + "_n2sub.txt";

    // black
    string tail_file1 = "output/tail_" + to_string(energy) + "_4mm_prec.dat";
    // red
    string tail_file2 = "output/tail_" + to_string(energy) + "_4mm_polsig.dat";
    // blue
    string tail_file3 = "output/tail_" + to_string(energy) + "_4mm_xy.dat";

    CAnalyzer c_ana;
    c_ana.ReadData("../../sys_study/tail_shifts.txt");
    auto energies = c_ana.GetColumn(0);
    auto shifts = c_ana.GetColumn(1);
    double shift = 0.;
    for(size_t i = 0; i < shifts.size(); ++i)
    {
        if(int(energies.at(i) + 0.5) == energy)
            shift = shifts.at(i);
    }

    fill_graph_error(g1, data_file, 2, 3, 4);
    fill_graph(g2a, tail_file1, 1, 2, 1.0, shift);
    fill_graph(g2b, tail_file2, 1, 2, 1.0, shift);
    fill_graph(g2c, tail_file3, 1, 2, 1.0, shift);

    TCanvas *c1 = new TCanvas("unpol tail","unpol tail",200,10,700,500);
    c1->SetGrid();

    double x1, x2, y1, y2;
    get_frame_range(x1, y1, x2, y2, {g1});
    c1->DrawFrame(x1, y1, x2, y2);

    // plot data
    g2a->SetLineColor(1);
    g2b->SetLineColor(2);
    g2c->SetLineColor(4);
    g2a->Draw("C");
    g2b->Draw("C");
    g2c->Draw("C");
    g1->SetLineColor(9);
    g1->SetMarkerColor(9);
    g1->SetMarkerStyle(8);
    g1->SetMarkerSize(0.7);
    g1->Draw("P");
}

void subtail_compare(int energy = 2135)
{
    // plots to show data
    TGraph *g1a = new TGraph();
    TGraph *g1b = new TGraph();
    TGraph *g1c = new TGraph();

    string cxsn_file1 = "data/" + to_string(energy) + "_tailsub_0mm.dat";
    string cxsn_file2 = "data/" + to_string(energy) + "_tailsub_2mm.dat";
    string cxsn_file3 = "data/" + to_string(energy) + "_tailsub_4mm.dat";

    fill_graph(g1a, cxsn_file1, 0, 1);
    fill_graph(g1b, cxsn_file2, 0, 1);
    fill_graph(g1c, cxsn_file3, 0, 1);

    TCanvas *c1 = new TCanvas("unpol tail","unpol tail",200,10,700,500);
    c1->SetGrid();

    double x1, x2, y1, y2;
    get_frame_range(x1, y1, x2, y2, {g1a, g1b, g1c});
    c1->DrawFrame(x1, y1, x2, y2);

    // plot data
    g1a->SetMarkerStyle(20);
    g1a->SetMarkerSize(0.7);
    g1a->SetMarkerColor(4);
    g1b->SetMarkerStyle(21);
    g1b->SetMarkerSize(0.7);
    g1b->SetMarkerColor(1);
    g1c->SetMarkerStyle(22);
    g1c->SetMarkerSize(0.7);
    g1c->SetMarkerColor(2);
    g1a->Draw("P");
    g1b->Draw("P");
    g1c->Draw("P");
}

void elas_test(double energy = 1147, double angle = 9.03, double nu_beg = 10, double nu_end = 700)
{
    TGraph *g1 = new TGraph();
    TGraph *g2 = new TGraph();

    CHe3Elas he3_model;
    he3_model.Configure("configs/elas_tail.conf");

    he3_model.SetConfigValue("Use X. Yan's Formula", true);
    he3_model.Configure();
    for(double nu = nu_beg; nu <= nu_end; nu += 1.0)
    {
        g1->SetPoint(g1->GetN(), nu, 1000.*he3_model.GetRadXS(energy, energy - nu, angle*cana::deg2rad, 0., 0.));
    }

    he3_model.SetConfigValue("Use X. Yan's Formula", false);
    he3_model.Configure();
    for(double nu = nu_beg; nu <= nu_end; nu += 1.0)
    {
        g2->SetPoint(g2->GetN(), nu, 1000.*he3_model.GetRadXS(energy, energy - nu, angle*cana::deg2rad, 0., 0.));
    }

    TCanvas *c1 = new TCanvas("unpol tail","unpol tail",200,10,700,500);
    c1->SetGrid();

    g1->SetLineColor(2);
    g1->Draw("AC");
    g2->Draw("C");
}

