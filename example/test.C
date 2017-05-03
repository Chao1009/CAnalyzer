#include "CAnalyzer.h"
#include "CMatrix.h"
#include "CEstimator.h"
#include "CExpData.h"
#include "CRadCorr.h"
#include "CNeuralNetwork.h"

using namespace cana;

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

double quad(const double &x)
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
void qe_compare(int energy = 2135, double scale = 1., double shift = 0.)
{
    // plots to show data
    TGraph *g1 = new TGraph();
    TGraph *g2 = new TGraph();
    TGraph *g3 = new TGraph();

    string data_file = "output/radcor_out.dat";

    string calc_file1 = "../../qe_calc/Hannover/s36.dat";
    string calc_file2 = "../../qe_calc/Golak/Golak.dat";


    fill_graph(g1, data_file, 1, 3);

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
    g1->Draw("CP");
    g2->SetLineColor(3);
    g2->SetLineWidth(1);
    g2->SetMarkerStyle(8);
    g2->SetMarkerColor(3);
    g2->SetMarkerSize(0.7);
    g2->Draw("CP");
    g3->SetMarkerStyle(8);
    g3->SetMarkerColor(2);
    g3->SetMarkerSize(0.7);
    g3->SetLineColor(2);
    g3->SetLineWidth(1);
    g3->Draw("CP");
}

