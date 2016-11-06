#include "CAnalyzer.h"

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
    estimator.Fit(5, 5);

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
    vector<double> gn = {0.9, 1.0, 1.1};
    for(auto &g : gn)
        cout << setw(10) << "gamma: "
             << setw(8) << g << ", "
             << setw(15) << CRadCorr::gamma(g) << endl;

    vector<double> sn = {1.1, 1.0, 0.8, 0.5, 0.3, 0.0, -0.9, -1.01, -1.1};
    for(auto &s : sn)
        cout << setw(10) << "spence: "
             << setw(8) << s << ", "
             << setw(15) << CRadCorr::spence(s) << endl;

    MyQuad myq(2.0);

    cout << CRadCorr::simpson(0, 10, &quad, 0.01, 1000) << endl;
    cout << CRadCorr::simpson(0, 10, &MyQuad::eval, &myq, 0.01, 1000) << endl;
}

void radcor_test()
{
    CRadCorr rad_cor;
    rad_cor.Configure("radcor.conf");

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

    vector<string> file_list = {"exp_9deg.dat", "model_9deg.dat"};
    rad_cor.ReadExpData(file_list);

    rad_cor.RadiativeCorrection(5);
    rad_cor.SaveResult("radcor_out.dat");
}
