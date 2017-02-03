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

void model_wrapper_test()
{
    double Z = 2.0, A = 3.0, Q2 = 1.0, W2 = 0.5;
    double f1, f2, rc;
    Bosted_f1f2in09(Z, A, Q2, W2, &f1, &f2, &rc);
    cout << f1 << ", " << f2 << endl;
    Bosted_f1f2qe09(Z, A, Q2, W2, &f1, &f2);
    cout << f1 << ", " << f2 << endl;

    double xs;
    for(double nu = 521.0; nu <= 550.0; nu += 1.0)
    {
        Bosted_xs(Z, A, 0.6, (600 - nu)/1000., 5.99/RADDEG, &xs);
        cout << nu << ", " << 1000.*xs << endl;
    }
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

    cout << "Peaking Approximation is ";
    if(rad_cor.GetConfig<bool>("Peaking Approximation"))
        cout << "ON" << endl;
    else
        cout << "OFF" << endl;

//    rad_cor.Radiate();
    rad_cor.RadiativeCorrection(1);
    rad_cor.SaveResult("radcor_out.dat");
}
