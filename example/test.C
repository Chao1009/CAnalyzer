#include "CAnalyzer.h"
#include "CMatrix.h"
#include "CEstimator.h"
#include "CExpData.h"
#include "CRadCorr.h"
#include "CElasTails.h"
#include "CNeuralNetwork.h"
#include "CMaidModel.h"
#include "../../sys_study/scripts/utils.C"


// combine polarized properties for helium-3
double pol_he3(double valp, double valn)
{
    return valp*(2.*-0.027*valp - 0.014) + valn*(1.*0.87 + 0.056);
}

void mm_test(double energy = 1147, double angle = 9.03)
{
    CMaidModel maid("maid_g1g2.dat");

    CMaidModel::MaidValue pro, neu;
    // convert MeV^3 to nb*MeV
    double unit = cana::hbarc2*1e7;

    TGraph *gr1 = new TGraph();
    for(double nu = 100.; nu < energy*0.8; nu += 1.0)
    {
        double w = get_W(energy, nu, angle);
        double q2 = get_Q2(energy, nu, angle);
        if(w < 1080. || w > 2000.)
            continue;

        maid.Interp(q2, w, pro, neu);

        double y = nu/energy;
        double c = (cana::amu*3.01493*q2)/4./cana::alpha/cana::alpha*y/(1 - y)/(2. - y);
        double coef1 = tan(angle*cana::deg2rad/2.);
        double coef2 = (1. + (1. - y)*cos(angle*cana::deg2rad))/(1. - y)/sin(angle*cana::deg2rad);
        double g1 = pol_he3(pro.g1, neu.g1);
        double g2 = pol_he3(pro.g2, neu.g2);

        double sig_t = (g1 + 2./y*g2)/c/(coef1 + coef2)*unit;
        gr1->SetPoint(gr1->GetN(), nu, sig_t);
    }

    TCanvas *c1 = new TCanvas("saGDH fit","fitting band",200,10,700,500);
    gr1->SetMarkerStyle(21);
    gr1->SetMarkerSize(0.7);
    gr1->Draw("AP");
}

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
    double eval(double x) {return c*x*x;};
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


    cout << cana::simpson(quad, 0., 10, 1000) << endl;
    MyQuad myq(1.0);
    cout << cana::simpson(&MyQuad::eval, &myq, 0, 10, 1000) << endl;
    auto lmd_quad = [](double x) {return x*x;};
    cout << cana::simpson(lmd_quad, 0, 10, 1000) << endl;

    auto nodes = cana::calc_legendre_nodes(96);
    cout << cana::gauss_quad(nodes, quad, 0., 10) << endl;
    cout << cana::gauss_quad(nodes, &MyQuad::eval, &myq, 0., 10) << endl;
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

void elas_test(double energy = 1147, double angle = 9.03, double nu_beg = 5, double nu_end = 700)
{
    TGraph *g1 = new TGraph();
    TGraph *g2 = new TGraph();
    TGraph *g3 = new TGraph();
    TGraph *g4 = new TGraph();

    CHe3Elas he3_MT("configs/elas_tail.conf");
    CHe3Elas he3_XY("configs/elas_tail.conf");
    CHe3Elas he3_BS("configs/elas_tail.conf");
    he3_MT.SetConfigValue("Radiative Approach", ConfigValue("MT Exact"));
    he3_MT.Configure();
    he3_XY.SetConfigValue("Radiative Approach", ConfigValue("MT Approx"));
    he3_XY.Configure();
    he3_BS.SetConfigValue("Radiative Approach", ConfigValue("BS"));
    he3_BS.Configure();

    TGraph *g1v2 = new TGraph();
    TGraph *g1v3 = new TGraph();

    double rlb = 2.073e-3;
    double rla = 4.492e-2;

    string tail_file = "energy_loss_BS.dat";
    fill_graph(g4, tail_file, 0, 1, 1.0, 0.0);
    ofstream outf("elas_test.txt");

    for(double nu = nu_beg; nu <= nu_end; nu += 1.0)
    {
        double xs_MT = 1000.*he3_MT.GetRadXS(energy, energy - nu, angle*cana::deg2rad, rlb, rla);
        double xs_XY = 1000.*he3_XY.GetRadXS(energy, energy - nu, angle*cana::deg2rad, rlb, rla);
        double xs_BS = 1000.*he3_BS.GetRadXS(energy, energy - nu, angle*cana::deg2rad, rlb, rla);
        g1->SetPoint(g1->GetN(), nu, xs_MT);
        g2->SetPoint(g2->GetN(), nu, xs_XY);
        g1v2->SetPoint(g1v2->GetN(), nu, (xs_XY/xs_MT - 1.)*100.);
        g3->SetPoint(g3->GetN(), nu, xs_BS);
        g1v3->SetPoint(g1v2->GetN(), nu, (xs_BS/xs_MT - 1.)*100.);
        outf << setw(8) << energy
             << setw(8) << nu
             << setw(15) << xs_MT
             << setw(15) << xs_XY
             << setw(15) << xs_BS
             << endl;
    }

    // blue
    g2->SetLineColor(4);
    g1v2->SetLineColor(4);
    // red
    g3->SetLineColor(2);
    g1v3->SetLineColor(2);

    TCanvas *c1 = new TCanvas("unpol tail","unpol tail",200,10,700,1000);
    c1->SetGrid();
    c1->Divide(1, 2);
    c1->cd(1);
    g1->Draw("AC");
    g2->Draw("C");
    g3->Draw("C");
    c1->cd(2);
    c1->DrawFrame(nu_beg, -10., nu_end, 10.);
    g1v3->Draw("AC");
    g1v2->Draw("C");
    /*
    g4->SetMarkerStyle(20);
    g4->Draw("P");
    */
}

void landau_test()
{
    TGraph *g1 = new TGraph();
    TGraph *g2 = new TGraph();
    for(double x = -5; x < 15; x += 0.1)
    {
        g1->SetPoint(g1->GetN(), x, cana::landau_straggle(x, 1, 0, false));
        g2->SetPoint(g2->GetN(), x, ROOT::Math::landau_pdf(x));
    }
    TCanvas *c1 = new TCanvas("Landau dist", "Landau dist", 200, 10, 1200, 500);
    g1->SetMarkerStyle(21);
    g1->Draw("AP");

    g2->SetMarkerStyle(22);
    g2->SetMarkerColor(2);
    g2->Draw("P");
}

#define N_LAMDA 1000
#define N_XS 1000
struct XSDist {double Ep, cdf;};
struct LamdaDist {double lamda, cdf; std::vector<XSDist> xs_dist;};
void calc_xs_dist(CHe3Elas &model, XSDist *dist, int Npoints, double Es, double angle)
{
    //bool MT_method = false;
    //rtails_init(false, false, 180);
    //rtails_set_radl(0., 0., 0., 0.);
    double theta = angle*cana::deg2rad;
    double elas_Ep = Es/(1. + 2.*Es*std::pow(sin(theta/2.), 2)/(3.0149322473*cana::amu));
    double Epmax = elas_Ep - 0.1;
    double Epmin = 0.20*elas_Ep;
    double step = (Epmax - Epmin)/(double)(Npoints - 1);
    dist[0].Ep = Epmin;
    dist[0].cdf = 0.;
    double prev_xs = model.GetRadXS(Es, Epmin, theta, 0., 0.);
    //double prev_xs = rtails_rad_cxsn(Es, Epmin, theta, MT_method);

    for(int i = 1; i < Npoints; ++i)
    {
        double Ep = Epmin + i*step;
        //double xs = rtails_rad_cxsn(Es, Ep, theta, MT_method);
        double xs = model.GetRadXS(Es, Ep, theta, 0., 0.);
        dist[i].Ep = Ep;
        dist[i].cdf = step*(xs + prev_xs)/2. + dist[i-1].cdf;
        prev_xs = xs;
    }
}

double loss_ext_brem(double Ei, double z, double t, double rnd)
{
    double b = 4./3.*(1. + (z + 1.)/9./(log(1194.*std::pow(z, -2./3.)) + z*log(184.15*std::pow(z, -1./3.))));
    return Ei*std::pow(rnd*0.999, 1./(b*t));
}

void energy_loss(double Es = 1147, double angle = 9.03)
{
    int nevents = 500000000;
    double xi1 = 7.417e-3;
    double xi2 = 0.1069;
    double avg_loss1 = 0.23;
    double avg_loss2 = 3.086;
    CHe3Elas model("configs/elas_tail.conf");
    model.SetConfigValue("Radiative Approach", ConfigValue("BS"));
//    model.SetConfigValue("Radiative Approach", ConfigValue("MT Exact"));
//    model.SetConfigValue("Radiative Approach", ConfigValue("MT Approx"));
    model.Configure();

    double m = cana::ele_mass;
    double gamma = Es/m;
    double beta = sqrt(1. - 1./gamma/gamma);
    double Emax = m*beta*beta*gamma*gamma/(1. + gamma);
    double lamda_avg1 = cana::euler - 1. - beta*beta - log(xi1/Emax);
    double lamda_avg2 = cana::euler - 1. - beta*beta - log(xi2/Emax);
    double lamda_max = 0.60715 + 1.1934*lamda_avg1 + (0.67794 + 0.052382*lamda_avg1)*exp(0.94753 + 0.74442*lamda_avg1);
    double lamda_min = lamda_avg1 - avg_loss1/xi1;

    double step = (lamda_max - lamda_min)/(double)N_LAMDA;
    std::vector<LamdaDist> lamda_dist(N_LAMDA + 1);
    LamdaDist fp;
    fp.lamda = lamda_min;
    fp.cdf = 0.;
    fp.xs_dist.resize(N_XS + 1);
    calc_xs_dist(model, &fp.xs_dist.at(0), N_XS + 1, Es, angle);
    lamda_dist.at(0) = fp;
    double luminosity = (double)nevents/fp.xs_dist.back().cdf;

    double plandau = cana::landau_fit(fp.lamda);

    for(int i = 1; i <= N_LAMDA; ++i)
    {
        double lamda = lamda_min + i*step;
        double landau = cana::landau_fit(lamda);
        LamdaDist &curr = lamda_dist.at(i);
        LamdaDist &prev = lamda_dist.at(i - 1);
        curr.lamda = lamda;
        curr.xs_dist.resize(N_XS + 1);
        double Es1 = Es - (lamda - lamda_avg1)*xi1 - avg_loss1;
        calc_xs_dist(model, &curr.xs_dist.at(0), N_XS + 1, Es1, angle);
        curr.cdf = step*(landau + plandau)/2. + prev.cdf;
        plandau = landau;
    }

    random_device rd;
    mt19937 rng{rd()};
    uniform_real_distribution<double> uni_dist(0., 1.);
    auto comp = [](LamdaDist p, double val) {return p.cdf - val;};
    auto comp2 = [](XSDist p, double val) {return p.cdf - val;};

    int counts[700];
    for(int i = 0; i < 700; ++i)
    {
        counts[i] = 0;
    }

    for(int i = 0; i < nevents; ++i)
    {
        double rnd = uni_dist(rng)*lamda_dist.back().cdf;
        auto itv = cana::binary_search_interval(lamda_dist.begin(), lamda_dist.end(), rnd, comp);
        double lamda = cana::linear_interp(itv.first->cdf, itv.first->lamda,
                                           itv.second->cdf, itv.second->lamda,
                                           rnd);
        double e_loss1 = (lamda - lamda_avg1)*xi1 + avg_loss1;
        double rn = uni_dist(rng);
        double rnd1 = rn*itv.first->xs_dist.back().cdf;
        double rnd2 = rn*itv.second->xs_dist.back().cdf;
        auto itv1 = cana::binary_search_interval(itv.first->xs_dist.begin(),
                                                 itv.first->xs_dist.end(),
                                                 rnd1,
                                                 comp2);
        auto itv2 = cana::binary_search_interval(itv.second->xs_dist.begin(),
                                                 itv.second->xs_dist.end(),
                                                 rnd2,
                                                 comp2);
        double Ep1 = cana::linear_interp(itv1.first->cdf, itv1.first->Ep,
                                         itv1.second->cdf, itv1.second->Ep,
                                         rnd1);
        double Ep2 = cana::linear_interp(itv2.first->cdf, itv2.first->Ep,
                                         itv2.second->cdf, itv2.second->Ep,
                                         rnd2);
        double Ep = cana::linear_interp(itv.first->cdf, Ep1,
                                        itv.second->cdf, Ep2,
                                        rnd);

        double rnd3 = uni_dist(rng)*lamda_dist.back().cdf;
        auto itvp = cana::binary_search_interval(lamda_dist.begin(), lamda_dist.end(), rnd3, comp);
        double lamdap = cana::linear_interp(itvp.first->cdf, itvp.first->lamda,
                                            itvp.second->cdf, itvp.second->lamda,
                                            rnd3);
        double ioni_loss2 = (lamdap - lamda_avg2)*xi2 + avg_loss2;
        double brem_loss2 = loss_ext_brem(Ep, 11.257, 4.492e-02, uni_dist(rng));

        int Epp = int(Es - (Ep - ioni_loss2 - brem_loss2) + 0.5);
        if(Epp < 700 && Epp >= 0)
        {
            counts[Epp]++;
        }
    }

    std::ofstream outf("energy_loss.dat");
    for(int i = 0; i < 700; ++i)
    {
        outf << std::setw(8) << i
             << std::setw(20) << ((double)counts[i]/luminosity)*1e3
             << std::endl;
    }
    outf.close();
}

void show_energy_loss()
{
    ConfigParser c_parser;
    c_parser.ReadFile("energy_loss_MT.dat");
    TGraph *gr1 = new TGraph();
    TGraph *gr2 = new TGraph();
    TGraph *gr3 = new TGraph();
    TGraph *gr2b = new TGraph();
    TGraph *gr3b = new TGraph();

    double nu, val;
    while(c_parser.ParseLine())
    {
        c_parser >> nu >> val;
        gr1->SetPoint(gr1->GetN(), nu, val);
    }

    c_parser.ReadFile("energy_loss_BS.dat");

    while(c_parser.ParseLine())
    {
        c_parser >> nu >> val;
        gr2->SetPoint(gr2->GetN(), nu, val);
    }

    c_parser.ReadFile("energy_loss_XY.dat");

    while(c_parser.ParseLine())
    {
        c_parser >> nu >> val;
        gr3->SetPoint(gr3->GetN(), nu, val);
    }

    for(int i = 0; i < gr1->GetN(); ++i)
    {
        gr2b->SetPoint(i, gr2->GetX()[i], (gr2->GetY()[i] - gr1->GetY()[i])/gr1->GetY()[i]*100.);
        gr3b->SetPoint(i, gr3->GetX()[i], (gr3->GetY()[i] - gr1->GetY()[i])/gr1->GetY()[i]*100.);
    }

    TCanvas *c1 = new TCanvas("unpol tail","unpol tail",200,10,700,1000);
    c1->SetGrid();
    c1->Divide(1, 2);

    gr1->SetMarkerStyle(20);
    gr2->SetMarkerStyle(21);
    gr2->SetMarkerColor(2);
    gr2b->SetMarkerStyle(21);
    gr2b->SetMarkerColor(2);
    gr3->SetMarkerStyle(22);
    gr3->SetMarkerColor(4);
    gr3b->SetMarkerStyle(22);
    gr3b->SetMarkerColor(4);


    c1->cd(1);
    gr1->Draw("AP");
    gr2->Draw("P");
    gr3->Draw("P");

    c1->cd(2);
    gr2b->Draw("AP");
    gr3b->Draw("P");
    gr2b->GetYaxis()->SetRangeUser(-10, 10);
}

