#include "CModelWrapper.h"
#include "canalib.h"
#include <iostream>

#define TARZ 2
#define TARA 3.0149322473

// define some operators for simplified coding
CModelWrapper::MPoint operator *(const CModelWrapper::MPoint &p, const double &val)
{
    CModelWrapper::MPoint newp;
    newp.F1_in = p.F1_in*val;
    newp.F2_in = p.F2_in*val;
    newp.F1_qe = p.F1_qe*val;
    newp.F2_qe = p.F2_qe*val;
    return newp;
}

CModelWrapper::MPoint operator *(const double &val, const CModelWrapper::MPoint &p)
{
    return p*val;
}

CModelWrapper::MPoint operator /(const CModelWrapper::MPoint &p, const double &val)
{
    CModelWrapper::MPoint newp;
    newp.F1_in = p.F1_in/val;
    newp.F2_in = p.F2_in/val;
    newp.F1_qe = p.F1_qe/val;
    newp.F2_qe = p.F2_qe/val;
    return newp;
}

CModelWrapper::MPoint operator +(const CModelWrapper::MPoint &p1, const CModelWrapper::MPoint &p2)
{
    CModelWrapper::MPoint newp;
    newp.F1_in = p1.F1_in + p2.F1_in;
    newp.F2_in = p1.F2_in + p2.F2_in;
    newp.F1_qe = p1.F1_qe + p2.F1_qe;
    newp.F2_qe = p1.F2_qe + p2.F2_qe;
    return newp;
}

CModelWrapper::MPoint operator -(const CModelWrapper::MPoint &p1, const CModelWrapper::MPoint &p2)
{
    CModelWrapper::MPoint newp;
    newp.F1_in = p1.F1_in - p2.F1_in;
    newp.F2_in = p1.F2_in - p2.F2_in;
    newp.F1_qe = p1.F1_qe - p2.F1_qe;
    newp.F2_qe = p1.F2_qe - p2.F2_qe;
    return newp;
}


CModelWrapper::CModelWrapper()
: norm(1.0)
{
    // place holder
}

CModelWrapper::~CModelWrapper()
{
    // place holder
}

void CModelWrapper::SetRange(double Qsq_min, double Qsq_max, int Qsq_bins,
                             double W_min, double W_max, int W_bins)
{
    // invalid
    if(Qsq_bins <= 0 || W_bins <= 0 || Qsq_max <= Qsq_min || W_max <= W_min) {
        std::cerr << "Failed to initialize CModelWrapper, invalid inputs." << std::endl;
        return;
    }

    norm = 1.0;

    double Qsq_step = (Qsq_max - Qsq_min)/(double)Qsq_bins;
    double W_step = (W_max - W_min)/(double)W_bins;

    model_sets.clear();
    model_sets.reserve(Qsq_bins);

    for(int i = 0; i <= Qsq_bins; ++i)
    {
        MSet new_set;
        new_set.Q2 = (Qsq_min + i*Qsq_step)*1e-6; // convert to GeV^2
        new_set.points.reserve(W_bins);
        for(int j = 0; j <= W_bins; ++j)
        {
            MPoint new_point;
            new_point.W = (W_min + j*W_step)*1e-3; // convert to GeV
            double W2 = new_point.W*new_point.W;
            double rc; // not used
            Bosted_f1f2in09(TARZ, TARA, new_set.Q2, W2,
                            &new_point.F1_in, &new_point.F2_in, &rc);
            Bosted_f1f2qe09(TARZ, TARA, new_set.Q2, W2,
                            &new_point.F1_qe, &new_point.F2_qe);
            new_set.points.push_back(new_point);
        }
        model_sets.push_back(new_set);
    }
}

double CModelWrapper::GetCrossSection(const double &Ei, const double &Ef, const double &theta)
const
{
    double sin2 = std::pow(sin(theta/2.), 2);
    double nu = (Ei - Ef)*1e-3;
    double Q2 = 4.0*Ei*Ef*sin2*1e-6;
    double W = sqrt(cana::proton_mass*1e-3*(cana::proton_mass*1e-3 + 2.*nu) - Q2);

    auto it_pair = cana::binary_search_interval(model_sets.begin(), model_sets.end(), Q2);
    MPoint interp;

    if(it_pair.second == model_sets.end() || it_pair.first == model_sets.end()) {
        std::cerr << "Cannot interpolate for Q2 = " << Q2 << std::endl;
        return 0.;
    }

    if(it_pair.first == it_pair.second) {
        interp = it_pair.first->Interp(W);
    } else {
        interp = (it_pair.first->Interp(W)*(it_pair.second->Q2 - Q2)
                  + it_pair.second->Interp(W)*(Q2 - it_pair.first->Q2))
                 / (it_pair.second->Q2 - it_pair.first->Q2);
    }

    double Mott = std::pow(2./Q2*cana::alpha*Ef*cos(theta/2.)*cana::hbarc, 2)*1e-6;
    double xs = 2./(cana::proton_mass*1e-3) * (interp.F1_qe + interp.F1_in)
                * std::pow(tan(theta/2.), 2)
                + (interp.F2_qe + interp.F2_in)/nu;

    return norm*xs*Mott/1.0e2;
}

CModelWrapper::MPoint CModelWrapper::MSet::Interp(const double &W)
const
{
    // search the position of w
    auto it_pair = cana::binary_search_interval(points.begin(), points.end(), W);

    // not found
    if(it_pair.second == points.end() || it_pair.first == points.end()) {
        std::cerr << "Cannot interpolate for Q2 = " << Q2 << ", W = " << W << std::endl;
        return MPoint();
    }

    // exact matched
    if(it_pair.first == it_pair.second)
        return *it_pair.first;

    const MPoint &p1 = *it_pair.first;
    const MPoint &p2 = *it_pair.second;

    return (p1*(p2.W - W) + p2*(W - p1.W))/(p2.W - p1.W);
}
