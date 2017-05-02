#ifndef C_MODEL_WRAPPER_H
#define C_MODEL_WRAPPER_H

// model wrappers
extern "C"
{
    // units
    // GeV, GeV^2
    void Bosted_f1f2in09(double Z, double A, double Q2, double W2,
                         double *F1, double *F2, double *rc);

    void Bosted_f1f2qe09(double Z, double A, double Q2, double W2,
                         double *F1, double *F2);

    // units
    // MeV, MeV^2, radian, ub/MeV/sr
    void Bosted_xs(double Z, double A, double Ei, double Ef, double theta,
                   double *xs);

    void QFS_xs(double Z, double A, double Ei, double Ef, double theta,
                double *xs);
}

#include <vector>

class CModelWrapper
{
public:
    struct MPoint
    {
        double W, F1_in, F2_in, F1_qe, F2_qe;

        MPoint() : W(0.), F1_in(0.), F2_in(0.), F1_qe(0.), F2_qe(0.) {}

        // for binary search
        bool operator <(const double &val) const {return W < val;}
        bool operator >(const double &val) const {return W > val;}
        bool operator ==(const double &val) const {return W == val;}
        bool operator !=(const double &val) const {return W != val;}
    };

    struct MSet
    {
        double Q2;
        std::vector<MPoint> points;

        MPoint Interp(const double &W) const;

        // for binary search
        bool operator <(const double &val) const {return Q2 < val;}
        bool operator >(const double &val) const {return Q2 > val;}
        bool operator ==(const double &val) const {return Q2 == val;}
        bool operator !=(const double &val) const {return Q2 != val;}
    };

public:
    CModelWrapper();
    ~CModelWrapper();

    // in MeV^2 and MeV
    void SetRange(double Qsq_min = 10, double Qsq_max = 2e5, int Q2_bins = 500,
                  double W_min = 900, double W_max = 1900, int W_bins = 1000);
    // in MeV and rad, out nb/MeV/sr
    double GetCrossSection(const double &Ei, const double &Ef, const double &theta) const;
    double GetNormalization() const {return norm;};
    void SetNormalization(const double &n) {norm = n;};

private:
    double Q2_bins, W_bins, norm;
    std::vector<MSet> model_sets;
};

#endif
