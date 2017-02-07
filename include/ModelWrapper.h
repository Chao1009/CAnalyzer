// model wrapper in C format
extern "C"
{
    // units
    // MeV, MeV^2, radian, ub/MeV/sr
    void Bosted_f1f2in09(double Z, double A, double Q2, double W2,
                         double *F1, double *F2, double *rc);

    void Bosted_f1f2qe09(double Z, double A, double Q2, double W2,
                         double *F1, double *F2);

    void Bosted_xs(double Z, double A, double Ei, double Ef, double theta,
                   double *xs);

    void QFS_xs(double Z, double A, double Ei, double Ef, double theta,
                double *xs);
}


