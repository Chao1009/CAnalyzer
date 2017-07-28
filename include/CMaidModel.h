#ifndef C_MAID_MODEL_H
#define C_MAID_MODEL_H

#include <vector>

class CMaidModel
{
public:
    struct MaidValue
    {
        double g1, g2;

        MaidValue() : g1(0.), g2(0.) {}
        MaidValue(double g1i, double g2i)
        : g1(g1i), g2(g2i) {}
    };

    struct DataPoint
    {
        double w;
        MaidValue p, n;

        DataPoint() : w(0.) {}
        DataPoint(double wi,
                  double g1p, double g2p,
                  double g1n, double g2n)
        : w(wi), p(g1p, g2p), n(g1n, g2n) {}

        bool operator == (double wi) const {return w == wi;}
        bool operator != (double wi) const {return w != wi;}
        bool operator > (double wi) const {return w > wi;}
        bool operator < (double wi) const {return w < wi;}
    };

    struct DataSet
    {
        double q2;
        std::vector<DataPoint> wset;

        DataSet() : q2(0.) {}
        DataSet(double q2i) : q2(q2i) {}

        bool operator == (double q2i) const {return q2 == q2i;}
        bool operator != (double q2i) const {return q2 != q2i;}
        bool operator > (double q2i) const {return q2 > q2i;}
        bool operator < (double q2i) const {return q2 < q2i;}
    };

public:
    CMaidModel(const char *data_path);
    virtual ~CMaidModel();

    void LoadData(const char *path);

    bool Interp(double Q2, double W, MaidValue &p, MaidValue &n) const;

private:
    std::vector<DataSet> data;
};

#endif //C_MAID_MODEL_H
