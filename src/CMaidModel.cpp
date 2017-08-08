#include "CMaidModel.h"
#include "ConfigParser.h"
#include "canalib.h"
#include <utility>



static double M_He3 = 3.0149322473*cana::amu;

CMaidModel::CMaidModel(const char *data)
{
    LoadData(data);
}

CMaidModel::~CMaidModel()
{
    // place holder
}

void CMaidModel::LoadData(const char *path)
{
    data.clear();

    ConfigParser c_parser;
    if(!c_parser.OpenFile(path)) {
        std::cerr << "Cannot open file " << path << std::endl;
        return;
    }

    double q2, w, g1p, g1n, g2p, g2n;
    while(c_parser.ParseLine())
    {
        c_parser >> q2 >> w >> g1p >> g1n >> g2p >> g2n;

        q2 *= 1e6; // convert from GeV^2 to MeV^2

        if(data.size() && q2 == data.back().q2)
        {
            data.back().wset.emplace_back(w, g1p, g2p, g1n, g2n);
        } else {
            auto itp = cana::binary_search_interval(data.begin(), data.end(), q2);
            // cannot find the q2
            if(itp.second == data.end() || itp.second != itp.first) {
                DataSet new_set(q2);
                new_set.wset.emplace_back(w, g1p, g2p, g1n, g2n);
                data.insert(itp.second, new_set);
            } else {
                itp.second->wset.emplace_back(w, g1p, g2p, g2n, g2n);
            }
        }
    }

    for(auto &dset : data)
    {
        auto lmda = [](const DataPoint &p1, const DataPoint &p2) {return p1.w < p2.w;};
        std::sort(dset.wset.begin(), dset.wset.end(), lmda);
    }
}

inline CMaidModel::DataPoint interp_point(const CMaidModel::DataPoint &p1,
                                          const CMaidModel::DataPoint &p2,
                                          double a, double b, double c)
{
    if(c < 1e-20) return p1;

    return CMaidModel::DataPoint(0.,
                                 (a*p2.p.g1 + b*p1.p.g1)/c,
                                 (a*p2.p.g2 + b*p1.p.g2)/c,
                                 (a*p2.p.g1 + b*p1.p.g1)/c,
                                 (a*p2.n.g2 + b*p1.n.g2)/c);
}

bool interp_wset(const CMaidModel::DataSet &dset, double W, CMaidModel::DataPoint &res)
{
    const auto &wset = dset.wset;
    auto it2 = cana::binary_search_interval(wset.begin(), wset.end(), W);

    if(it2.first == wset.end()) {
        std::cerr << "W = " << W << " exceeds the range min = " << wset.front().w
                  << ", max = " << wset.back().w
                  << std::endl;
        return false;
    }

    double a = W - it2.first->w;
    double b = it2.second->w - W;
    double c = it2.second->w - it2.first->w;
    res = interp_point(*it2.first, *it2.second, a, b, c);
    return true;
}


bool CMaidModel::Interp(double Q2, double W, MaidValue &p, MaidValue &n)
const
{
    if(data.empty()) {
        std::cerr << "No data read-in." << std::endl;
        return false;
    }

    auto it = cana::binary_search_interval(data.begin(), data.end(), Q2);
    if(it.second == data.end()) {
        std::cerr << "Q2 = " << Q2 << " exceeds the range min = " << data.front().q2
                  << ", max = " << data.back().q2
                  << std::endl;
        return false;
    }

    if(it.first == it.second) {
        DataPoint point;
        if(!interp_wset(*it.first, W, point))
            return false;
        p = point.p;
        n = point.n;
    } else {
        DataPoint point, point1, point2;
        if(!interp_wset(*it.first, W, point1) || !interp_wset(*it.second, W, point2))
            return false;

        double a = Q2 - it.first->q2;
        double b = it.second->q2 - Q2;
        double c = it.second->q2 - it.first->q2;
        point = interp_point(point1, point2, a, b, c);
        p = point.p;
        n = point.n;
    }

    return true;
}

inline double momentum_square(const double *p)
{
    return p[1]*p[1] + p[2]*p[2] + p[3]*p[3];
}

// p is four vector p[0] = energy, (p[1], p[2], p[3]) = p->
void calc_DN(const double *p, double gamma, double Fs, double Ft, double DN[][2])
{
    double gamma2 = gamma*gamma, gamma4 = gamma2*gamma2;
    double p2 = momentum_square(p);
    double npz = p[3]/sqrt(p2), npz2 = npz*npz, npz4 = npz2*npz2;
    double c1 = p[3]/gamma/M_He3;
    double c2 = p2/M_He3/M_He3;

    DN[0][0] = Fs + (3. - gamma2)/6./gamma2*(3.*npz2 - 1.)*Ft + c1*(Fs + 2./3.*Ft)
               + c2*((3. - gamma2)*npz2 - 1. - gamma2)/12./gamma2*(3.*Fs - Ft);

    DN[0][1] = (-(3.*npz2 - 1.)/2./gamma2*Ft + c1*(Fs + (3./2.*npz2 - 5./6.)*Ft)
                - c2*((1. + npz2*(4.*gamma2 - 3.))/4./gamma2*Fs
                      + (5. + 18.*npz4*gamma2 - 5.*npz2*(3. + 2.*gamma2))/12./gamma2*Ft))*(gamma2 - 1.);

    DN[1][0] = (3.*npz2 - 1.)/2./gamma2*Ft - c1*(Fs + 2./3.*Ft) - c2*(3.*npz2 - 1.)/12./gamma2*(3.*Fs - Ft);
    DN[1][1] = Fs + (2.*gamma2 - 3.)/6./gamma2*(3.*npz2 - 1.)*Ft
               + c1*((1. - gamma2)*Fs + (-5./6. + gamma2/3. + npz2*(3./2. - gamma2))*Ft)
               + c2*((npz2*(3. - 6.*gamma2 + 4.*gamma4)- 1. - 2.*gamma2)/4./gamma2*Fs
                     + (5. - 2.*gamma2*(1. + 3.*npz2) + 4.*npz2*gamma4)/12./gamma2*(3.*npz2 - 1.)*Ft);
}

static double twopi4 = std::pow(2.*cana::pi, 4);
static auto nodes = cana::calc_legendre_nodes(96);

double f00dpxdpydpz(double pz, double py, double px, double y, double gamma)
{
    double eps = (y - 1)*M_He3 - gamma*pz;
    double p[4] = {eps + M_He3, px, py, pz};
    double D[2][2];
    calc_DN(p, gamma, 0.856, 0.013, D);
    return D[0][0]/twopi4;
}

double f00dpxdpy(double py, double px, double y, double gamma)
{
    return cana::gauss_quad(nodes, f00dpxdpydpz, 0., 10., py, px, y, gamma);
}

double f00dpx(double px, double y, double gamma)
{
    return cana::gauss_quad(nodes, f00dpxdpy, 0., 10., px, y, gamma);
}

double f00(double y, double gamma)
{
    return cana::gauss_quad(nodes, f00dpx, 0., 10., y, gamma);
}

inline double P1N_EPA(double Fs0, double Fs2, double Ft2)
{
    return Fs0 - 1./3.*(Fs2 - Ft2/3.);
}

inline double P2N_EPA(double Fs0, double Fs2, double Ft2)
{
    return Fs0 - 2./3.*(Fs2 - Ft2/15.);
}

bool CMaidModel::GetHe3gg_EPA(double Q2, double W, double &g1, double &g2)
const
{
    MaidValue p, n;
    if(!Interp(Q2, W, p, n)) return false;

    // KPSV
    double P1p = -0.028, P1n = 0.851, P2p = -0.028, P2n = 0.844;
    // SS
    //double P1p = -0.021, P1n = 0.884, P2p = -0.021, P2n = 0.878;

    g1 = p.g1*2.*P1p + n.g1*P1n;
    g2 = p.g2*2.*P2p + n.g2*P2n;

    return true;
}
