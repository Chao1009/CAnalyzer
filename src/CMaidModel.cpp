#include "CMaidModel.h"
#include "ConfigParser.h"
#include "canalib.h"
#include <utility>



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
