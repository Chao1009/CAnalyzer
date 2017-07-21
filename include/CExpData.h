#ifndef C_EXP_DATA_H
#define C_EXP_DATA_H

#include <vector>
#include <string>

class CExpData
{
public:
    struct Settings
    {
        double angle, targetZ, targetA;
        double radl_factor_before, radl_factor_after;
        std::string data_dir, accpt_dir, coll_dir;
    };

    struct DataPoint
    {
        // read from data file
        double nu;       // nu (MeV)
        double stat;     // statistical error
        double syst;     // systematic error

        // calculated
        double Ep;       // final energy (MeV)
        double W;        // invariant mass (MeV)
        double Q2;       // four momentum transfer square (MeV^2)
        double rad;      // Radiated cross section (nb/MeV/sr)
        double born;     // Born cross section (nb/MeV/sr)
        double last;     // Save info from last iteration

        // constructors
        DataPoint() {};
        DataPoint(double n, double c, double st, double sy)
        : nu(n), stat(st), syst(sy), rad(c), born(c), last(c)
        {};

        // for binary search
        bool operator <(const double &val) const {return W < val;}
        bool operator >(const double &val) const {return W > val;}
        bool operator ==(const double &val) const {return W == val;}
        bool operator !=(const double &val) const {return W != val;}
    };

    struct DataSet
    {
        // read from data file
        double energy;        // incident electron energy
        double angle;         // scattering angle (degree)
        double radl_before;   // radiation length before target
        double radl_after;    // radiation length after target
        double coll_before;   // collision thickness before target
        double coll_after;    // collision thickness after target
        double radl_wall;     // radiation length of the target cell wall
        double error;         // relative error for RC
        double normalization; // normalization factor
        bool non_rad;         // non radiated means born cross section from file
        std::vector<DataPoint> data;

        // information related
        std::string data_file, data_label, accpt_file, coll_file;

        // constructors
        DataSet()
        : energy(0.), angle(0.),
          radl_before(0.), radl_after(0.), coll_before(0.), coll_after(0.),
          radl_wall(0.), error(0.), normalization(1.), non_rad(false)
        {};

        void ReadConfig(const std::string &conf_str, const Settings &gset);
        void ReadData(const std::string &path, const std::string &label);
        double Interp(const double &v) const;
        void CopySettings(const DataSet &ref);

        bool operator <(const double &val) const {return energy < val;}
        bool operator >(const double &val) const {return energy > val;}
        bool operator ==(const double &val) const {return energy == val;}
        bool operator !=(const double &val) const {return energy != val;}
    };


    CExpData(const std::string &config = "");
    virtual ~CExpData();

    void ReadConfigFile(const std::string &path, bool verbose = true);
    void ReadSettings(const std::string &conf_str, const std::string &path);
    void ReadDataPoints(DataSet &dset, const std::string &path, const std::string &label);

    double GetCrossSection(const double &E0, const double &Eb) const;
    void DataUpdate();
    void SaveResult(const std::string &path, bool save_model = true) const;

    inline const Settings &GetSettings() const {return settings;}
    inline double Angle() const {return settings.angle;}
    inline double TargetZ() const {return settings.targetZ;}
    inline double TargetA() const {return settings.targetA;}
    inline size_t Size() const {return data_sets.size();}
    inline bool Empty() const {return data_sets.empty();}
    inline const std::vector<DataSet> &GetSets() const {return data_sets;}
    inline const DataSet &GetSet(size_t i) const {return data_sets.at(i);}
    inline std::vector<DataSet> &GetSets() {return data_sets;}
    inline DataSet &GetSet(size_t i) {return data_sets.at(i);}
    inline bool InRange(const double &E0)
    const
    {
        if(data_sets.empty())
            return false;
        return (E0 >= data_sets.front().energy && E0 <= data_sets.back().energy);
    }

private:
    Settings settings;
    double targetM, sin2;
    std::vector<DataSet> data_sets;
};

#endif
