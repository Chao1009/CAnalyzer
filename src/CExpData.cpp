#include "CExpData.h"
#include "ConfigParser.h"
#include "ConfigObject.h"
#include "canalib.h"
#include <fstream>
#include <iostream>
#include <iomanip>



// constructor
CExpData::CExpData()
{
    // place holder
}

// destructor
CExpData::~CExpData()
{
    // place holder
}

// helper function to determine the file path
std::string combine_path(const std::string &dir, const std::string &path)
{
    // invalid input or absolute path, do not combine
    if(dir.empty() || (path.size() && path.front() == '/'))
        return path;

    if(dir.back() != '/')
        return dir + "/" + path;
    return dir + path;
}

void CExpData::ReadConfigFile(const std::string &path, bool verbose)
{
    data_sets.clear();

    // read in file
    std::string buffer = ConfigParser::file_to_string(path);

    // remove comments
    ConfigParser::comment_between(buffer, "/*", "*/");
    ConfigParser::comment_line(buffer, "//", "\n");
    ConfigParser::comment_line(buffer, "#", "\n");

    // break into blocks
    auto blocks = ConfigParser::break_into_blocks(buffer, "{", "}");

    // read configuration blocks
    for(auto &block : blocks)
    {
        if(ConfigParser::case_ins_equal(block.label, "SETTING")) {

            ReadSettings(block.content, path);

        } else if(ConfigParser::case_ins_equal(block.label, "DATASET")) {

            DataSet new_set;
            new_set.ReadConfig(block.content);

            if((new_set.energy > 0.) &&
               !cana::is_in(new_set.energy, data_sets.begin(), data_sets.end())) {
                data_sets.emplace_back(new_set);
            } else {
                std::cerr << "Error: Invalid energy information or existing a data"
                          << "set with the same energy = " << new_set.energy
                          << std::endl;
            }
        }
    }

    // read-in data points according to the set data file
    for(auto &dset : data_sets) {
        dset.ReadData(combine_path(settings.data_dir, dset.data_file), dset.data_label);
    }

    // initialize all data sets and sort them in energy transcendent order
    DataUpdate();

    if(verbose) {
        // show how many data sets are read-in
        std::cout << "Read " << data_sets.size() << " data sets." << std::endl;

        for(auto &dset : data_sets)
        {
            std::cout << "Energy: " << std::setw(8) << dset.energy
                      << ", DataPoints:" << std::setw(8) << dset.data.size()
                      << ", Type: " << ((dset.non_rad) ? "Born" : "Rad")
                      << std::endl;
        }
    }
}

// save result to the path
void CExpData::SaveResult(const std::string &path)
const
{
    std::ofstream output(path);
    if(!output.is_open()) {
        std::cerr << "Cannot open output file "
                  << "\"" << path << "\""
                  << std::endl;
        return;
    }

    output << "#"
           << std::setw(7) << "ENERGY"
           << std::setw(8) << "NU"
           << std::setw(15) << "SIGRAD"
           << std::setw(15) << "SIGBORN"
           << std::setw(15) << "STAT. ERR."
           << std::setw(15) << "SYST. ERR."
//           << std::setw(15) << "ITER. CHANGE"
           << std::endl;
    for(auto &dset : data_sets)
    {
        if(dset.non_rad)
            continue;

        for(auto &point : dset.data)
        {
            double xsr = point.rad;
            double xsb = point.born;
//            double xsl = point.last;
            double scale = xsb/xsr;
            double stat_err = point.stat*scale;
            double syst_err = point.syst*scale;
            double rc_err = dset.error*std::fabs(xsb - xsr);
            syst_err = sqrt(syst_err*syst_err + rc_err*rc_err);
            output << std::setw(8) << dset.energy
                   << std::setw(8) << point.nu
                   << std::setw(15) << xsr
                   << std::setw(15) << xsb
                   << std::setw(15) << stat_err
                   << std::setw(15) << syst_err
//                   << std::setw(15) << xsb - xsl
                   << std::endl;
        }
    }
}

// helper function
inline void warn_setting_miss(bool warn, std::string term)
{
    if(warn)
        std::cerr << "Error: Missing (" << term << ") in global settings." << std::endl;
}

// read settings
void CExpData::ReadSettings(const std::string &config, const std::string &path)
{
    ConfigObject conf_obj;
    conf_obj.ReadConfigString(config);

    // some thing must be set
    warn_setting_miss(!conf::update_config(conf_obj, "Angle", settings.angle), "Angle");
    warn_setting_miss(!conf::update_config(conf_obj, "Target Z", settings.targetZ), "Target Z");
    warn_setting_miss(!conf::update_config(conf_obj, "Target A", settings.targetA), "Target A");

    std::string conf_dir = ConfigParser::decompose_path(path).dir;
    conf::update_config(conf_obj, "Data Folder", settings.data_dir);
    conf::update_config(conf_obj, "Collimator Folder", settings.coll_dir);
    conf::update_config(conf_obj, "Acceptance Folder", settings.accpt_dir);
    settings.data_dir = combine_path(conf_dir, settings.data_dir);
    settings.coll_dir = combine_path(conf_dir, settings.coll_dir);
    settings.accpt_dir = combine_path(conf_dir, settings.accpt_dir);
}

// update data, sorting data points and sets
void CExpData::DataUpdate()
{
    for(auto &dset : data_sets)
    {
        // sort data points by nu
        std::sort(dset.data.begin(), dset.data.end(),
                  [] (const DataPoint &p1, const DataPoint &p2)
                  {
                      return p1.nu < p2.nu;
                  });
    }

    // sort data sets by energy
    std::sort(data_sets.begin(), data_sets.end(),
             [] (const DataSet &set1, const DataSet &set2)
             {
                 return set1.energy < set2.energy;
             });
}

// interpolates or extrapolates
// scaling for mott cross section under the same scattering angle
double CExpData::GetCrossSection(const double &E0, const double &Eb)
const
{
    double res;
    double w = (E0 - Eb)/E0;

    auto inter = cana::binary_search_interval(data_sets.begin(), data_sets.end(), E0);

    // extrapolation, may result in a huge error
    if(inter.first == data_sets.end() || inter.second == data_sets.end()) {
        if(E0 < data_sets.front().energy) {
            res = data_sets.front().Interp(w)*std::pow(data_sets.front().energy/E0, 2);
        } else {
            res = data_sets.back().Interp(w)*std::pow(data_sets.back().energy/E0, 2);
        }
    // exact match
    } else if(inter.first == inter.second) {
        res = inter.first->Interp(w)*std::pow(inter.first->energy/E0, 2);
    // interpolation between two sets
    } else {
        double E1 = inter.first->energy;
        double E2 = inter.second->energy;
        res = (inter.first->Interp(w)*std::pow(inter.first->energy/E0, 2)*(E2 - E0)
               + inter.second->Interp(w)*std::pow(inter.first->energy/E0, 2)*(E0 - E1))
              / (E2 - E1);
    }

    // scale by Mott
    return res;
}

// read configuration for the data set
void CExpData::DataSet::ReadConfig(const std::string &config)
{
    ConfigObject conf_obj;

    conf_obj.ReadConfigString(config);

    // connect data set variables to configurations
    conf::update_config(conf_obj, "Energy", energy);
    conf::update_config(conf_obj, "Radiation Length Before", radl_before);
    conf::update_config(conf_obj, "Radiation Length After", radl_after);
    conf::update_config(conf_obj, "Collisional Loss Before", coll_before);
    conf::update_config(conf_obj, "Collisional Loss After", coll_after);
    conf::update_config(conf_obj, "Ice Before", ice_before);
    conf::update_config(conf_obj, "Ice After", ice_after);
    conf::update_config(conf_obj, "RC Error", error);
    conf::update_config(conf_obj, "Model", non_rad);
    conf::update_config(conf_obj, "Normalization", normalization);
    conf::update_config(conf_obj, "Data File", data_file);
    conf::update_config(conf_obj, "Data Label", data_label);
    conf::update_config(conf_obj, "Acceptance File", accpt_file);
    conf::update_config(conf_obj, "Collimator File", coll_file);
}

// read experimental data in the format line by line
// comment marks are # and //
// splitter can be tab, space and comma
// additional white spaces (space and tab) will be trimmed
// expected 4 ~ 5 columns (label can be neglected if input data_label is empty)
// *label*, nu, cxsn, stat. error, syst. error
void CExpData::DataSet::ReadData(const std::string &path, const std::string &label)
{
    ConfigParser c_parser;

    if(!c_parser.ReadFile(path)) {
        std::cerr << "Cannot open file \"" << path << "\", no data points read "
                  << "for data set at energy = " << energy
                  << std::endl;
        return;
    }

    // clean up old data
    data.clear();

    while(c_parser.ParseLine())
    {
        // check if label agrees
        if(c_parser.NbofElements() == 5) {
            std::string plabel = c_parser.TakeFirst().String();
            if(plabel != label)
                continue;
        }

        // check format
        if(!c_parser.CheckElements(4))
            continue;

        // new data points
        double nu, cxsn, stat, syst;
        c_parser >> nu >> cxsn >> stat >> syst;
        // apply normalization
        cxsn *= normalization;
        stat *= normalization;

        DataPoint new_point(nu, cxsn, stat, syst);
        // calculate ep
        new_point.Ep = energy - nu;
        new_point.v = nu/energy;

        data.push_back(std::move(new_point));
    }
}

// interpolation between a data set according to (1 - Ep/E0)
// return cross section normalized to mott factor
double CExpData::DataSet::Interp(const double &w)
const
{
    // for low nu region, it is impossible to extrapolate
    if(w < data.front().v)
        return 0.;
    // assuming uniform distribution for high nu region (DIS)
    if(w >= data.back().v)
        return data.back().born;

    // search the position of w
    auto it_pair = cana::binary_search_interval(data.begin(), data.end(), w);

    // not found
    if(it_pair.second == data.end() || it_pair.first == data.end())
        return 0;

    // exact matched
    if(it_pair.first == it_pair.second)
        return it_pair.first->born;

    // only have 2 points, do a straight line interpolation
    if(it_pair.first == data.begin()) {
        const DataPoint &p1 = *it_pair.first;
        const DataPoint &p2 = *it_pair.second;

        double _interp = p1.born*(p2.v - w) + p2.born*(w - p1.v);

        // LATEST CHANGE: removed unknown constant 0.00001 here, 
        // probably some protection for two same points
        // return _interp/(p2.v - p1.v + 0.00001);
        return _interp/(p2.v - p1.v);
    }

    // 3 points parabolic fit
    const DataPoint &p1 = *(it_pair.first - 1);
    const DataPoint &p2 = *it_pair.first;
    const DataPoint &p3 = *it_pair.second;
    double x, xp, xm, a, b, c;
    x = w - p2.v;
    xp = p3.v - p2.v;
    xm = p1.v - p2.v;
    a = p2.born;
    c = (p3.born - p1.born)*(xp + xm)/(xp - xm) - (p3.born + p1.born - 2*p2.born);
    c /= xp*xm*2;
    b = (p3.born - p1.born)/(xp - xm) - c*(xp + xm);
    return (a + b*x + c*x*x);
}

