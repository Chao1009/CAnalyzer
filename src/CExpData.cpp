#include "CExpData.h"
#include "ConfigParser.h"
#include "ConfigObject.h"
#include "canalib.h"
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

    for(auto &block : blocks)
    {
        if(ConfigParser::case_ins_equal(block.label, "SETTING"))
            globalSetting(block.content);
        if(ConfigParser::case_ins_equal(block.label, "DATASET"))
            createDataSet(block.content, path);
    }

    // initialize all data sets and sort them in energy transcendent order
    dataInit();

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

void CExpData::globalSetting(const std::string &config)
{
    ConfigObject conf_obj;
    conf_obj.ReadConfigString(config);

    if(!conf::update_config(conf_obj, "Angle", angle)) {
        std::cerr << "Error: Missing angle in global settings." << std::endl;
    }

}

// create a new data set based on configuration text
void CExpData::createDataSet(const std::string &config, const std::string &path)
{
    ConfigObject conf_obj;

    conf_obj.ReadConfigString(config);

    // create a data set with default values
    DataSet new_set;

    if(!conf::update_config(conf_obj, "Energy", new_set.energy) ||
       cana::is_in(new_set.energy, data_sets.begin(), data_sets.end())) {

        std::cerr << "Error: Missing energy information or there exists a data"
                  << "set with the same energy in configuration: " << config
                  << std::endl;
        return;
    }

    // connect data set variables to configurations
    conf::update_config(conf_obj, "Radiation Length Before", new_set.radl_before);
    conf::update_config(conf_obj, "Radiation Length After", new_set.radl_after);
    conf::update_config(conf_obj, "Collisional Loss Before", new_set.coll_before);
    conf::update_config(conf_obj, "Collisional Loss After", new_set.coll_after);
    conf::update_config(conf_obj, "Ice Before", new_set.ice_before);
    conf::update_config(conf_obj, "Ice After", new_set.ice_after);
    conf::update_config(conf_obj, "RC Error", new_set.error);
    conf::update_config(conf_obj, "Born Level", new_set.non_rad);
    conf::update_config(conf_obj, "Normalization", new_set.normalization);

    // read-in the data points in this set
    std::string data_file, data_label;
    conf::update_config(conf_obj, "Data File", data_file);
    conf::update_config(conf_obj, "Data Label", data_label);

    std::string dir = ConfigParser::decompose_path(path).dir;
    readData(new_set, ConfigParser::form_path(dir, data_file), data_label);

    // sort data points by nu
    std::sort(new_set.data.begin(), new_set.data.end(),
              [] (const DataPoint &p1, const DataPoint &p2)
              {
                  return p1.nu < p2.nu;
              });

    data_sets.emplace_back(std::move(new_set));
}

// read experimental data in the format line by line
// comment marks are # and //
// splitter can be tab, space and comma
// additional white spaces (space and tab) will be trimmed
// expected 4 ~ 5 columns (label can be neglected if input data_label is empty)
// *label*, nu, cxsn, stat. error, syst. error
void CExpData::readData(DataSet &dset, const std::string &path, const std::string &label)
{
    ConfigParser c_parser;

    if(!c_parser.ReadFile(path)) {
        std::cerr << "Cannot open file \"" << path << "\", no data points read "
                  << "for data set at energy = " << dset.energy
                  << std::endl;
        return;
    }

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
        cxsn *= dset.normalization;
        stat *= dset.normalization;

        DataPoint new_point(nu, cxsn, stat, syst);
        // calculate ep
        new_point.Ep = dset.energy - nu;
        new_point.v = nu/dset.energy;

        dset.data.push_back(std::move(new_point));
    }
}

void CExpData::dataInit()
{
    double theta = angle*cana::deg2rad;
    double sin2 = std::pow(sin(theta/2.), 2);

    // mott cross section
    double mott_factor = std::pow(cana::hbarc*cana::alpha*cos(theta/2)/2/sin2, 2)*1e7; // MeV^2*nb/sr

    for(auto &dset : data_sets)
    {
        dset.weight_mott = mott_factor/std::pow(dset.energy, 2);
    }

    // sort data sets by energy
    std::sort(data_sets.begin(), data_sets.end(),
             [] (const DataSet &set1, const DataSet &set2)
             {
                 return set1.energy < set2.energy;
             });
}

// interpolates or extrapolates
double CExpData::GetCrossSection(const double &E0, const double &Eb)
const
{
    double w = (E0 - Eb)/E0;
    double scale = mott_factor/E0/E0;

    auto inter = cana::binary_search_interval(data_sets.begin(), data_sets.end(), E0);

    // extrapolation
    if(inter.first == data_sets.end() || inter.second == data_sets.end()) {
        if(E0 < data_sets.front().energy)
            return data_sets.front().Interp(w)*scale;
        else
            return data_sets.back().Interp(w)*scale;
    // exact match
    } else if(inter.first == inter.second) {
        return inter.first->Interp(w)*scale;
    // interpolation between two sets
    } else {
        double E1 = inter.first->energy;
        double E2 = inter.second->energy;
        return (inter.first->Interp(w)*(E2 - E0) + inter.second->Interp(w)*(E0 - E1))/(E2 - E1)*scale;
    }
}

// interpolation between a data set according to (1 - Ep/E0)
double CExpData::DataSet::Interp(const double &w)
const
{
    // for low nu region, it is impossible to extrapolate
    if(w < data.front().v)
        return 0.;
    // assuming uniform distribution for high nu region (DIS)
    if(w >= data.back().v)
        return data.back().born/weight_mott;

    // search the position of w
    auto it_pair = cana::binary_search_interval(data.begin(), data.end(), w);

    // not found
    if(it_pair.second == data.end() || it_pair.first == data.end())
        return 0;

    // exact matched
    if(it_pair.first == it_pair.second)
        return it_pair.first->born/weight_mott;

    // only have 2 points, do a straight line interpolation
    if(it_pair.first == data.begin()) {
        const DataPoint &p1 = *it_pair.first;
        const DataPoint &p2 = *it_pair.second;

        double _interp = p1.born*(p2.v - w) + p2.born*(w - p1.v);

        // LATEST CHANGE: removed unknown constant 0.00001 here, 
        // probably some protection for two same points
        // return _interp/(p2.v - p1.v + 0.00001);
        return _interp/(p2.v - p1.v)/weight_mott;
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
    return (a + b*x + c*x*x)/weight_mott;
}

