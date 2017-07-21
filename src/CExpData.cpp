#include "CExpData.h"
#include "ConfigParser.h"
#include "ConfigObject.h"
#include "canalib.h"
#include <fstream>
#include <iostream>
#include <iomanip>


// some inlines
inline double __Q2(double E0, double nu, double sinsq)
{
    return 4.*E0*(E0 - nu)*sinsq;
}


inline double __W(double Mt, double nu, double Qsq)
{
    return sqrt(Mt*(Mt + 2*nu) - Qsq);
}

// calculate ice radiation length, input is thickness in mm
inline double __ice_radl(double thickness)
{
    // hard coded for ice ONLY
    // http://pdg.lbl.gov/2008/AtomicNuclearProperties/HTML_PAGES/325.html
    return thickness/393.1;
}

// calculate ice collisional loss, input is thickness in mm
// it assumes beta is close to 1, so for electron beam, the beam energy should
// be far greater than electron's mass
inline double __ice_coll(double thickness)
{
    // hard coded for ice ONLY
    // See J. Singh's technical note for details
    // http://hallaweb.jlab.org/experiment/E97-110/tech/radlength_sagdhv130.pdf
    // Z/A = 0.55509 mol/g
    // a = 0.15353747 MeV*(cm^2/mol)
    // rho = 0.918 g/cm^3
    return 0.55509*0.15353747*0.918*(thickness/10.);

}



// constructor
CExpData::CExpData(const std::string &config)
{
    if(!config.empty())
        ReadConfigFile(config);
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

    if(buffer.empty()) {
        std::cerr << "Cannot open configuration file \"" << path << "\", "
                  << "or the file has no contents."
                  << std::endl;
        return;
    }

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
            new_set.ReadConfig(block.content, settings);

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
void CExpData::SaveResult(const std::string &path, bool save_model)
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
        if(!save_model && dset.non_rad)
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

    // default value
    settings.radl_factor_before = 1.0;
    settings.radl_factor_after = 1.0;
    conf::update_config(conf_obj, "RadL Factor Before", settings.radl_factor_before);
    conf::update_config(conf_obj, "RadL Factor After", settings.radl_factor_after);

    // set directories
    settings.data_dir = combine_path(conf_dir, settings.data_dir);
    settings.coll_dir = combine_path(conf_dir, settings.coll_dir);
    settings.accpt_dir = combine_path(conf_dir, settings.accpt_dir);
}

// update data, sorting data points and sets
void CExpData::DataUpdate()
{
    // pre-calculate variable to save time
    sin2 = std::pow(sin(settings.angle*cana::deg2rad/2.), 2);
    // inelastic part, electron-nucleon reactions
    targetM = cana::proton_mass;

    for(auto &dset : data_sets)
    {
        // sort data points by nu
        std::sort(dset.data.begin(), dset.data.end(),
                  [] (const DataPoint &p1, const DataPoint &p2)
                  {
                      return p1.nu < p2.nu;
                  });

        // calculate the invariant mass for each point
        for(auto &p : dset.data)
        {
            p.Q2 = __Q2(dset.energy, p.nu, sin2);
            p.W = __W(targetM, p.nu, p.Q2);
        }
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
    double Qsq = __Q2(E0, E0 - Eb, sin2);
    double w = __W(targetM, E0 - Eb, Qsq);

    auto inter = cana::binary_search_interval(data_sets.begin(), data_sets.end(), E0);

    // extrapolation, may result in a huge error
    if(inter.first == data_sets.end() || inter.second == data_sets.end()) {
        if(E0 < data_sets.front().energy) {
            res = data_sets.front().Interp(w);
        } else {
            res = data_sets.back().Interp(w);
        }
    // exact match
    } else if(inter.first == inter.second) {
        res = inter.first->Interp(w);
    // interpolation between two sets
    } else {
        double E1 = inter.first->energy;
        double E2 = inter.second->energy;
        res = (inter.first->Interp(w)*(E2 - E0) + inter.second->Interp(w)*(E0 - E1))
              / (E2 - E1);
    }

    // scale by Q2
    return res/Qsq;
}

// read configuration for the data set
void CExpData::DataSet::ReadConfig(const std::string &config, const CExpData::Settings &gset)
{
    ConfigObject conf_obj;

    conf_obj.ReadConfigString(config);

    // by design, the angle should be read from the global setting
    angle = gset.angle;

    // connect data set variables to configurations
    conf::update_config(conf_obj, "Energy", energy);
    conf::update_config(conf_obj, "Radiation Length Before", radl_before);
    conf::update_config(conf_obj, "Radiation Length After", radl_after);
    conf::update_config(conf_obj, "Collisional Loss Before", coll_before);
    conf::update_config(conf_obj, "Collisional Loss After", coll_after);
    conf::update_config(conf_obj, "RC Error", error);
    conf::update_config(conf_obj, "Model", non_rad);
    conf::update_config(conf_obj, "Normalization", normalization);
    conf::update_config(conf_obj, "Data File", data_file);
    conf::update_config(conf_obj, "Data Label", data_label);
    conf::update_config(conf_obj, "Acceptance File", accpt_file);
    conf::update_config(conf_obj, "Collimator File", coll_file);
    conf::update_config(conf_obj, "Radiation Length Wall", radl_wall);

    // update radiation length from global factors
    // this is for systematics study
    radl_before *= gset.radl_factor_before;
    radl_after *= gset.radl_factor_after;
    coll_before *= gset.radl_factor_before;
    coll_after *= gset.radl_factor_after;

    // update ice thickness and change corresponding values
    double ice_before = 0., ice_after = 0.;
    conf::update_config(conf_obj, "Ice Before", ice_before);
    conf::update_config(conf_obj, "Ice After", ice_after);
    radl_before += __ice_radl(ice_before);
    radl_after += __ice_radl(ice_after);
    coll_before += __ice_coll(ice_before);
    coll_after += __ice_coll(ice_after);

    // update data files path
    data_file = combine_path(gset.data_dir, data_file);
    coll_file = combine_path(gset.coll_dir, coll_file);
    accpt_file = combine_path(gset.accpt_dir, accpt_file);
    // read data points
    ReadData(data_file, data_label);
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
        if(!label.empty() && c_parser.NbofElements() > 4) {
            std::string plabel = c_parser.TakeFirst().String();
            if(plabel != label)
                continue;
        }

        // new data points
        double nu, cxsn, stat, syst;
        c_parser >> nu >> cxsn >> stat >> syst;
        // apply normalization
        cxsn *= normalization;
        stat *= normalization;

        DataPoint new_point(nu, cxsn, stat, syst);
        // calculate ep
        new_point.Ep = energy - nu;

        data.push_back(std::move(new_point));
    }
}

// interpolation between data sets by invariant mass W
// it returns the cross section scaled by E0*(E0 - Ep)
double CExpData::DataSet::Interp(const double &w)
const
{
    // out of range
    // try to extrapolate, but it is very dangerous and may result in a big error
    if(w < data.front().W) {
        double xs = (data.size() < 3) ? 0. :
                    cana::parabolic_interp(data[0].W, data[0].born*data[0].Q2,
                                           data[1].W, data[1].born*data[0].Q2,
                                           data[2].W, data[2].born*data[0].Q2, w);
        if(xs < 0.) xs = 0.;
        return xs;
    // assuming uniform distribution for high nu region (DIS)
    } else if(w >= data.back().W) {
        return data.back().born*data.back().Q2;
    }


    // search the position of w
    auto it_pair = cana::binary_search_interval(data.begin(), data.end(), w);

    // exact matched
    if(it_pair.first == it_pair.second)
        return it_pair.first->born*it_pair.first->Q2;

    // only have 2 points, do a straight line interpolation
    if(it_pair.first == data.begin()) {
        const DataPoint &p1 = *it_pair.first;
        const DataPoint &p2 = *it_pair.second;
        return cana::linear_interp<double>(p1.W, p1.born*p1.Q2,
                                           p2.W, p2.born*p2.Q2, w);
    // 3 points parabolic fit
    } else {
        const DataPoint &p1 = *(it_pair.first - 1);
        const DataPoint &p2 = *it_pair.first;
        const DataPoint &p3 = *it_pair.second;
        return cana::parabolic_interp<double>(p1.W, p1.born*p1.Q2,
                                              p2.W, p2.born*p2.Q2,
                                              p3.W, p3.born*p3.Q2, w);
    }
}

void CExpData::DataSet::CopySettings(const CExpData::DataSet &ref)
{
    energy = ref.energy;
    radl_before = ref.radl_before;
    radl_after = ref.radl_after;
    coll_before = ref.coll_before;
    coll_after = ref.coll_after;
    error = ref.error;
    normalization = ref.normalization;
    non_rad = ref.non_rad;
}
