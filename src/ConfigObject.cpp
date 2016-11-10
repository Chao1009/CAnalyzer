//============================================================================//
// A class based on the support from ConfigParser and ConfigValue             //
// It provides a simple way to read text file as configuration file, and read //
// or modify a parameter in the inherited class                               //
// The Configure() function should be overloaded according to specialized     //
// requirements, and be called after the parameters being configured          //
//                                                                            //
// 10/31/2016                                                                 //
//============================================================================//

#include <iostream>
#include <iomanip>
#include <algorithm>
#include "ConfigObject.h"


// constructor
ConfigObject::ConfigObject(const std::string &splitter, const std::string &ignore)
: split_chars(splitter), ignore_chars(ignore), __empty_value("")
{
    // place holder
}

// destructor
ConfigObject::~ConfigObject()
{
    // place holder
}

// configure the cluster method
void ConfigObject::Configure(const std::string & /*path*/)
{
    // to be overloaded
}

const ConfigValue &ConfigObject::GetConfigValue(const std::string &var_name)
const
{
    // convert to lower case and remove uninterested characters
    std::string key = ConfigParser::str_lower(ConfigParser::str_remove(var_name, ignore_chars));

    auto it = config_map.find(key);
    if(it == config_map.end())
        return __empty_value;
    else
        return it->second;
}

void ConfigObject::SetConfigValue(const std::string &var_name, const ConfigValue &c_value)
{
    // convert to lower case and remove uninterested characters
    std::string key = ConfigParser::str_lower(ConfigParser::str_remove(var_name, ignore_chars));

    config_map[key] = c_value;
}

void ConfigObject::SaveConfig(const std::string &path)
{
    std::string save_path;

    if(path.empty())
        save_path = config_path;
    if(save_path.empty())
        save_path = "new.conf";

    std::ofstream save(save_path);
    for(auto &it : config_map)
    {
        save << it.first
             << " = "
             << it.second
             << std::endl;
    }
}

// read configuration file and build the configuration map
void ConfigObject::readConfigFile(const std::string &path)
{
    ConfigParser c_parser(split_chars); // self-defined splitters

    if(!c_parser.OpenFile(path)) {
        std::cerr << "PRad HyCal Cluster Error: Cannot open file "
                  << "\"" << path << "\""
                  << std::endl;
        return;
    }

    // save the path
    config_path = path;

    // clear the map
    config_map.clear();

    while(c_parser.ParseLine())
    {
        if (c_parser.NbofElements() != 2)
            continue;

        std::string var_name, key;
        ConfigValue var_value;
        c_parser >> var_name >> var_value;
        // convert to lower case and remove uninterested characters
        key = ConfigParser::str_lower(ConfigParser::str_remove(var_name, ignore_chars));
        config_map[key] = var_value;
    }
}

// get configuration value from the map
ConfigValue ConfigObject::getConfigValue(const std::string &name,
                                         const ConfigValue &def_value,
                                         bool verbose)
{
    std::string key = ConfigParser::str_lower(ConfigParser::str_remove(name, ignore_chars));

    auto it = config_map.find(key);
    if(it == config_map.end())
    {
        if(def_value.IsEmpty())
            return __empty_value;

        if(verbose) {
            std::cout << name << " (key: " << key << ")"
                      << " not defined in configuration file, set to default value "
                      << def_value
                      << std::endl;
        }
        config_map[key] = def_value;
        return def_value;
    }

    return it->second;
}

