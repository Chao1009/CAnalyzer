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
ConfigObject::ConfigObject(const std::string &splitter,
                           const std::string &ignore)
: split_chars(splitter), ignore_chars(ignore), __empty_value("")
{
    // set default replace bracket
    replace_pair = std::make_pair("{", "}");
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

void ConfigObject::ClearConfig()
{
    config_path = "";
    config_map.clear();
}

bool ConfigObject::HasKey(const std::string &var_name)
const
{
    std::string key = ConfigParser::str_lower(ConfigParser::str_remove(var_name, ignore_chars));
    if(config_map.find(key) != config_map.end())
        return true;
    return false;
}

void ConfigObject::ListKeys()
const
{
    for(auto &it : config_map)
    {
        std::cout << it.first << std::endl;
    }
}

std::vector<std::string> ConfigObject::GetKeyList()
const
{
    std::vector<std::string> res;

    for(auto &it : config_map)
        res.push_back(it.first);

    return res;
}

void ConfigObject::SaveConfig(const std::string &path)
const
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

ConfigValue ConfigObject::GetConfigValue(const std::string &var_name)
const
{
    // convert to lower case and remove uninterested characters
    std::string key = ConfigParser::str_lower(ConfigParser::str_remove(var_name, ignore_chars));

    auto it = config_map.find(key);
    if(it == config_map.end())
        return __empty_value;
    else
        return form(it->second, replace_pair.first, replace_pair.second);
}

void ConfigObject::SetConfigValue(const std::string &var_name, const ConfigValue &c_value)
{
    // convert to lower case and remove uninterested characters
    std::string key = ConfigParser::str_lower(ConfigParser::str_remove(var_name, ignore_chars));

    config_map[key] = c_value;
}

// read configuration file and build the configuration map
void ConfigObject::readConfigFile(const std::string &path)
{
    ConfigParser c_parser(split_chars); // self-defined splitters

    if(!c_parser.ReadFile(path)) {
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

    return form(it->second, replace_pair.first, replace_pair.second);
}

ConfigValue ConfigObject::form(const std::string &input,
                               const std::string &op,
                               const std::string &cl)
const
{
    ConfigValue res(input);

    auto pairs = ConfigParser::find_pair(res._value, op, cl);

    for(auto &pair : pairs)
    {
        int beg_pos = pair.first + op.size();
        int end_pos = pair.second - cl.size();
        int size = end_pos - beg_pos + 1;
        ConfigValue val = GetConfigValue(res._value.substr(beg_pos, size));
        res._value.replace(pair.first, pair.second - pair.first + 1, val.c_str());
    }

    return res;
}
