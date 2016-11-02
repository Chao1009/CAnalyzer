#ifndef CONFIG_OBJECT_H
#define CONFIG_OBJECT_H

#include <string>
#include <unordered_map>
#include "ConfigParser.h"

class ConfigObject
{
public:
    ConfigObject(const std::string &spliiter = ":=", const std::string &ignore = " _\t");
    virtual ~ConfigObject();

    void SetConfigValue(const std::string &var_name, const ConfigValue &c_value);
    void SetIgnoreChars(const std::string &ignore) {ignore_chars = ignore;};
    void SetSplitChars(const std::string &splitter) {split_chars = splitter;};
    void SaveConfig(const std::string &path = "");

    const ConfigValue &GetConfigValue(const std::string &var_name) const;
    const std::string &GetConfigPath() const {return config_path;};
    const std::string &GetSplitChars() const {return split_chars;};
    const std::string &GetSpaceChars() const {return ignore_chars;};

    template<typename T>
    T GetConfig(const std::string &var_name)
    const
    {
        return GetConfigValue(var_name).Convert<T>();
    }

    // functions that to be overloaded
    virtual void Configure(const std::string &path = "");

protected:
    void readConfigFile(const std::string &path);
    ConfigValue getConfigValue(const std::string &var_name,
                               const ConfigValue &def_value,
                               bool verbose = true);
    template<typename T>
    T getConfig(const std::string &var_name, const T &val, bool verbose = true)
    {
        return getConfigValue(var_name, ConfigValue(val), verbose).Convert<T>();
    }

protected:
    std::string split_chars;
    std::string ignore_chars;
    std::string config_path;
    std::unordered_map<std::string, ConfigValue> config_map;

    // return this reference when there is no value found in the map
    const ConfigValue __empty_value;
};

#endif
