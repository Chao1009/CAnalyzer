#ifndef CONFIG_PARSER_H
#define CONFIG_PARSER_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <queue>
#include <list>
#include <utility>
#include <typeinfo>

// config value class
class ConfigValue;

// config parser class
class ConfigParser
{
public:
    ConfigParser(const std::string &s = " ,\t",  // splitters
                 const std::string &w = " \t",  // white_space
                 const std::vector<std::string> &c = {"#", "//"}, // comment_mark
                 const std::pair<std::string, std::string> &p = std::make_pair("/*", "*/")); // comment_pair
    virtual ~ConfigParser();
    void SetSplitters(const std::string &s) {splitters = s;};
    void SetWhiteSpace(const std::string &w) {white_space = w;};
    void SetCommentMarks(const std::vector<std::string> &c) {comment_marks = c;};
    void SetCommentPair(const std::string &o, const std::string &c) {comment_pair = std::make_pair(o, c);};
    void AddCommentMark(const std::string &c);
    void RemoveCommentMark(const std::string &c);
    void EraseCommentMarks();

    bool OpenFile(const std::string &path);
    bool ReadFile(const std::string &path);
    void ReadBuffer(const char *);
    void CloseFile();
    void Clear();

    bool ParseLine();
    void ParseLine(const std::string &line);

    int NbofElements() const {return elements.size();};
    int LineNumber() const {return line_number;};
    const std::string &CurrentLine() const {return current_line;};
    std::string TakeLine();
    ConfigValue TakeFirst();
    std::vector<ConfigValue> TakeAll();

    const std::string &GetSplitters() const {return splitters;};
    const std::string &GetWhiteSpace() const {return white_space;};
    const std::vector<std::string> &GetCommentMarks() const {return comment_marks;};
    const std::pair<std::string, std::string> &GetCommentPair() const {return comment_pair;};

private:
    std::string splitters;
    std::string white_space;
    std::vector<std::string> comment_marks;
    std::pair<std::string, std::string> comment_pair;
    std::queue<std::string> lines;
    std::string current_line;
    int line_number;
    bool in_comment_pair;
    std::queue<std::string> elements;
    std::ifstream infile;

private:
    void buffer_process(std::string &buffer);
    bool parse_file();
    bool parse_buffer();
    std::string comment_out(const std::string &str, size_t index = 0);
    bool comment_between(std::string &str, const std::string &open, const std::string &close);

public:
    static std::string comment_out(const std::string &str, const std::string &c);
    static std::string trim(const std::string &str, const std::string &w);
    static std::queue<std::string> split(const std::string &str, const std::string &s);
    static std::string str_remove(const std::string &str, const std::string &ignore);
    static std::string str_replace(const std::string &str, const std::string &ignore, const char &rc = ' ');
    static std::string str_lower(const std::string &str);
    static std::string str_upper(const std::string &str);
    static std::vector<std::pair<int, int>> find_pair(const std::string &str,
                                                      const std::string &open,
                                                      const std::string &close);
    static bool strcmp_case_insensitive(const std::string &str1, const std::string &str2);
    static int find_integer(const std::string &str, const size_t &pos = 0);
    static std::vector<int> find_integers(const std::string &str);
    static void find_integer_helper(const std::string &str, std::vector<int> &result);
};

ConfigParser &operator >> (ConfigParser &c, bool &v);
ConfigParser &operator >> (ConfigParser &c, std::string &v);
ConfigParser &operator >> (ConfigParser &c, char &v);
ConfigParser &operator >> (ConfigParser &c, unsigned char &v);
ConfigParser &operator >> (ConfigParser &c, short &v);
ConfigParser &operator >> (ConfigParser &c, unsigned short &v);
ConfigParser &operator >> (ConfigParser &c, int &v);
ConfigParser &operator >> (ConfigParser &c, unsigned int &v);
ConfigParser &operator >> (ConfigParser &c, long &v);
ConfigParser &operator >> (ConfigParser &c, unsigned long &v);
ConfigParser &operator >> (ConfigParser &c, long long &v);
ConfigParser &operator >> (ConfigParser &c, unsigned long long &v);
ConfigParser &operator >> (ConfigParser &c, float &v);
ConfigParser &operator >> (ConfigParser &c, double &v);
ConfigParser &operator >> (ConfigParser &c, long double &v);
ConfigParser &operator >> (ConfigParser &c, const char *&v);
ConfigParser &operator >> (ConfigParser &c, ConfigValue &v);

// demangle type name
#ifdef __GNUG__
#include <cstdlib>
#include <memory>
#include <cxxabi.h>
// gnu compiler needs to demangle type info
static std::string demangle(const char* name)
{

    int status = 0;

    //enable c++11 by passing the flag -std=c++11 to g++
    std::unique_ptr<char, void(*)(void*)> res {
        abi::__cxa_demangle(name, NULL, NULL, &status),
        std::free
    };

    return (status==0) ? res.get() : name ;
}
#else
// do nothing if not gnu compiler
static std::string demangle(const char* name)
{
    return name;
}
#endif

// this helps template specialization in class
template <typename T>
struct __cv_identifier { typedef T type; };

class ConfigValue
{
public:
    std::string _value;

    ConfigValue() {};

    ConfigValue(const std::string &value);
    ConfigValue(const bool &value);
    ConfigValue(const int &value);
    ConfigValue(const long &value);
    ConfigValue(const long long &value);
    ConfigValue(const unsigned &value);
    ConfigValue(const unsigned long &value);
    ConfigValue(const unsigned long long &value);
    ConfigValue(const float &value);
    ConfigValue(const double &value);
    ConfigValue(const long double &value);

    template<typename T>
    T Convert()
    const
    {
        return convert( __cv_identifier<T>());
    };

    bool Bool() const;
    char Char() const;
    unsigned char UChar() const;
    short Short() const;
    unsigned short UShort() const;
    int Int() const;
    unsigned int UInt() const;
    long Long() const;
    long long LongLong() const;
    unsigned long ULong() const;
    unsigned long long ULongLong() const;
    float Float() const;
    double Double() const;
    long double LongDouble() const;
    const char *c_str() const;
    std::string String() const {return _value;};
    bool IsEmpty() const {return _value.empty();};

    operator std::string() const
    {
        return _value;
    };

    bool operator ==(const std::string &rhs) const
    {
        return _value == rhs;
    }

private:
    template<typename T>
    T convert(__cv_identifier<T> &&)
    const
    {
        std::stringstream iss(_value);
        T _cvalue;

        if(!(iss >> _cvalue)) {
            std::cerr << "Config Value Warning: Undefined value returned, failed to convert "
                      <<  _value
                      << " to "
                      << demangle(typeid(T).name())
                      << std::endl;
        }

        return _cvalue;
    }

    bool convert(__cv_identifier<bool> &&)
    const
    {
        return (*this).Bool();
    }

    std::string convert(__cv_identifier<std::string> &&)
    const
    {
        return (*this)._value;
    }
};

// show string content of the config value to ostream
std::ostream &operator << (std::ostream &os, const ConfigValue &b);

#endif
