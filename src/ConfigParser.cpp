//============================================================================//
// A simple parser class to read text file                                    //
//                                                                            //
// Chao Peng                                                                  //
// 06/07/2016                                                                 //
//============================================================================//

#include "ConfigParser.h"
#include <cstring>
#include <climits>
#include <algorithm>

using namespace std;

//============================================================================//
// Config Parser                                                              //
//============================================================================//
ConfigParser::ConfigParser(const string &s,
                           const string &w,
                           const vector<string> &c)
: splitters(s), white_space(w), comment_marks(c), line_number(0)
{

}

ConfigParser::~ConfigParser()
{
    CloseFile();
}

void ConfigParser::AddCommentMark(const string &c)
{
    auto it = find(comment_marks.begin(), comment_marks.end(), c);
    if(it == comment_marks.end())
        comment_marks.push_back(c);
}

void ConfigParser::RemoveCommentMark(const string &c)
{
    auto it = find(comment_marks.begin(), comment_marks.end(), c);
    if(it != comment_marks.end())
        comment_marks.erase(it);
}

void ConfigParser::EraseCommentMarks()
{
    comment_marks.clear();
}

bool ConfigParser::OpenFile(const string &path)
{
    infile.open(path);
    line_number = 0;
    return infile.is_open();
}

void ConfigParser::CloseFile()
{
    infile.close();
}

void ConfigParser::OpenBuffer(char *buf)
{
    CloseFile(); // close file first
    line_number = 0;

    string buffer = buf;

    string line;
    for(auto c : buffer)
    {
        if(c != '\n') {
            line += c;
        } else {
            lines.push(line);
            line = "";
        }
    }
}

void ConfigParser::ClearBuffer()
{
    queue<string>().swap(lines);
}

string ConfigParser::TakeLine()
{
    if(lines.size()) {
        string out = lines.front();
        lines.pop();
        return out;
    }

    return "";
}

bool ConfigParser::ParseLine()
{
    queue<string>().swap(elements);

    if(infile.is_open()) { // if file is open, parse file
        while(elements.empty())
        {
            if(!getline(infile, current_line))
                return false; // end of file

            ParseLine(current_line);
        }
    } else { // if file is not closed, parse buffer
        while(elements.empty())
        {
            if(lines.empty())
                return false; // end of buffer
            current_line = std::move(lines.front());
            lines.pop();
            ParseLine(current_line);
        }
    }
    return true; // parsed a line
}

void ConfigParser::ParseLine(const string &line)
{
    ++line_number;
    string trim_line = trim(comment_out(line), white_space);
    queue<string> eles = split(trim_line, splitters);

    while(eles.size())
    {
        string ele = trim(eles.front(), white_space);
        if(ele.size())
            elements.push(ele);
        eles.pop();
    }
}

ConfigValue ConfigParser::TakeFirst()
{
    if(elements.empty()) {
        cout << "Config Parser: WARNING, trying to take elements while there is nothing, 0 value returned." << endl;
        return ConfigValue("0");
    }

    string output = elements.front();
    elements.pop();

    return ConfigValue(output);
}

vector<ConfigValue> ConfigParser::TakeAll()
{
    vector<ConfigValue> output;

    while(elements.size())
    {
        output.push_back(elements.front());
        elements.pop();
    }

    return output;
}

string ConfigParser::comment_out(const string &str, size_t index)
{
    if(index >= comment_marks.size())
        return str;
    else {
        string str_co = comment_out(str, comment_marks.at(index));
        return comment_out(str_co, ++index);
    }
}

string ConfigParser::comment_out(const string &str, const string &c)
{
    if(c.empty()) {
        return str;
    }

    const auto commentBegin = str.find(c);
    return str.substr(0, commentBegin);
}


string ConfigParser::trim(const string &str, const string &w)
{

    const auto strBegin = str.find_first_not_of(w);
    if (strBegin == string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(w);

    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

queue<string> ConfigParser::split(const string &str, const string &s)
{
    queue<string> eles;

    char *cstr = new char[str.length() + 1];

    strcpy(cstr, str.c_str());

    char *pch = strtok(&cstr[0], s.c_str());
    string ele;

    while(pch != nullptr)
    {
        ele = pch;
        eles.push(pch);
        pch = strtok(nullptr, s.c_str());
    }

    delete cstr;

    return eles;
}

int ConfigParser::find_integer(const string &str, const size_t &pos)
{
    vector<int> integers = find_integers(str);
    if(pos >= integers.size())
    {
        cerr << "Config Parser: Cannot find " << pos + 1 << " integers from "
             << "\"" << str << "\"."
             << endl;
        return 0;
    }

    return integers.at(pos);
}

vector<int> ConfigParser::find_integers(const string &str)
{
    vector<int> result;

    find_integer_helper(str, result);

    return result;
}

void ConfigParser::find_integer_helper(const string &str, vector<int> &result)
{
   if(str.empty())
       return;

   int negative = 1;
   auto numBeg = str.find_first_of("-0123456789");
   if(numBeg == string::npos)
       return;

   // check negative sign
   string str2 = str.substr(numBeg);

   if(str2.at(0) == '-')
   {
       negative = -1;
       int num_check;

       do {
           str2.erase(0, 1);

           if(str2.empty())
               return;

           num_check = str2.at(0) - '0';
       } while (num_check > 9 || num_check < 0);
   }

   auto numEnd = str2.find_first_not_of("0123456789");
   if(numEnd == string::npos)
       numEnd = str2.size();

   int num = 0;
   size_t i = 0;

   for(; i < numEnd; ++i)
   {
       if( (num > INT_MAX/10) ||
           (num == INT_MAX/10 && ((str2.at(i) - '0') > (INT_MAX - num*10))) )
       {
           ++i;
           break;
       }

       num = num*10 + str2.at(i) - '0';
   }

   result.push_back(negative*num);
   find_integer_helper(str2.substr(i), result);
}

// return the lower case of this string
std::string ConfigParser::str_lower(const std::string &str)
{
    string res = str;
    for(auto &c : res)
    {
        c = tolower(c);
    }
    return res;
}

// return the upper case of this string
std::string ConfigParser::str_upper(const std::string &str)
{
    string res = str;
    for(auto &c : res)
    {
        c = toupper(c);
    }
    return res;
}

// remove characters in ignore list
string ConfigParser::str_remove(const string &str, const string &ignore)
{
    string res = str;

    for(auto &c : ignore)
    {
        res.erase(remove(res.begin(), res.end(), c), res.end());
    }
    return res;
}

// replace characters in the list with certain char
string ConfigParser::str_replace(const string &str, const string &list, const char &rc)
{
    if(list.empty())
        return str;

    string res = str;

    for(auto &c : res)
    {
        if(list.find(c) != string::npos)
            c = rc;
    }

    return res;
}

// compare two strings, can be case insensitive
bool ConfigParser::strcmp_case_insensitive(const string &str1, const string &str2)
{
    if(str1.size() != str2.size()) {
        return false;
    }

    for(auto c1 = str1.begin(), c2 = str2.begin(); c1 != str1.end(); ++c1, ++c2)
    {
        if(tolower(*c1) != tolower(*c2)) {
            return false;
        }
    }

    return true;
}

//============================================================================//
// trivial funcs                                                              //
//============================================================================//

ConfigParser &operator >> (ConfigParser &c, std::string &v)
{
    v = c.TakeFirst().String();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, bool &v)
{
    v = c.TakeFirst().Bool();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, char &v)
{
    v = c.TakeFirst().Char();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, unsigned char &v)
{
    v = c.TakeFirst().UChar();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, short &v)
{
    v = c.TakeFirst().Short();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, unsigned short &v)
{
    v = c.TakeFirst().UShort();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, int &v)
{
    v = c.TakeFirst().Int();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, unsigned int &v)
{
    v = c.TakeFirst().UInt();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, long &v)
{
    v = c.TakeFirst().Long();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, unsigned long &v)
{
    v = c.TakeFirst().ULong();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, long long &v)
{
    v = c.TakeFirst().LongLong();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, unsigned long long &v)
{
    v = c.TakeFirst().ULongLong();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, float &v)
{
    v = c.TakeFirst().Float();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, double &v)
{
    v = c.TakeFirst().Double();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, long double &v)
{
    v = c.TakeFirst().LongDouble();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, const char *&v)
{
    v = c.TakeFirst().c_str();
    return c;
}

ConfigParser &operator >> (ConfigParser &c, ConfigValue &v)
{
    v = c.TakeFirst();
    return c;
}

//============================================================================//
// Config Value                                                               //
//============================================================================//
ostream &operator << (ostream &os, const ConfigValue &b)
{
    return  os << b._value;
};

ConfigValue::ConfigValue(const string &value)
: _value(value)
{}

ConfigValue::ConfigValue(const bool &value)
{
    if(value)
        _value = "1";
    else
        _value = "0";
}

ConfigValue::ConfigValue(const int &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const long &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const long long &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const unsigned &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const unsigned long &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const unsigned long long &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const float &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const double &value)
: _value(to_string(value))
{}

ConfigValue::ConfigValue(const long double &value)
: _value(to_string(value))
{}

bool ConfigValue::Bool()
const
{
    if((_value == "1") ||
       (ConfigParser::strcmp_case_insensitive(_value, "T")) ||
       (ConfigParser::strcmp_case_insensitive(_value, "True")))
        return true;

    if((_value == "0") ||
       (ConfigParser::strcmp_case_insensitive(_value, "F")) ||
       (ConfigParser::strcmp_case_insensitive(_value, "False")))
        return false;

    cout << "Config Value: Failed to convert "
         << _value << " to bool type. Return false."
         << endl;
    return false;
}

char ConfigValue::Char()
const
{
    try {
       int value = stoi(_value);
       if(value > CHAR_MAX)
           cout << "Config Value: Limit exceeded while converting "
                << _value << " to char." << endl;
       return (char) value;
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to char. 0 returned." << endl;
             return 0;
    }
}

unsigned char ConfigValue::UChar()
const
{
    try {
       unsigned long value = stoul(_value);
       if(value > UCHAR_MAX)
           cout << "Config Value: Limit exceeded while converting "
                << _value << " to unsigned char." << endl;
       return (unsigned char) value;
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to unsigned char. 0 returned." << endl;
             return 0;
    }
}

short ConfigValue::Short()
const
{
    try {
       int value = stoi(_value);
       if(value > SHRT_MAX)
           cout << "Config Value: Limit exceeded while converting "
                << _value << " to short." << endl;
       return (short) value;
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to short. 0 returned." << endl;
             return 0;
    }
}

unsigned short ConfigValue::UShort()
const
{
    try {
       unsigned long value = stoul(_value);
       if(value > USHRT_MAX)
           cout << "Config Value: Limit exceeded while converting "
                << _value << " to unsigned short." << endl;
       return (unsigned short) value;
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to unsigned short. 0 returned." << endl;
             return 0;
    }
}

int ConfigValue::Int()
const
{
    try {
       return stoi(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to int. 0 returned." << endl;
             return 0;
    }
}

unsigned int ConfigValue::UInt()
const
{
    try {
        return (unsigned int)stoul(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to unsigned int. 0 returned." << endl;
             return 0;
    }
}

long ConfigValue::Long()
const
{
    try {
        return stol(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to long. 0 returned." << endl;
             return 0;
    }
}

long long ConfigValue::LongLong()
const
{
    try {
        return stoll(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to long long. 0 returned." << endl;
             return 0;
    }
}

unsigned long ConfigValue::ULong()
const
{
    try {
        return stoul(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to unsigned long. 0 returned." << endl;
             return 0;
    }
}

unsigned long long ConfigValue::ULongLong()
const
{
    try {
        return stoull(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to unsigned long long. 0 returned." << endl;
             return 0;
    }
}

float ConfigValue::Float()
const
{
    try {
        return stof(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to float. 0 returned." << endl;
             return 0;
    }
}

double ConfigValue::Double()
const
{
    try {
        return stod(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to double. 0 returned." << endl;
             return 0;
    }
}

long double ConfigValue::LongDouble()
const
{
    try {
        return stold(_value);
    } catch (exception &e) {
        cerr << e.what() << endl;
        cerr << "Config Value: Failed to convert "
             << _value << " to long double. 0 returned." << endl;
             return 0;
    }
}

const char *ConfigValue::c_str()
const
{
    return _value.c_str();
}

