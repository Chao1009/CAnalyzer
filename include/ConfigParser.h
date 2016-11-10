#ifndef CONFIG_PARSER_H
#define CONFIG_PARSER_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include "ConfigValue.h"

// config parser class
class ConfigParser
{
public:
    ConfigParser(const std::string &s = " ,\t",  // splitters
                 const std::string &w = " \t",  // white_space
                 const std::vector<std::string> &c = {"#", "//"}); // comment_mark
    virtual ~ConfigParser();

    // Set members
    void SetSplitters(const std::string &s) {splitters = s;};
    void SetWhiteSpace(const std::string &w) {white_space = w;};
    void SetCommentMarks(const std::vector<std::string> &c) {comment_marks = c;};
    void SetCommentPair(const std::string &o, const std::string &c) {comment_pair = std::make_pair(o, c);};
    void AddCommentMark(const std::string &c);
    void RemoveCommentMark(const std::string &c);
    void EraseCommentMarks();

    // source manipulation
    bool OpenFile(const std::string &path);
    bool ReadFile(const std::string &path);
    void ReadBuffer(const char *);
    void CloseFile();
    void Clear();

    // parse line, return false if no more line to parse
    bool ParseLine();
    // parse a string as a line, by default it will count it into line_number
    void ParseLine(const std::string &line, const bool &count = true);

    // get current parsing status
    int NbofElements() const {return elements.size();};
    int LineNumber() const {return line_number;};
    const std::string &CurrentLine() const {return current_line;};

    // take the lines/elements
    std::string TakeLine();
    ConfigValue TakeFirst();
    template<class BidirIt>
    int Take(BidirIt first, BidirIt last)
    {
        int count = 0;
        for(auto it = first; it != last; ++it, ++count)
        {
            if(elements.empty())
                break;

            *it = elements.front();
            elements.pop();
        }
        return count;
    }

    template<template<class, class> class Container>
    Container<ConfigValue, std::allocator<ConfigValue>> TakeAll()
    {
        Container<ConfigValue, std::allocator<ConfigValue>> res;
        while(elements.size())
        {
            res.emplace_back(std::move(elements.front()));
            elements.pop();
        }
        return res;
    }

    template<template<class, class> class Container, class T>
    Container<T, std::allocator<T>> TakeAll()
    {
        Container<T, std::allocator<T>> res;
        while(elements.size())
        {
            ConfigValue tmp(std::move(elements.front()));
            elements.pop();
            res.emplace_back(tmp.Convert<T>());
        }
        return res;
    }

    // get members
    const std::string &GetSplitters() const {return splitters;};
    const std::string &GetWhiteSpace() const {return white_space;};
    const std::vector<std::string> &GetCommentMarks() const {return comment_marks;};
    const std::pair<std::string, std::string> &GetCommentPair() const {return comment_pair;};

private:
    // private functions
    void buffer_process(std::string &buffer);
    bool parse_file();
    bool parse_buffer();
    std::string comment_out(const std::string &str, size_t index = 0);
    bool comment_between(std::string &str, const std::string &open, const std::string &close);

private:
    // private members
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

public:
    // static functions
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

#endif
