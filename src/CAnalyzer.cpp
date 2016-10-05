#include "CAnalyzer.h"

CAnalyzer::CAnalyzer()
{
}

CAnalyzer::~CAnalyzer()
{
}

void CAnalyzer::ReadData(const std::string &path,
                         const int &cri_col,
                         const double &cri_val)
{
    Erase();
    ReadData(path, columns, cri_col, cri_val);
}

void CAnalyzer::ReadData(const std::string &path,
                         std::vector< std::vector<double> > &cols,
                         const int &cri_col,
                         const double &cri_val)
{
    ConfigParser c_parser;

    if(!c_parser.OpenFile(path))
    {
        std::cout << "Cannot open file " << path << std::endl;
    }

    while(c_parser.ParseLine())
    {
        int eles = c_parser.NbofElements();

        if(!eles)
            continue;

        if(eles > (int)cols.size())
        {
            int newcols = eles - (int)cols.size();

            size_t fillups = 0;
            if(cols.size() > 0)
                fillups = cols.at(0).size();

            for(int i = 0; i < newcols; ++i)
            {
                std::vector<double> newcol(fillups, 0.);
                cols.push_back(newcol);
            }
        }

        auto row = c_parser.TakeAll();

        // check data taking criteria
        if(cri_col >= 0 &&
           cri_col < (int)row.size() &&
           row.at(cri_col).Double() != cri_val)
            continue;

        for(size_t i = 0; i < row.size(); ++i)
        {
            cols.at(i).push_back(row.at(i).Double());
        }
    }

    c_parser.CloseFile();
}
