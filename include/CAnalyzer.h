#ifndef C_ANALYZER_H
#define C_ANALYZER_H

#include "ConfigParser.h"
#include "CMatrix.h"
#include "CEstimator.h"
#include "CRadCorr.h"
#include "canalib.h"
#include <vector>
#include <string>

class CAnalyzer
{
public:
    CAnalyzer();
    virtual ~CAnalyzer();

    void ReadData(const std::string &path,
                  const int &cri_col = -1,
                  const double &cri_val = 0,
                  const double &diff = 0);
    void ReadData(const std::string &path,
                  std::vector< std::vector<double> > &cols,
                  const int &cri_col = -1,
                  const double &cri_val = 0,
                  const double &diff = 0);
    size_t NbofCols() {return columns.size();};
    std::vector<double> GetColumn(int i) {return columns.at(i);};
    void Erase() {std::vector< std::vector<double> > empty; empty.swap(columns);};

private:
    std::vector< std::vector<double> > columns;
};

#endif
