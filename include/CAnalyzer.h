#ifndef CHAO_ANALYZER_H
#define CHAO_ANALYZER_H

#include "ConfigParser.h"
#include "CMatrix.h"
#include "CEstimator.h"
#include "CRadCorr.h"
#include <vector>
#include <string>
#include <iterator>
#include <algorithm>

class CAnalyzer
{
public:
    CAnalyzer();
    virtual ~CAnalyzer();

    void ReadData(const std::string &path,
                  const int &cri_col = -1,
                  const double &cri_val = 0);
    void ReadData(const std::string &path,
                  std::vector< std::vector<double> > &cols,
                  const int &cri_col = -1,
                  const double &cri_val = 0);
    size_t NbofCols() {return columns.size();};
    std::vector<double> GetColumn(int i) {return columns.at(i);};
    void Erase() {std::vector< std::vector<double> > empty; empty.swap(columns);};

private:
    std::vector< std::vector<double> > columns;
};

namespace cana
{
    // the function is based on c++ source code
    // it adds permutation parity track
    template<class BidirIt>
    bool permutate(BidirIt first, BidirIt last, int &parity)
    {
        if (first == last) return false;
        BidirIt i = last;
        if (first == --i) return false;

        while (true) {
            BidirIt i1, i2;

            i1 = i;
            if (*--i < *i1) {
                i2 = last;
                while (!(*i < *--i2))
                    ;
                std::iter_swap(i, i2);
                std::reverse(i1, last);
                size_t swap = std::distance(i1, last)/2 + 1;
                // odd number of swaps
                if(swap&1)  parity *= -1;
                // even number of swaps, no change needed
                return true;
            }
            if (i == first) {
                std::reverse(first, last);
                size_t swap = std::distance(first, last)/2;
                // odd number of swaps
                if(swap&1)  parity *= -1;
                // even number of swaps, no change needed
                return false;
            }
        }
    }
};

#endif
