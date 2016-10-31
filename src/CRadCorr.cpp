#include "CRadCorr.h"
#include "ConfigParser.h"
#include <iostream>

#define ALPHA 7.297352568E-3    // 1./137.03599911
#define PI 3.1415926535897932   // pi
#define ELECM 0.510998918       // MeV
#define HBARC 197.326968        // hbar*c (MeV*fm)
#define AMUMEV 931.494043       // MeV per amu

CRadCorr::CRadCorr()
{
    // place holder
}

CRadCorr::~CRadCorr()
{
    // place holder
}

void CRadCorr::ReadExpData(const std::string &path)
{
    ConfigParser c_parser;
    c_parser.OpenFile(path);

}

#define __GAMMA_NUM 9
#define __GAMMA_G 7.0
static double __gamma_c[] = {0.9999999999998099,
                             6.7652036812188510E2,
                            -1.2591392167224028E3,
                             7.7132342877765313E2,
                            -1.7661502916214059E2,
                             1.2507343278686905E1,
                            -1.3857109526572012E-1,
                             9.9843695780195716E-6,
                             1.5056327351493116E-7};

// gamma function
double CRadCorr::gamma(const double &z)
{

    if(z < 1) {
        return PI*(1 - z)/sin(PI*(1 - z))/gamma(2 - z);
    } else if(z == 1) {
        return 1;
    } else {
        double ag = __gamma_c[0];
        for(int k = 1; k < __GAMMA_NUM; ++k)
            ag += __gamma_c[k]/(z - 1. + k);

        double output = 0.5*log(2*PI)
                        + (z - 0.5)*log(z - 0.5 + __GAMMA_G)
                        - (z - 0.5 + __GAMMA_G)
                        + log(ag);
        return exp(output);
    }
}

#define __SPENCE_NUM 8
#define __SPENCE_NMAX 50
#define __SPENCE_NMAX_PREC 1000
static double __spence_c[] = {-1.1741940560772957946600E-1,
                              -2.7618966846029390643791E-2,
                              -8.0493987190845793511240E-3,
                              -2.7095568666150792944136E-3,
                              -7.1455906877666711465857E-4,
                               4.1757495974272487715417E-2,
                              -4.9028996486663818655E-2,
                              -9.08073640732783360};

// spence function
double CRadCorr::spence(const double &z, const double &res)
{
    if(z > 1) {
        return 2*PI*PI/6 - log(z)*log(z)/2 - spence(1/z, res);
    } else if (z == 1) {
        return PI*PI/6;
    } else if (z > 0.5) {
        return PI*PI/6 - log(1-z)*log(z) - spence(1-z, res);
    } else if (z == 0.5) {
        return PI*PI/6/2 - log(0.5)*log(0.5)/2;
    } else if (z > 0) {
        return spence_tr(z, res, __SPENCE_NMAX); // do nothing, fall into the bottom session
    } else if (z == 0) {
        return 0;
    } else if (z > -0.95) {
        return spence_tr(z, res, __SPENCE_NMAX_PREC);
    } else if (z > -1.05) {
        // poly fit
        double output = 0;
        double dz = z + 1;
        for(int i = 0; i < __SPENCE_NUM; ++i)
            output += __spence_c[i]*pow(dz, i);
        return -(1 + output*dz*dz)*PI*PI/6/2 + dz*log(2);
    } else {
        return -PI*PI/6 - log(-z)*log(-z)/2 - spence(1/z, res);
    }
}

// truncation of spence function
double CRadCorr::spence_tr(const double &z, const double &res, const int &nmax)
{
    // calculate spence until res is reached
    double output = 0.;
    int n = 0;
    while(++n <= nmax)
    {
        double nth = pow(z, n)/n/n;
        output += nth;

        if(std::abs(nth) < res)
            break;
    }

    return output;
}
