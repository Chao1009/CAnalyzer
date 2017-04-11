#ifndef RTAILS_CWRAPPER_H
#define RTAILS_CWRAPPER_H

// a C++ wrapper for the fortran subroutines in rtails.f
extern "C"
{
    // initialize rtails, must be called before using rtails
    // xi : use user defined XI or not
    // xib : user defined XI before, has no effect if xi is false
    // xia : user defined XI after, has no effect if xi is false
    // pol : polarized or unpolarized
    // theta : polarization angle, only accepts 180 deg (para) or 270 deg (perp)
    void rtails_init(bool xi, double xib, double xia, bool pol, double theta);

    // call the rtails to get the cross sections at specific kinematics
    // es : initial energy in MeV
    // ep : final energy in MeV
    // ang : scattering angle in degree
    // rlin : radiation length before
    // rlout : radiation length after
    // sigrad : OUTPUT, cross section in nb/MeV/sr
    void rtails_rad_cxsn(double es, double ep, double ang, double rlin, double rlout,
                         double *sigrad);
}

#endif
