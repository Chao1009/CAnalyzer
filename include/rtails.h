#ifndef RTAILS_CWRAPPER_H
#define RTAILS_CWRAPPER_H

// a C++ wrapper for the fortran subroutines in rtails.f
extern "C"
{
    // initialize rtails, must be called before using rtails
    // xi : use user defined XI or not
    // pol : polarized or unpolarized
    // theta : polarization angle, only accepts 180 deg (para) or 270 deg (perp)
    void rtails_init(bool xi, bool pol, double theta);

    // initialize the radiation length before get radiative cross section
    // rlb : radiation length of all materials before
    // rla : radiation length of all materials after
    // xib : user defined XI before, has no effect if xi is false
    // xia : user defined XI after, has no effect if xi is false
    void rtails_set_radl(double rlb, double rla, double xib = 0., double xia = 0.);

    // call the rtails to get the cross sections at specific kinematics
    // es : initial energy in MeV
    // ep : final energy in MeV
    // ang : scattering angle in degree
    // unpol_MT : using MT approach or BS approach on unpolarized internal radiative effects
    // output: cross section in ub/MeV/sr
    double rtails_rad_cxsn(double es, double ep, double ang, bool unpol_MT = true);
}

#endif
