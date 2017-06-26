/* KYNrefrev - ionised reflection (lamp-post Compton reflection) reverberation
 *             model
 *
 * ref. Dovciak M., Karas V., Martocchia A., Matt G., Yaqoob T. (2004)
 * --> this is very old reference on KY package of models when the reverberation
 *     code has not been developed yet but was developed from this package
 * -----------------------------------------------------------------------------
 * OTHER REFERENCES:
 * 
 * Dovciak M., Karas V., Martocchia A., Matt G. & Yaqoob T. (2004). XSPEC model
 * to explore spectral features from black hole sources.
 * In Proc. of the workshop on processes in the vicinity of black holes
 * and neutron stars. S.Hledik & Z.Stuchlik, Opava. In press. [astro-ph/0407330]
 * 
 * Dovciak M., Karas V. & Yaqoob, T. (2004). An extended scheme
 * for fitting X-ray data with accretion disk spectra
 * in the strong gravity regime. ApJS, 153, 205.
 * 
 * Dovciak M. (2004). Radiation of accretion discs in strong gravity. Faculty of
 * Mathematics and Physics, Charles University, Prague. PhD thesis.
 * [astro-ph/0411605]
 * -----------------------------------------------------------------------------
 * 
 * This code computes the emission from an acrretion disc that is
 * illuminated from the primary power-law source located on the axis above the
 * central black hole with a flash. All relativistic effects are taken into 
 * account (in all three parts of the light path - from the primary source to 
 * the observer, to disc and from the disc to the observer). The code calls 
 * subroutine ide() for integrating local emission over the disc and uses the 
 * fits file 'KBHtables80.fits' defining the transfer functions needed for the 
 * integration. For details on ide() and the fits file see the subroutine ide(). 
 * This model also uses KBHlamp_qt.fits file with transfer functions for the 
 * light coming from primary source and hitting the accretion disc (see the 
 * description of this file below). The reflection is taken from Ross & Fabian 
 * tables reflionx.mod.
 *
 * par1  ... a/M     - black hole angular momentum (-1 <= a/M <= 1)
 * par2  ... theta_o - observer inclination in degrees (0-pole, 90-disc)
 * par3  ... rin - inner edge of non-zero disc emissivity (in GM/c^2 or in 
 *                 r_mso)
 * par4  ... ms  - 0 - we integrate from inner edge = par3 
 *                 1 - if the inner edge of the disc is below marginally stable
 *                     orbit then we integrate emission above MSO only
 *                 2 - we integrate from inner edge given in units of MSO, i.e.
 *                     inner edge = par3 * r_mso (the same applies for outer 
 *                     edge)
 * par5  ... rout  - outer edge of non-zero disc emissivity (in GM/c^2 or in 
 *                   r_mso)
 * par6  ... phi   - lower azimuth of non-zero disc emissivity (deg)
 * par7  ... dphi  - (phi + dphi) is upper azimuth of non-zero disc emissivity
 *                   0 <= dphi <= 360  (deg)
 * par8 ... M/M8   - black hole mass in units of 10^8 solar masses
 * par9 ... height - height on the axis (measured from the center) at which
 *                   the primary source is located (GM/c^2)
 * par10 ... PhoIndex - power-law energy index of the primary flux
 * par11 ... L/Ledd   - dE/dt, the intrinsic local (if negative) or the 
 *                      observed (if positive) primary isotropic flux in the 
 *                      X-ray energy range 2-10keV in units of Ledd
 * par12 ... Np:Nr  - ratio of the primary to the reflected normalization
 *                    1 - self-consistent model for isotropic primary source
 *                    0 - only reflection, primary source is hidden
 *                  - if negative then L/Ledd (par11) means the luminosity 
 *                    towards the disc
 *                  - if positive then L/Ledd (par11) means the luminosity 
 *                    towards the observer
 * par13 ... density  - density profile normalization in 10^15 cm^(-3)
 * par14 ... den_prof - radial power-law density profile
 * par15 ... abun    - Fe abundance (in solar abundance)
 * par16 ... therm   - fraction of thermalised flux from the overal incident 
 *                     flux illuminating the disc
 *                     = 0 - only the reverberation of reflected radiation is 
 *                           computed
 *                     < 0 - only the reverberation of thermal radiation is 
 *                           computed
 *                     > 0 - both the thermal and reflection reverberation is
 *                           included
 *                     abs(par16) > 1 - the fraction of thermalisation is 
 *                                      computed from difference between the 
 *                                      incident and reflected fluxes
 * par17 ... arate  - accretion rate in units of Ledd if positive or in 
 *                    Solar mass per Julian year (365.25days) if negative
 * par18 ... f_col  - spectral hardening factor
 * par19 ... alpha  - position of the cloud centre in GM/c^2 in alpha coordinate
 *                    (alpha being the impact parameter in phi direction, 
 *                     positive for approaching side of the disc)
 * par20 ... beta   - position of the cloud centre in GM/c^2 in beta coordinate
 *                    (beta being the impact parameter in theta direction, 
 *                     positive in up direction, i.e. away from the disc)
 * par21 ... rcloud - radius of the obscuring cloud
 *                  - the meaning of cloud is inverted for negative values of 
 *                    rcloud, i.e. only the radiation behind the cloud is 
 *                    computed
 * par22 ... zshift - overall Doppler shift
 * par23 ... limb   - limb darkening/brightening law (emission directionality)
 *                  - if = 0 the local emisivity is not multiplied by anything
 *                  - if = 1 the local emisivity is multiplied by 1+2.06*mu
 *                    (limb darkening)
 *                  - if = 2 the local emisivity is multiplied by ln(1+1/mu)
 *                    (limb brightening)
 * par24 ... tab - which reflion table to use 
 *                 1 -> reflion (the old one, lower cut-off energy at 1eV,
 *                      not good for PhoIndex > 2)
 *                 2 -> reflionx (the new one, lower cut-off energy at 100eV)
 * par25 ... sw  - swich for the way how to compute the refl. spectra
 *                 1 -> use the computed Xi for the interpolation in reflion,
 *                      i.e. use proper total incident intensity
 *                      with the shifted cut-offs
 *                 2 -> use the Xi correspondent to the computed normalization
 *                      of the incident flux, i.e. do not shift the cut-offs
 *                      when computing the total incident intensity
 * par26 ... ntable - table of relativistic transfer functions used in the model
 *                    (defines fits file with tables), 0<= ntable <= 99
 * par27 ... nrad   - number of grid points in radius
 *                  - if negative than the number of radial grid points is 
 *                    dependent on height as -nradius / height^0.66
 * par28 ... division - type of division in radial integration
 *                      0 -> equidistant radial grid (constant linear step)
 *                      1 -> exponential radial grid (constant logarithmic step)
 *                      >1 -> mixed radial grid with a constant logarithmic step
 *                            in the inner region and with a constant linear 
 *                            step in the outer region; the total nradius 
 *                            (par27) number of points is divided in the 3:2 
 *                            ratio in these regions; the value of par28 gives 
 *                            the transition radius between these regions 
 *                            (in GM/c^2)
 *                      -1 -> mixed radial grid with the transition radius at 
 *                            2*height
 * par29 ... nphi     - number of grid points in azimuth
 * par30 ... deltaT   - length of the time bin (GM/c^3)
 * par31 ... nt       - number of time subbins per one time bin
 * par32 ... t1/f1/E1 - the time to be used in xspec for the spectrum
 *                      (0 means average spectrum, i.e. divided by the 
 *                       flare duration)
 *                    - the frequency to be used in XSPEC for the energy 
 *                      dependent Fourier transform
 *                      (0 means average values in the range of 0 to the first
 *                       wrapping frequency)
 *                    - positive values are in sec or Hz
 *                    - negative values are in GM/c^3 or (GM/c^3)^(-1)
 *                    - if different than par33, the value gives the lower end
 *                      of the time/frequency interval of interest
 *                    - if same as par33, then the functions are computed for 
 *                      this value of the time/frequency of interest
 *                    - in case of frequency dependent lags it defines the lower
 *                      value of the energy band of interest in keV
 * par33 ... t2/f2/E2 - used only if different than par32 and if par32 is
 *                      nonzero
 *                    - its value gives the upper end of the time/frequency
 *                      interval of interest
 *                    - positive values are in sec or Hz
 *                    - negative values are in GM/c^3 or (GM/c^3)^(-1)
 *                    - in case of frequency dependent lags it defines the upper
 *                      value of the energy band of interest in keV
 * par34 ... Eref1   - it defines the lower value of the reference energy band
 *                     for lag or amplitude energy dependence as well as in 
 *                     case of frequency dependent lags and amplitudes
 *                   - if zero no reference band is used
 *                   - if negative:
 *                     * for lag-energy spectra, the whole energy band is used 
 *                       as a reference band, always excluding the current 
 *                       energy bin
 *                     * for lag-frequency dependence, the energy reference band 
 *                       is abs(par34) to abs(par35) excluding overlaping part 
 *                       with energy band of interest abs(par32) to abs(par33)
 * par35 ... Eref2   - it defines the upper value of the reference energy band
 *                     for lag-energy dependence as well as in case of 
 *                     frequency dependent lags
 * par36 ... dt/Af   - lag shift for lag-energy dependence in case of 
 *                     par38=+6
 *                   - multiplicative factor in case of adding empirical hard
 *                     lags Af*f^(qf), used for par38=+16 and par38=+18; 
 *                     if par36=-1 then the following hard lags prescription
 *                     is used (see Epitropakis & Papadakis, 2017):
 *                     100 * log10(Eref/E) * (f/1e-4)^(-1) s
 *                     with Eref being middle of the reference energy band
 *                     and E middle of the energy band of interest
 * par37 ... Amp/qf  - multiplicative factor for the amplitude-energy 
 *                     dependence in case of par38=+5
 *                   - powerlaw index in case of adding empirical hard 
 *                     lags Af*f^(qf), used for par38=+16 and par38=+18
 * par38 ... xsw - function to be stored in the XSPEC photar array
 *                  0 -> spectrum at time defined by par32 and par33,
 *                 the following values correspond to energy
 *                 dependent Fourier transform at the frequency band 
 *                 defined by par32 and par33:
 *                 -1 -> real part of FT of the relative reflection
 *                 -2 -> imaginary part of FT of the relative reflection
 *                 -3 -> amplitude of FT of the relative reflection
 *                 -4 -> phase of FT of the relative reflection
 *                 -5 -> amplitude  for the relative reflection
 *                       divided by amplitude in the reference energy band
 *                       defined by par34 and par35 (integration in frequencies
 *                       is done in real and imaginary parts first and then 
 *                       the amplitudes are computed)
 *                 -6 -> lag for the relative reflection with respect
 *                       to reference energy band defined by par34 and 
 *                       par35 (integration in frequencies is done in real and
 *                       imaginary parts first and then the lags are computed
 *                       with frequency at half of the wrapping frequency or 
 *                       middle of the frequency band)
 *                 -7 -> amplitude  for the relative reflection
 *                       divided by amplitude in the reference energy band
 *                       defined by par34 and par35 (integration in frequencies
 *                       here is done in amplitudes directly)
 *                 -8 -> lag for the relative reflection with respect
 *                       to reference energy band defined by par34 and 
 *                       par35 (integration in frequencies here is done in 
 *                       lags directly)
 *                  1 -> real part of FT including primary radiation
 *                  2 -> imaginary part of FT including primary radiation
 *                  3 -> amplitude of FT including primary radiation
 *                  4 -> phase of FT including primary radiation
 *                  5 -> amplitude including the primary radiation
 *                       divided by amplitude in the reference energy band
 *                       defined by par34 and par35 (integration in frequencies
 *                       is done in real and imaginary parts first and then 
 *                       the amplitudes are computed)
 *                  6 -> lag diluted by primary radiation with respect
 *                       to reference energy band defined by par34 and 
 *                       par35 (integration in frequencies is done in real and
 *                       imaginary parts first and then the lags are computed  
 *                       with frequency at half of the wrapping frequency or 
 *                       middle of the frequency band)
 *                  7 -> amplitude including the primary radiation
 *                       divided by amplitude in the reference energy band
 *                       defined by par34 and par35 (integration in frequencies
 *                       here is done in amplitudes directly)
 *                  8 -> lag diluted by primary radiation with respect
 *                       to reference energy band defined by par34 and 
 *                       par35 (integration in frequencies here is done in 
 *                       lags directly)
 *                 the following values correspond to frequency dependent
 *                 Fourier transform for the energy band of interest
 *                 defined by par32 and par33:
 *                 -11 -> real part of FT of the relative reflection
 *                 -12 -> imaginary part of FT of the relative reflection
 *                 -13 -> amplitude of FT of the relative reflection
 *                 -14 -> phase of FT of the relative reflection
 *                 -15 -> amplitude  for the relative reflection
 *                        divided by amplitude in the reference energy
 *                        band defined by par34 and par35
 *                        (rebinning here is done in real and imaginary parts 
 *                         first and then the amplitudes are computed)
 *                 -16 -> lag for the relative reflection with respect
 *                        to reference energy band defined by par34 and 
 *                        par35 (rebinning here is done in real and imaginary
 *                        parts first and then the lags are computed)
 *                 -17 -> amplitude  for the relative reflection
 *                        divided by amplitude in the reference energy
 *                        band defined by par34 and par35
 *                        (rebinning here is done in amplitudes directly)
 *                 -18 -> lag for the relative reflection with respect
 *                        to reference energy band defined by par34 and 
 *                        par35 (rebinning here is done in lags directly)
 *                  11 -> real part of FT including primary radiation
 *                  12 -> imaginary part of FT including primary radiation
 *                  13 -> amplitude of FT including primary radiation
 *                  14 -> phase of FT including primary radiation
 *                  15 -> amplitude including the primary radiation
 *                        divided by amplitude in the reference energy
 *                        band defined by par34 and par35
 *                        (rebinning here is done in real and imaginary parts 
 *                         first and then the amplitudes are computed)
 *                  16 -> lag diluted by primary radiation with respect
 *                        to reference energy band defined by par34 and 
 *                        par35 (rebinning here is done in real and imaginary
 *                        parts first and then the lags are computed)
 *                  17 -> amplitude including the primary radiation
 *                        divided by amplitude in the reference energy
 *                        band defined by par34 and par35
 *                        (rebinning here is done in amplitudes directly)
 *                  18 -> lag diluted by primary radiation with respect
 *                        to reference energy band defined by par34 and 
 *                        par35 (rebinning here is done in lags directly)
 * 
 * par39 ... nthreads - how many threads should be used for computations
 *
 * NOTES:
 *  -> accuracy vs. speed trade off depends mainly on: nradius, nphi, nt, deltaT
 *  -> the normalization of this model has to be unity!
 *
 * -------------------------------------------
 * KBHlamp_qt.fits
 *
 * This file contains pre-calculated values of the functions needed for the
 * lamp-post models. It is supposed that a primary source of emission is placed
 * on the axis at a height h from the centre above the Kerr black hole.
 * The matter in the disc rotates on stable circular (free) orbits above the
 * marginally stable orbit and it is freely falling below this orbit
 * where it has the same energy and angular momentum as the matter
 * which is on the marginally stable orbit. It is assumed that the corona
 * between the source and the disc is optically thin, therefore ray-tracing
 * in the vacuum Kerr space-time could be used for computing the functions.
 *
 * There are six functions stored in the KBHlamp_qt.fits file as columns. 
 * They are parametrized by the value of the horizon of the black hole (FITS 
 * table extensions), height of the primary source (rows) and either the
 * inclination angles of the observer (elements) or the radius in GM/c^2
 * at which a photon strikes the disc (elements). All these are defined
 * as vectors at the beginning of the FITS file. The tables are defined for
 * r_horizon = 1.00, 1.05, ..., 1.95, 2.00,
 * h-r_horizon = 0.1 - 100 (100 values with exponentially growing step),
 * inclination = 0.1, 1, 5, 10, 15, ..., 80, 85, 89 and
 * r-r_horizon = 0.01 - 1000 (100 values with exponentially growing step).
 *
 * The functions included are:
 * dWadWo - amplification of the primary source flux:
 *        = sin(theta_axis_local) / sin(theta_observer) *
 *          dtheta_axis_local/dtheta_observer
 * delay_a - delay of photon arrival time from the axis to the observer
 * q  - constant of motion defining the photon angular momentum
 * pr - the radial component of the photon momentum at the disc
 * dWadSd - part of the amplification of the incident flux on the disc from the
 *          primary source:
 *        = sin(theta_axis_local) * dtheta_axis_local / dtheta_fake, with the
 *          theta_fake defined as tan(theta_fake) = radius_incident / height
 *          one has to multiply by h/(r^2+h^2) to get the full amplification
 * delay - delay of photon arrival time from the axis to the disc
 *
 * The rest of the functions below are computed from the above ones:
 * g-factor - the ratio of the energy of a photon hitting the disc to the energy
 *            of the same photon when emitted from a primary source,
 * cosine of the incident angle - an absolute value of the cosine of the local
 *                                incident angle between the incident light ray
 *                                and local disc normal
 * azimuthal incident angle - the angle (in radians) between the projection of
 *                            the three-momentum of the incident photon into the
 *                            disc (in the local rest frame co-moving with the
 *                            disc) and the radial tetrad vector.
 *
 * The definition of the file KBHlamp_qt.fits:
 * 0. All of the extensions defined below are binary.
 * 1. The first extension contains a vector of the horizon values in GM/c^2
 *    (1.00 <= r_horizon <= 2.00).
 * 2. The second extension contains a vector of the values of heights h of
 *    a primary source in GM/c^2 (0.1 <= h-r_horizon <= 100).
 * 3. The third extension contains a vector of the values of the observer
 *    inclinations (0.1 <= inclination <= 89).
 * 4. The fourth extension contains a vector of the values of the incident
 *    radius (0.01 <= radius-r_horizon <= 1000).
 * 5. In the following extensions the functions are defined, each extension is
 *    for a particular value of r_horizon, each row is for a particular value of
 *    height and each element is either for a particular value of inclination
 *    (dWadWo, delay_a) or incident radius (the rest of the functions).
 * 6. Each of these extensions has six columns. In each column, a particular
 *    function is stored - dWadWo, delay_a, q, pr, dWadSd and delay, 
 *    respectively.
 * 
 *==============================================================================
 *
 * non-revreberation code:
 * 27. 8.2009  changed to work with new tables with azimuthal dependence
 * 12.11.2009  changed to work with new lamp-post tables, we can fit for height
 *             as well now
 * 
 * reverberation code:
 * 21. 7.2012  time variation (box flare) added
 * 30. 9.2015  code transformed into the C-language
 * 23. 2.2016  major update so that the code is suitable for fitting lag-energy
 *             as well as lag-frequency dependences inside XSPEC
 *
 ******************************************************************************/

/*******************************************************************************
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"

//definition of some parameters
#define DURATION -10. // duration of the flare (GM/c^3), 
                      // if negative, then duration = -duration * deltaT
#define TMAX     256. // maximum time for computation of light curves
                      // and dynamic spectra (for response echo)
#define NN       8192 // number of time bins for computing FFT
#define SMOOTH   0

/*******************************************************************************
*******************************************************************************/

#ifdef OUTSIDE_XSPEC

#define IFL    1
#define NPARAM 40
/*
// for the energy dependence
#define NE     15
#define E_MIN  0.3
#define E_MAX  10.
#define NBANDS 5
*/
/*
// for a nice energy dependence
#define NE     200
#define E_MIN  0.1
#define E_MAX  80.
#define NBANDS 5
*/

// for the frequency dependence
#define NE     80
#define E_MIN  3e-5
#define E_MAX  0.004
#define NBANDS 2


/* Let's declare variables that are common for the main and KYNrefrev
   subroutines */
static double ener_low[NBANDS], ener_high[NBANDS];

int main() {

void KYNrefrev(double *ear, int ne, const double *param, int ifl, 
               double *photar, double *photer, const char* init);

double ear[NE + 1], param[NPARAM], photar[NE], photer[NE];
char   initstr[0] = "";
int    ie, ia, iinc, ih;

//definition of energy band of interest, reference band is defined as the last 
//one, usually the whole energy range
ener_low[0] = 0.3;
ener_high[0] = 1.;
ener_low[1] = 1.2;
ener_high[1] = 4.;
/*
ener_low[0] = 0.3;
ener_high[0] = 0.9;
ener_low[1] = 1.;
ener_high[1] = 3.;
ener_low[2] = 3.;
ener_high[2] = 9.;
ener_low[3] = 12.;
ener_high[3] = 40.;
ener_low[4] = E_MIN;
ener_high[4] = E_MAX;
*/
//definition of the KYNrefrev parameters
param[ 0] = 1.;       // a/M
param[ 1] = 30.;      // thetaO
param[ 2] = 1.;       // rin
param[ 3] = 1.;       // ms
param[ 4] = 1000.;    // rout
param[ 5] = 0.;       // phi0
param[ 6] = 360.;     // dphi
param[ 7] = 0.1;      // M/M8
param[ 8] = 3.;       // height
param[ 9] = 2.;       // PhoIndex
param[10] = 0.001;    // L/Ledd
param[11] = 1.;       // Np:Nr
param[12] = 1.;       // density
param[13] = 0.;       // den_prof
param[14] = 1.;       // abun
param[15] = 0.;       // thermalisation
param[16] = 0.1;      // arate
param[17] = 2.4;      // f_col
param[18] = -6.;      // alpha
param[19] = 0.;       // beta
param[20] = 0.;       // rcloud
param[21] = 0.;       // zshift
param[22] = 0.;       // limb
param[23] = 2.;       // tab
param[24] = 2.;       // sw
param[25] = 80.;      // ntable
param[26] = -4488.;   // nrad
param[27] = -1.;      // division
param[28] = 180.;     // nphi
param[29] = 1.;       // deltaT
param[30] = 1.;       // nt
/*
param[31] = 0.;       // t1/f1/E1
param[32] = 8.3e-4;   // t2/f2/E2
param[33] = -1.;      // Eref1
param[34] = 3.;       // Eref2
*/
// the following should be used only for debugging purposes for the case of 
// abs(photar_sw) > 10
// the energy bands above should be then re-defined to consist of just 2 bands..
// for energy band definitions the param[] values are used, while ener_low[]
// and ener_high[] are ignored later on!

param[31] = ener_low[0];       // t1/f1/E1
param[32] = ener_high[0];      // t2/f2/E2
param[33] = ener_low[1];       // Eref1
param[34] = ener_high[1];      // Eref2

param[35] = 0.;       // dt/Af
param[36] = 1.;       // Amp/qf
param[37] = 16.;      // xsw
param[38] = 4.;       // nthreads
param[39] = 1.;       // norm

for (ie = 0; ie <= NE; ie++) {
//  ear[ie] = E_MIN + ie * (E_MAX-E_MIN) / NE;
  ear[ie] = E_MIN * pow(E_MAX / E_MIN, ((double) ie) / NE);
}

//for (ia=0;ia<=1;ia++){
for (ia=1;ia<=1;ia++){
//  param[0] = (double) ia;
//  for (iinc=30;iinc<=60;iinc+=30){
  for (iinc=30;iinc<=30;iinc+=30){
//    param[1] = (double) iinc;
//    for (ih=1;ih<=5;ih++){
    for (ih=2;ih<=2;ih++){
      if(ih==1)param[8]=1.5;
//      else if(ih==2)param[8]=3.;
      else if(ih==3)param[8]=6.;
      else if(ih==4)param[8]=15.;
      else if(ih==5)param[8]=30.;
//    for (ih=1;ih<=20;ih++){
//       param[8] = 1.5 * (100./1.5)**((ih-1.)/19.);
       if( param[8] >= (1.1+sqrt(1.-param[0]*param[0])) )
         KYNrefrev(ear, NE, param, IFL, photar, photer, initstr);
    }
  }
}

return(0);
}

#endif
/*******************************************************************************
*******************************************************************************/
//Mpc in cm
#define MPC      3.0856776e+24
#define ERG      6.241509e8
// rg is gravitational unit GM/c^2 in cm^2 for M = 10^8*Msun and rg2=rg^2
#define RG2      2.1819222e26
// if Ledd is changed one need to recalculate also logxi_norm0 below!!!
// Ledd is in erg (not W) and multiplied by 10^8 due to (M / (10^8*Msun)) scale
#define LEDD     1.26e46
#define HUBBLE   70.
//speed of light in km/s
#define CLIGHT   299792.458
// sec = G*10^8*Msun/c^3
#define SEC      492.568
// logxi_norm=alog10(Ledd*(M/(10^8*Msun))/(rg*m8)^2/nH)=
//            alog10(Ledd/RG2/1e15)+alog10(mass)=
//            alog10(1.26e46/(1.477e5*1e8)^2/1e15)+alog10(mass)
//            logxi_norm0+alog10(mass)
#define LOGXI_NORM0   4.761609554
#define PI       3.14159265359
#define PI2      6.2831853071795865
//kev in SI units (not in ergs!)
#define KEV      1.6022e-16
#define H_KEVS   4.13566743e-18
#define K_KEVK   8.6174e-8
//speed of light in cm/s
#define C_CMS    2.99792458e10
#define MSOLAR   1.989e+30
#define SIGMA    5.6704e-8
#define YEAR     31557600.0
#define LAMP     "KBHlamp_qt.fits"
#define REFLION1 "reflion.mod"
#define REFLION2 "reflionx.mod"
#define REFLIONX_NORM 1.e20
#define EC       300.

// Let's declare variables that can be seen from ide and FFT routine
double    *del=NULL, t0, fwrap=0.;
float     *radius=NULL;
long int  nradius;
int       exclude_energy=0;

/* Let's declare variables that are common for the main and emissivity 
   subroutines */
static float  *xi=NULL, *logxi=NULL;
static double *gfac=NULL, *transf_d=NULL, *energy1=NULL, *flux1=NULL;
//static double *cosin=NULL, *phiph=NULL;
static double *flx=NULL;
static double h, gamma0, nH0, qn, mass, mass2, am2, r_plus, Np, dt,
              flare_duration_rg, flare_duration_sec;
static double thermalisation, arate, x0, x1, x2, x3, Tnorm, am, f_col;
static long   nxi;
static int    sw, limb, nt_ratio;

extern char*  FGMODF(void);
extern char*  FGMSTR(char* dname);
extern void   FPMSTR(const char* value1, const char* value2);
extern int    xs_write(char* wrtstr, int idest);
//extern double incgamma(double a, double x);
//extern void   cutoffpl(double *ear, int ne, double *param, double *photar);

void KYNrefrev(double *ear_xspec, int ne_xspec, const double *param, 
               int ifl, double *photar, double *photer, const char* init) {
  
extern int ide(const double *ear, const int ne, const long nt, double *far, 
               double *qar, double *uar, double *var, 
               const double *ide_param, void (*emissivity)(), 
               const int ne_loc);

extern void fft_reverberation(const double *ear, int ne, double *photar,
                              double frequency1, double frequency2, 
                              int photar_sw, double *time, 
                              long nn, long nt, double *far, double *far_prim, 
                              double *flux_bands, double *flux_bands_prim, 
                              int nbands, double tot_time_sec, 
                              double flare_duration_sec, char *cparam);


void emis_KYNrefrev(double** ear_loc, const int ne_loc, const long nt, 
                    double *far_loc, double *qar_loc, double *uar_loc, 
                    double *var_loc, const double r, const double phi, 
                    const double cosmu, const double phiphoton, 
                    const double alpha_o, const double beta_o, 
                    const double delay, const double g);

double incgamma(double a, double x);


/* Let's declare static variables whose values should be remembered for next
   run in XSPEC */
static char   kydir[255] = "";
static char   pname[128] = "KYDIR", pkyLxLamp[128] = "KYLxLamp";
static char   pkyxiin[128] = "KYXIin", pkyxiout[128] = "KYXIout";
static char   pkyRefl[128] = "KYRefl", pkyfwrap[128] = "KYfwrap";
static long   nrh, nh, nincl, nener;
static float  *r_horizon, *height, *incl, *dWadWo, *delay_a, *q, *pr, *dWadSd, 
              *delay, *abun, *gam, *emission, *energy0;
static long   ne_loc, nabun, ngam;
static double *energy2, *flux0, transf_o, del_a;
static double h_rh_old = -1., gam_old = -1., abun_old = -1., thetaO_old = -1.,
              am_old = -1.;
static int    rix_old = -1, first_h = 1, first_rix = 1;

FILE   *fw;
char   errortxt[255], filename[255];
char   reflion[255], text[255];
char   kyxiin[32], kyxiout[32], kyLxLamp[32], kyRefl[32], kyfwrap[32];
char   cparam[12];
double ide_param[25];
double ear_short[4];
double *time=NULL;
double *ear=NULL, *far=NULL, *qar=NULL, *uar=NULL, *var=NULL;
double *far_final=NULL;
//double far_prim[ne_xspec];
//double flux_prim=0.,spectrum_prim[ne_xspec];
double *flux=NULL, *spectrum=NULL, *far_prim=NULL;
double flux_prim=0., *spectrum_prim=NULL;
#ifndef OUTSIDE_XSPEC
#define NBANDS 2
double ener_low[NBANDS], ener_high[NBANDS];
#endif
double *flux_bands=NULL, flux_bands_prim[NBANDS];
int    ne=0., nbands=0;
long   nn=NN, ntmax=0;
double en1=0., en2=0., en3=0., en4=0.;
double lag_shift=0., ampl_ampl=1.;
double ttmp, ttmp1, utmp, utmp1, vtmp, vtmp1, y1, y2, y3, y4, y5, y6, y7, y8;
double pr_final, pom, pom1, pom2, pom3;
double r, r2, delta, ULt, rms, tmp1, Ut, U_r, UrU_r, Lms, Ems, Ut_rms;
//double q_final, U_phi, Ur;
double thetaO, rin=0., rout, h_rh, elow, ehigh;
double Anorm, Dnorm, g_L, E0, Lx;
double flux_prim_tot, flux_refl_tot, refl_ratio, NpNr;
double zzshift;
double abundance, lensing, gfactor0, ionisation;
double arcosa3;
double time1=0., time1_rg=0., time2=0., time2_rg=0., frequency1=0., 
       frequency2=0., tot_time_sec, deltaT, Af=0., qf=0., fmin, fmax, dt0;
double Tmax, Tmax_final, freq_min, freq_max;
double sumt, sumf, sumt2, sumtf, qtime, Atime;
long   nt, iabun, igam, ixi, it;
int    imin, imax, irh0, ih0, ith0, ir0, iabun0, igam0;
int    nt0, ms, nspec, polar, stokes, rix, photar_sw, npoints;
int    i, ie, je, quit, it0, itn, iband;
//int  itmin;
//double photar1_xspec[ne_xspec]; photar1_short[3], *photar1=NULL;
//double gf_xspec[ne_xspec+1], gf_short[4], *gf=NULL;
double *photar1=NULL, *gf=NULL;

// the following are needed for cut-off power-law taken from XSPEC
// //double ear1[ne + 1], param1[2];
//double *ear1=NULL, param1[2];

// these are needed to work with a fits file...
fitsfile *fptr;
char     tables_file[255];
int      hdutype = 2;
int      colnum = 1;
long     frow = 1, felem = 1, nelems, nrow;
float    float_nulval = 0.;
int      nelements, nelements1, nelements2;
int      ihorizon, irow, anynul, status = 0;

// Let's initialize parameters for subroutine ide()
sprintf(cparam, "%3ld_%02ld_%04ld",
        lround(100.*(1.+sqrt(1.-param[0]*param[0]))),
        lround(param[1]),lround(10.*param[8]));
// am - black hole angular momentum
ide_param[0] = param[0];
am = ide_param[0];
am2 = am * am;
pom1 = pow(1. + am, 1. / 3.);
pom2 = pow(1. - am, 1. / 3.);
pom3 = pow(1. - am2, 1. / 3.);
pom = 1. + pom3 * (pom1 + pom2);
pom1 = sqrt(3. * am2 + pom * pom);
if (am >= 0) rms = 3. + pom1 - sqrt((3. - pom) * (3. + pom + 2. * pom1));
else rms = 3. + pom1 + sqrt((3. - pom) * (3. + pom + 2. * pom1));
r_plus = 1. + sqrt(1. - am2);
Ut_rms=(rms*rms-2.*rms+am*sqrt(rms))/rms/sqrt(rms*rms-3.*rms+2.*am*sqrt(rms));
// thetaO - observer inclination
ide_param[1] = param[1];
thetaO = ide_param[1];
// rin - inner edge of non-zero disc emissivity
ide_param[2] = param[2];
// ms - whether to integrate from rin or rms
ide_param[3] = param[3];
// rout - outer edge of non-zero disc emissivity
ide_param[4] = param[4];
rout = param[4];
//set up inner and outer radii:
ms = (int) ide_param[3];
if(ms!=0 && ms!=1 && ms!=2){
  xs_write("kynrefrev: ms must be 0, 1 or 2",5);
  for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
  goto error;
}
if(ms==1){
 if(ide_param[2]<rms)rin=rms-r_plus;
 else rin=ide_param[2]-r_plus;
 rout=ide_param[4]-r_plus;
}else if(ms==2){
 rin=ide_param[2]*rms-r_plus;
 rout=ide_param[4]*rms-r_plus;
}else if(ms==0){
 rin=ide_param[2]-r_plus;
 rout=ide_param[4]-r_plus;
}
if(rin<0.)rin=0.;
if(rout<0.)rout=0.;
// phi  - lower azimuth of non-zero disc emissivity (deg)
ide_param[5] = param[5];
// dphi - (phi+dphi) is upper azimuth of non-zero disc emissivity (deg)
ide_param[6] = param[6];
// nphi - number of grid points in azimuth
ide_param[9] = param[28];
// smooth - whether to smooth the resulting spectrum (0-no, 1-yes)
ide_param[10] = SMOOTH;
// normal - how to normalize the final spectrum
ide_param[11] = -1.;
// ntable - table model (defines fits file with tables)
ide_param[13] = param[25];
// M/M8 - black hole mass in units of 10^8 solar masses
mass = param[7];
mass2 = mass * mass;
// height - height on the axis (measured from the center) at which the primary
//          source is located (GM/c^2)
h = param[8];
if (h >= 0.)
  if (h < r_plus + 0.1){
    h_rh = 0.1; h=r_plus+0.1;
      xs_write("kynrefrev: too low height, we set it to 0.1 above horizon", 5);
  }else h_rh = h - r_plus;
else {
  xs_write("kynrefrev: height has to be positive.", 5);
  for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
  goto error;
}
// PhoIndex - power-law energy index of the lamp emission
gamma0 = param[9];
// L/Ledd - dE/dt primary isotropic flux in Eddington luminosity as seen by the 
//          disc
Np = param[10];
// Np:Nr - ratio of the primary to the reflected normalization
NpNr = param[11];
if( NpNr > 0. ) Np /= NpNr;
// nH0 - density profile normalization in 10^15 cm^(-3)
nH0 = param[12];
// q_n - radial power-law density profile
qn = param[13];
// Fe abundance (in solar abundance)
abundance = param[14];
// zshift - overall Doppler shift
if (param[21] > 0.) {
  ide_param[12] = param[21];
  Dnorm = pow( HUBBLE / MPC / CLIGHT / param[21], 2 );
}else if (param[21] < 0.) {
  ide_param[12] = 0.;
  Dnorm = pow( - HUBBLE / MPC / CLIGHT / param[21], 2 );
}else {
  ide_param[12] = 0.;
  Dnorm = 1. / ( MPC * MPC );
}
// zzshift - multiplication factor for gfac from zshift needed for primary
zzshift=1.0/(1.0+ide_param[12]);
// limb - table model (defines fits file with tables)
limb = (int) param[22];
if ((limb < 0) || (limb > 2)) {
  xs_write("kynrefrev: limb must be >= 0 and <= 2.", 5);
  for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
  goto error;
}
// tab - which reflion table to use 
rix = (int) param[23];
if (rix == 1) {
  E0 = 0.001;
  sprintf(reflion, REFLION1);
  nener = 500;
  elow = 0.01;
}
else{
  E0 = 0.1;
  sprintf(reflion, REFLION2);
  nener = 375;
  elow = 0.1;
}
ehigh = EC;
// sw  - swich for the way how to compute the refl. spectra
sw = (int) param[24];
// edivision - type of division in local energies (0-equidistant, 1-exponential)
ide_param[14] = 1.;
// variability type
ide_param[15]=1.;
//photar_sw
photar_sw = (int) param[37];
// let's set up the deltaT, ntbin and nn params according to the frequency range 
// needed
dt0 = param[29];
//let's keep the T_tail constant, i.e. nt0 * dt0 = TMAX
nt0 = (int) ceil( TMAX / dt0 );
if(param[4] != 1000.){
//here the last term is for inner echo, specially needed for low rout and low h
  nt0 = (int) ceil( 1.2 / dt0 * 
        ( rout*(1+sin(thetaO/180.*PI)) + h*cos(thetaO/180.*PI) + h*h/rout/2. +
          321./(h+1.) ) );
  if(nt0 > nn) nn = (int) exp2( ceil( log2( nt0 ) ) );
}
Tmax = nt0 * dt0;
Tmax_final = nn * dt0;
fmin = 1. / (nn * dt0 * SEC * mass);
fmax = 1. / (2. * dt0 * SEC * mass);
if(photar_sw){
  if(abs(photar_sw) <= 10){
    freq_min = ( param[31] >=0. ? param[31] : param[31] / ( - SEC * mass));
    freq_max = ( param[32] >=0. ? param[32] : param[32] / ( - SEC * mass));
    if(freq_max < freq_min)freq_max=freq_min;
  }else if(abs(photar_sw)>10){
    if( ear_xspec[0] <= 0. ){
      sprintf(errortxt,
        "kynrefrev: lower frequency boundary has to be larger than zero.");
      xs_write(errortxt, 5);
      for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
      goto error;
    }
    freq_min = ear_xspec[0];
    freq_max = ear_xspec[ne_xspec];
  }
  if( freq_min > 0. ){
    if( log10( freq_max / freq_min) > 5.){
      sprintf(errortxt,
        "kynrefrev: frequency cannot span more than 5 orders of magnitude!");
      xs_write(errortxt, 5);
      for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
      goto error;
    }
    if( freq_max > fmax ){
      dt0 = 1. / freq_max / (2. * SEC * mass);
      if(dt0 < 0.1){
        sprintf(errortxt,
          "kynrefrev: frequency is too high, much above %lg Hz, please, set it below %lg Hz.", 
          fmax, 1. / (0.2 * SEC * mass));
        xs_write(errortxt, 5);
        for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
        goto error;
      }
      sprintf(errortxt,
        "kynrefrev: frequency is quite high, above %lg Hz. Do you really want it?", fmax);
      xs_write(errortxt, 5);
      sprintf(errortxt, "           time step changed to %lg Rg/c.", dt0);
      xs_write(errortxt, 5);
      sprintf(errortxt, "           computation will be longer!");
      xs_write(errortxt, 5);
      nt0 = (int) ceil(Tmax / dt0);
      nn = (int) exp2( ceil( log2( Tmax_final / dt0 ) ) );
    }
    if( freq_min < 1. / (nn * dt0 * SEC * mass) ){
      if(freq_min < 0.1 * fmin){
        sprintf(errortxt,
          "kynrefrev: frequency is too low, much below %lg Hz, please, set it above %lg Hz.", 
          fmin, 0.1 * fmin);
        xs_write(errortxt, 5);
        for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
        goto error;
      }
      nn = (int) exp2( ceil( log2( 1. / freq_min / (dt0 * SEC * mass) ) ) );
      sprintf(errortxt,
        "kynrefrev: frequency is quite low, below %lg Hz. Do you really want it?", fmin);
      xs_write(errortxt, 5);
    }
  }
}
// duration of the flare (if negative then 10 times the step in time)
if( DURATION < 0)flare_duration_rg = - DURATION * dt0;
else flare_duration_rg = DURATION;
flare_duration_sec = flare_duration_rg * SEC * mass;
// nt_ratio (nt)
nt_ratio= (int) param[30];
// nt
nt = ((long) nt_ratio) * ((long) nt0);
if(nt==1){
  xs_write("kynrefrev:  nt = nt_final * nt_ratio must be larger than 1", 5);
  for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
  goto error;
}
// dt
dt=dt0/nt_ratio;
ide_param[16]=dt;
// time/frequency
if (photar_sw == 0){
  time1 = param[31];
  if(time1 >= 0.) time1_rg = time1 / SEC / mass;
  else{
    time1_rg = -time1;
    time1 = -time1 * SEC * mass;
  }
  if (time1_rg > (nt-1)*dt) {
    sprintf(errortxt,"kynrefrev: time has to be smaller than %lg sec or %lg Rg.",
            (nt-1)*dt*SEC*mass, (nt-1)*dt);
    xs_write(errortxt, 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
  time2 = param[32];
  if(time2 >= 0.) time2_rg = time2 / SEC / mass;
  else{
    time2_rg = -time2;
    time2 = -time2 * SEC * mass;
  }
  if (time2_rg > (nt-1)*dt) {
    sprintf(errortxt,"kynrefrev: time has to be smaller than %lg sec or %lg Rg.",
            (nt-1)*dt*SEC*mass, (nt-1)*dt);
    xs_write(errortxt, 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
}else if(abs(photar_sw) <= 10){
  fmin = 1. / (nn * dt0 * SEC * mass);
  fmax = 1. / (2. * dt0 * SEC * mass);
  frequency1 = param[31];
  if(frequency1 < 0.) frequency1 /= ( - SEC * mass);
  if (frequency1 != 0. && frequency1 < fmin ) {
    sprintf(errortxt,"kynrefrev: frequency has to be larger than %lg Hz or %lg / Rg.",
            fmin, 1. / (nn * dt0));
    xs_write(errortxt, 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
  if (frequency1 != 0. && frequency1 > fmax ) {
    sprintf(errortxt,"kynrefrev: frequency has to be smaller than %lg Hz or %lg / Rg.",
            fmax, 1. / (2. * dt0));
    xs_write(errortxt, 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
  frequency2 = param[32];
  if(frequency2 < 0.) frequency2 /= ( - SEC * mass);
  if (frequency1 != 0. && frequency2 > fmax ) {
    sprintf(errortxt,"kynrefrev: frequency has to be smaller than %lg Hz or %lg / Rg.",
            fmax, 1. / (2. * dt0));
    xs_write(errortxt, 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
}else{
//definition of the energy band of interest for frequency dependent Fourier transform
  en1 = param[31];
  en2 = param[32];
  if ( en1 < 0. || en1 >= en2 ) {
    sprintf(errortxt,"kynrefrev: wrong definition of energy band of interest (0 <= E1 < E2).");
    xs_write(errortxt, 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
//check the frequency range
  fmin = 1. / (nn * dt0 * SEC * mass);
  fmax = 1. / (2. * dt0 * SEC * mass);
  if( ear_xspec[0] < fmin || ear_xspec[0] > fmax || 
      ear_xspec[ne_xspec] < fmin || ear_xspec[ne_xspec] > fmax ){
    sprintf(errortxt,"kynrefrev: frequency has to be in the range of %lg to %lg Hz.",
            fmin, fmax);
    xs_write(errortxt, 5);
  }
}
//definition of the reference energy band for both frequency dependent lags and 
//amplitudes as well as energy dependent lags and amplitudes
en3 = param[33];
en4 = param[34];
if ( ( ( abs(photar_sw) == 5 || abs(photar_sw) == 6 ||
         abs(photar_sw) == 7 || abs(photar_sw) == 8 ) && 
       ( en3 > 0. && en3 >= en4 ) ) ||
     ( ( abs(photar_sw) == 15 || abs(photar_sw) == 16 || 
         abs(photar_sw) == 17 || abs(photar_sw) == 18 ) && ( fabs(en3) >= fabs(en4) ) ) ) {
  sprintf(errortxt,"kynrefrev: wrong definition of the reference energy band.");
  xs_write(errortxt, 5);
  for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
  goto error;
}
exclude_energy=0;
#ifdef OUTSIDE_XSPEC
if( abs(photar_sw) > 10 ){
  nbands=1;
  ener_low[0]=en1;
  ener_high[0]=en2;
}else{
  nbands = NBANDS;
  if(ener_low[nbands-1] == ear_xspec[0] &&
     ener_high[nbands-1] == ear_xspec[ne_xspec])exclude_energy=1;
}
if( abs(photar_sw) == 15 || abs(photar_sw) == 16 || 
    abs(photar_sw) == 17 || abs(photar_sw) == 18 ){
  nbands++;
  ener_low[nbands-1] = fabs(en3);
  ener_high[nbands-1] = fabs(en4);
  if(en3<0.)exclude_energy=1.;
}
#else
if( abs(photar_sw) > 10 ){
  nbands=1;
  ener_low[0]=en1;
  ener_high[0]=en2;
}else nbands = 0;
if( abs(photar_sw) == 5 || abs(photar_sw) == 6 ||
    abs(photar_sw) == 7 || abs(photar_sw) == 8 ){
  if(en3 < 0){
    nbands++;
    ener_low[nbands-1] = ear_xspec[0];
    ener_high[nbands-1] = ear_xspec[ne_xspec];
    exclude_energy = 1;
  }else if( en3 > 0 && en3 < en4 ){
    nbands++;
    ener_low[nbands-1] = en3;
    ener_high[nbands-1] = en4;
    exclude_energy = 0;
  }
}
if( abs(photar_sw) == 15 || abs(photar_sw) == 16 ||
    abs(photar_sw) == 17 || abs(photar_sw) == 18 ){
  if( fabs(en3) < fabs(en4) ){
    nbands++;
    ener_low[nbands-1] = fabs(en3);
    ener_high[nbands-1] = fabs(en4);
    if( en3 < 0 )exclude_energy = 1;
    else exclude_energy = 0;
  }
}
#endif
//let's define the energy for ide integration for energy dependent lags as well
//as for frequency dependent lags of energy bands
if(abs(photar_sw) <= 10){
  ear=ear_xspec;
  ne=ne_xspec;
}else{
  ear=ear_short;
//let's define ne and new ear energy array
  if( (en1==fabs(en3) && en2==fabs(en4)) || abs(photar_sw) < 15){
    ne=1;
    ear[0]=en1;
    ear[1]=en2;
  }else if(en1==fabs(en3)){
    ne=2;
    ear[0]=en1;
    if(en2<fabs(en4)) {ear[1]=en2;ear[2]=fabs(en4);}
    else {ear[1]=fabs(en4);ear[2]=en2;}
  }else if(en2==fabs(en4)){
    ne=2;
    ear[2]=en2;
    if(en1<fabs(en3)) {ear[0]=en1;ear[1]=fabs(en3);}
    else {ear[0]=fabs(en3);ear[1]=en1;}
  }else if(en2==fabs(en3)){
    ne=2;
    ear[0]=en1;
    ear[1]=en2;
    ear[2]=fabs(en4);
  }else if(en1==fabs(en4)){
    ne=2;
    ear[0]=fabs(en3);
    ear[1]=en1;
    ear[2]=en2;
  }else if(en1<fabs(en3)){
    ne=3;
    ear[0]=en1;
    if(en2<fabs(en3)){ear[1]=en2;ear[2]=fabs(en3);ear[3]=fabs(en4);}
    else if(en2<fabs(en4)){ear[1]=fabs(en3);ear[2]=en2;ear[3]=fabs(en4);}
    else {ear[1]=fabs(en3);ear[2]=fabs(en4);ear[3]=en2;}
  }else if(fabs(en3) < en1){
    ne=3;
    ear[0]=fabs(en3);
    if(fabs(en4)<en1){ear[1]=fabs(en4);ear[2]=en1;ear[3]=en2;}
    else if(fabs(en4)<en2){ear[1]=en1;ear[2]=fabs(en4);ear[3]=en2;}
    else {ear[1]=en1;ear[2]=en2;ear[3]=fabs(en4);}
  }
}
//Create arrays of different dimensions according to photar_sw, i.e. different
//dimensions if ear_xspec contains energy bins or frequency bins!
ntmax=( nt > nn ? nt : nn);
if ( (time = (double *) malloc( ntmax * sizeof(double)) )  == NULL ||
     (far  = (double *) malloc( ne * ntmax * sizeof(double)) )       == NULL ||
     (far_final = (double *) malloc( ne * ntmax * sizeof(double)))   == NULL ||
     (flux = (double *) malloc( ntmax * sizeof(double)))   == NULL ||
     (flux_bands = (double *) malloc( nbands * nn * sizeof(double))) == NULL ||
     (spectrum = (double *) malloc( ne * sizeof(double)))  == NULL ||
     (far_prim = (double *) malloc( ne * sizeof(double)))  == NULL ||
     (spectrum_prim = (double *) malloc( ne * sizeof(double))) == NULL ||
     (photar1 = (double *) malloc( ne * sizeof(double)))   == NULL ||
     (gf = (double *) malloc( (ne+1) * sizeof(double)))    == NULL ) {
  xs_write("kynrefrev: Failed to allocate memory for tmp arrays.", 5);
  for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
  goto error;
}
//if ((ear1 = (double *) malloc( (ne+1) * sizeof(double))) == NULL) {
//  xs_write("kynrefrev: Failed to allocate memory for tmp arrays.", 5);
//  for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
//  goto error;
//}
if( photar_sw == 5 || photar_sw == 6 || photar_sw == 7 || photar_sw == 8 ){
  lag_shift = param[35];
  ampl_ampl = param[36];
}
if(photar_sw == 16 || photar_sw == 18){
  if(param[35]==-1.){
    Af = 350.e-4*log10((fabs(en3)+fabs(en4))/(en1+en2));
    qf = 1.;
  }else{
    Af = param[35];
    qf = param[36];
  }
}
// polar - whether we need value of change in polarization angle (0-no,1-yes)
//stokes = (int) param[39];
stokes = 0;
if ((stokes < 0) || (stokes > 6)) {
  xs_write("kynrefrev: stokes has to be 0-6.", 5);
  for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
  goto error;
}
polar = 0;
if (stokes > 0) polar = 1;
ide_param[17] = polar;
//delay_r
ide_param[18]=-1.;
// delay_phi is defined later on
//ide_param[19]=del_a
// number of threads for multithread computations
ide_param[20] = param[38];
// thermalisation - thermalisation fraction
thermalisation = param[15];
// arate - accretion rate (note that LEDD is in erg/s not W and CLIGHT in km/s
//          not m/s and 1e-13 = erg/s/W / (km/m)^2 = 1e-7 / (1e3)^2)
if(param[16] < 0.) arate = fabs(param[16]);
else arate = param[16]*LEDD*1e-13*mass/MSOLAR*YEAR/CLIGHT/CLIGHT/(1-Ut_rms);
// f_col - colour correction factor
f_col = param[17];
if(thermalisation != 0.){
  Tnorm = f_col * K_KEVK * sqrt( C_CMS ) *
          pow( 3. * MSOLAR / (8. * PI * RG2 * YEAR * SIGMA), 0.25);
  x0 = sqrt(rms);
  arcosa3 = acos(am) / 3.;
  x1 = 2 * cos(arcosa3 - PI / 3.);
  x2 = 2 * cos(arcosa3 + PI / 3.);
  x3 = -2 * cos(arcosa3);
}
// alpha - position of the cloud in alpha impact parameter (in GM/c^2)
ide_param[21] = param[18];
// beta - position of the cloud in beta impact parameter (in GM/c^2)
ide_param[22] = param[19];
// rcloud - radius of the cloud (in GM/c^2)
ide_param[23] = param[20];
//whether the flux defined in emissivity subroutine is local one (0) or the 
//observed one (1)
ide_param[24] = 0.;

// check if normalization parameter is equal to 1.
if ((param[21] != 0.) && (param[39] != 1.)) {
  xs_write("kynrefrev: the normalisation parameter par34 should be frozen to unity.", 5);
}
// Let's clear the flux array
for (ie = 0; ie < ne; ie++)
  for (it = 0; it < nt; it++)
    far[ie + ne * it] = 0;

/******************************************************************************/
// Let's read the lamp post tables
if(first_h) {
// The status parameter must always be initialized.
  status = 0;
// Open the FITS file for readonly access
// - if set try KYDIR directory, otherwise look in the working directory
//   or in the xspec directory where tables are usually stored...
  sprintf(kydir, "%s", FGMSTR(pname));
  if (strlen(kydir) == 0) sprintf(tables_file, "./%s", LAMP);
  else if (kydir[strlen(kydir) - 1] == '/') sprintf(tables_file, "%s%s",
                                                    kydir, LAMP);
  else sprintf(tables_file, "%s/%s", kydir, LAMP);
// Let's read the 'KBHlamp_q' fits file
// The status parameter must always be initialized.
  status = 0;
  ffopen(&fptr, tables_file, READONLY, &status);
  if (status) {
    sprintf(tables_file, "%s%s", FGMODF(), LAMP);
    status = 0;
    ffopen(&fptr, tables_file, READONLY, &status);
  }
  if (status) {
    if (status) ffrprt(stderr, status);
    ffclos(fptr, &status);
    xs_write("\nkynrefionx: set the KYDIR to the directory with the KY tables.", 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
// Let's read tables (binary tables => hdutype=2)
// Move to the extension 'r_horizon' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &nrh, &status);
//******************************************************************************
//  fprintf(stdout,"nrh = %ld\n",nrh);
//******************************************************************************   
// Allocate memory for r_horizon...
  if ((r_horizon = (float *) malloc(nrh * sizeof(float))) == NULL) {
    xs_write("kynrefrev: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
// Read the data in the 'r_horizon' table
  nelems = nrh;
// FFGCV reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, r_horizon,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nrh; i++)fprintf(stdout,"%f\n",r_horizon[i]);
//******************************************************************************   
// Move to the extension 'height' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &nh, &status);
//******************************************************************************
//  fprintf(stdout,"nh = %ld\n",nh);
//******************************************************************************   
// Allocate memory for height...
  if ((height = (float *) malloc(nh * sizeof(float))) == NULL) {
    xs_write("kynrefrev: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
// Read the data in the 'height' table
  nelems = nh;
// FFGCV reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, height,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nh; i++)fprintf(stdout,"%f\n",height[i]);
//******************************************************************************   
// Move to the extension 'inclination' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &nincl, &status);
//******************************************************************************
//  fprintf(stdout,"nincl = %ld\n",nincl);
//******************************************************************************   
// Allocate memory for inclination...
  if ((incl = (float *) malloc(nincl * sizeof(float))) == NULL) {
    xs_write("kynrefrev: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
// Read the data in the 'inclination' table
  nelems = nincl;
// FFGCV reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, incl,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nincl; i++)fprintf(stdout,"%f\n",incl[i]);
//******************************************************************************   
// Move to the extension 'r_rh' and read its values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &nradius, &status);
//******************************************************************************
//  fprintf(stdout,"nradius = %ld\n",nradius);
//******************************************************************************   
// Allocate memory for r_rh...
  if ((radius = (float *) malloc(nradius * sizeof(float))) == NULL) {
    xs_write("kynrefrev: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
// Read the data in the 'r_rh' table
  nelems = nradius;
// FFGCV reads the VALUES from the first column.
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, radius,
        &anynul, &status);
//******************************************************************************
//  for ( i=0; i<nradius; i++)fprintf(stdout,"%f\n",radius[i]);
//******************************************************************************   
// Let's read the tables for dWadWo, q, p^r and dWadSd
// allocate memory for the arrays
  if ((dWadWo  = (float *) malloc(nincl * nh * nrh * sizeof(float))) == NULL ||
      (delay_a = (float *) malloc(nincl * nh * nrh * sizeof(float))) == NULL ||
      (q     = (float *) malloc(nradius * nh * nrh * sizeof(float))) == NULL ||
      (pr    = (float *) malloc(nradius * nh * nrh * sizeof(float))) == NULL ||
      (dWadSd= (float *) malloc(nradius * nh * nrh * sizeof(float))) == NULL ||
      (delay = (float *) malloc(nradius * nh * nrh * sizeof(float))) == NULL ){
    xs_write("kynrefrev: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
// read the tables
  for (ihorizon = 0; ihorizon < nrh; ihorizon++) {
    ffmrhd(fptr, 1, &hdutype, &status);
/*  to read the file only once we have to read in blocks (all columns
    from the extension are put to buffer together)
    let's find out how many rows are going to be read into the buffer */
    ffgrsz(fptr, &nrow, &status);
    nelements1 = nrow * nincl;
    nelements2 = nrow * nradius;
    for (irow = 0; irow < nh; irow += nrow) {
//    the last block to read may be smaller:
      if ((nradius - irow) < nrow) {
        nelements1 = (nradius - irow) * nincl;
        nelements2 = (nradius - irow) * nradius;
      }
      ffgcv(fptr, TFLOAT, 1, irow + 1, 1, nelements1, &float_nulval, 
            &dWadWo[irow * nincl + nh * nincl * ihorizon],
            &anynul, &status);
      ffgcv(fptr, TFLOAT, 2, irow + 1, 1, nelements1, &float_nulval, 
            &delay_a[irow * nincl + nh * nincl * ihorizon],
            &anynul, &status);
      ffgcv(fptr, TFLOAT, 3, irow + 1, 1, nelements2, &float_nulval, 
            &q[irow * nradius + nh * nradius * ihorizon],
            &anynul, &status);
      ffgcv(fptr, TFLOAT, 4, irow + 1, 1, nelements2, &float_nulval, 
            &pr[irow * nradius + nh * nradius * ihorizon],
            &anynul, &status);
      ffgcv(fptr, TFLOAT, 5, irow + 1, 1, nelements2, &float_nulval, 
            &dWadSd[irow * nradius + nh * nradius * ihorizon],
            &anynul, &status);
      ffgcv(fptr, TFLOAT, 6, irow + 1, 1, nelements2, &float_nulval, 
            &delay[irow * nradius + nh * nradius * ihorizon],
            &anynul, &status);
    }
  }
// The FITS file must always be closed before exiting the program.
  ffclos(fptr, &status);
/*******************************************************************************
  irh=20;
  ih=99;
  for ( i=0; i<nincl; i++ ) fprintf(stdout,"%d\t%f\t%f\n",i,incl[i],
    dWadWo[i+nincl*ih+nincl*nh*irh]);
  for ( i=0; i<nradius; i++) fprintf(stdout,"%d\t%f\t%f\t%f\t%f\n",i,radius[i],
    q[i+nradius*ih+nradius*nh*irh],pr[i+nradius*ih+nradius*nh*irh],
    dWadSd[i+nradius*ih+nradius*nh*irh]);
*******************************************************************************/
// Firstly we have to free allocated memory for the arrays gfac,
// cosin, phiph, transf_d
  if((gfac    = (double *) malloc(nradius * sizeof(double))) == NULL ||
//     (cosin   = (double *) malloc(nradius * sizeof(double))) == NULL ||
//     (phiph   = (double *) malloc(nradius * sizeof(double))) == NULL ||
     (transf_d= (double *) malloc(nradius * sizeof(double))) == NULL ||
     (del     = (double *) malloc(nradius * sizeof(double))) == NULL ) {
    xs_write("kynrefrev: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
  first_h = 0;
}
/******************************************************************************/
if (h_rh > height[nh-1]) {
  sprintf(errortxt, "kynrefrev: the height must be lower than or equal to %f.",
          height[nh - 1] + r_plus);
  xs_write(errortxt, 5);
  for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
  goto error;
}
if (h_rh < height[0]) {
  sprintf(errortxt, "kynrefrev: the height is too low, we set it to %f.", 
          height[0] + r_plus);
  xs_write(errortxt, 5);
  h_rh = height[0];
  h = h_rh + r_plus;
}
// nrad - number of grid points in radius
if(param[26] > 0) ide_param[7] = param[26];
else ide_param[7] = - param[26] * pow(h,-0.66);
// division - type of division in r integration (0-equidistant, 1-exponential)
if(param[27]==-1.)ide_param[8]=2.*h;
else ide_param[8]=param[27];
// Let's interpolate the tables to desired spin and height
if ((am != am_old) || (h_rh != h_rh_old) || (thetaO != thetaO_old)) {
// given am->r_plus, find the corresponding index in r_horizon[]:
  imin = 0;
  imax = nrh;
  irh0 = nrh / 2;
  while ((imax - imin) > 1) {
    if (r_plus >= r_horizon[irh0 - 1]) imin = irh0;
    else imax = irh0;
    irh0 = (imin + imax) / 2;
  }
  if (irh0 == 0) irh0 = 1;
//if ((imax == nrh) && (r_plus > r_horizon[nrh - 1])) irh0 = nrh;
  ttmp = (r_plus - r_horizon[irh0 - 1]) / 
         (r_horizon[irh0] - r_horizon[irh0 - 1]);
  ttmp1 = 1. - ttmp;
// given h, find the corresponding index in height[]:
  imin = 0;
  imax = nh;
  ih0 = nh / 2;
  while ((imax - imin) > 1) {
    if (h_rh >= height[ih0 - 1]) imin = ih0;
    else imax = ih0;
    ih0 = (imin + imax) / 2;
  }
  if (ih0 == 0) ih0 = 1;
//if ((imax == nh) && (h_rh > height[nh - 1])) ih0 = nh;
  utmp = (h_rh - height[ih0 - 1]) / (height[ih0] - height[ih0 - 1]);
  utmp1 = 1. - utmp;
// given thetaO, find the corresponding index in incl[]:
  imin = 0;
  imax = nincl;
  ith0 = nincl / 2;
  while ((imax - imin) > 1) {
    if (thetaO >= incl[ith0 - 1]) imin = ith0;
    else imax = ith0;
    ith0 = (imin + imax) / 2;
  }
  if (ith0 == 0) ith0 = 1;
//if ((imax == nincl) && (thetaO > incl[nincl - 1])) ith0 = nincl;
  vtmp = (thetaO - incl[ith0 - 1]) / (incl[ith0] - incl[ith0 - 1]);
  vtmp1 = 1. - vtmp;
// transfer function from the axis to the observer
  y1 = dWadWo[ith0 - 1 + nincl * (ih0 - 1) + nincl * nh * (irh0 - 1)];
  y2 = dWadWo[ith0 - 1 + nincl * (ih0 - 1) + nincl * nh * irh0];
  y3 = dWadWo[ith0 - 1 + nincl * ih0 + nincl * nh * irh0];
  y4 = dWadWo[ith0 - 1 + nincl * ih0 + nincl * nh * (irh0 - 1)];
  y5 = dWadWo[ith0 + nincl * (ih0 - 1) + nincl * nh * (irh0 - 1)];
  y6 = dWadWo[ith0 + nincl * (ih0 - 1) + nincl * nh * irh0];
  y7 = dWadWo[ith0 + nincl * ih0 + nincl * nh * irh0];
  y8 = dWadWo[ith0 + nincl * ih0 + nincl * nh * (irh0 - 1)];
  transf_o = (vtmp1 * (utmp1 * (ttmp1 * y1 + ttmp * y2) + 
                       utmp *  (ttmp * y3 + ttmp1 * y4)) + 
              vtmp * (utmp1 * (ttmp1 * y5 + ttmp * y6) + 
                      utmp *  (ttmp * y7 + ttmp1 * y8)));
// delay from the axis to the observer
  y1 = delay_a[ith0 - 1 + nincl * (ih0 - 1) + nincl * nh * (irh0 - 1)];
  y2 = delay_a[ith0 - 1 + nincl * (ih0 - 1) + nincl * nh * irh0];
  y3 = delay_a[ith0 - 1 + nincl * ih0 + nincl * nh * irh0];
  y4 = delay_a[ith0 - 1 + nincl * ih0 + nincl * nh * (irh0 - 1)];
  y5 = delay_a[ith0 + nincl * (ih0 - 1) + nincl * nh * (irh0 - 1)];
  y6 = delay_a[ith0 + nincl * (ih0 - 1) + nincl * nh * irh0];
  y7 = delay_a[ith0 + nincl * ih0 + nincl * nh * irh0];
  y8 = delay_a[ith0 + nincl * ih0 + nincl * nh * (irh0 - 1)];
  del_a = (vtmp1 * (utmp1 * (ttmp1 * y1 + ttmp * y2) + 
                    utmp *  (ttmp * y3 + ttmp1 * y4)) + 
           vtmp * (utmp1 * (ttmp1 * y5 + ttmp * y6) + 
                   utmp *  (ttmp * y7 + ttmp1 * y8)));
  if ((am != am_old) || (h_rh != h_rh_old)) {
    for (i = 0; i < nradius; i++) {
/*
// q from the axis to the disc
      y1 = q[i + nradius * (ih0 - 1) + nradius * nh * (irh0 - 1)];
      y2 = q[i + nradius * (ih0 - 1) + nradius * nh * irh0];
      y3 = q[i + nradius * ih0 + nradius * nh * irh0];
      y4 = q[i + nradius * ih0 + nradius * nh * (irh0 - 1)];
      q_final = utmp1 * (ttmp1 * y1 + ttmp * y2) + 
                utmp * (ttmp * y3 + ttmp1 * y4);
*/
// pr at the disc
      y1 = pr[i + nradius * (ih0 - 1) + nradius * nh * (irh0 - 1)];
      y2 = pr[i + nradius * (ih0 - 1) + nradius * nh * irh0];
      y3 = pr[i + nradius * ih0 + nradius * nh * irh0];
      y4 = pr[i + nradius * ih0 + nradius * nh * (irh0 - 1)];
      pr_final = utmp1 * (ttmp1 * y1 + ttmp * y2) + 
                 utmp * (ttmp * y3 + ttmp1 * y4);
// temporary variables
      r = r_plus + radius[i];
      r2 = r * r;
      delta = r2 - 2. * r + am2;
      ULt = sqrt((h * h + am2) / (h * h - 2. * h + am2));
      if (r >= rms) {
        tmp1 = sqrt(r2 - 3. * r + 2. * am * sqrt(r));
        Ut = (r2 + am * sqrt(r)) / r / tmp1;
//      U_phi = (r2 + am2 - 2. * am * sqrt(r)) / sqrt(r) / tmp1;
        U_r = 0.;
//      Ur = 0.;
        UrU_r = 0.;
      }
      else {
        tmp1 = sqrt(rms * (rms - 3.) + 2. * am * sqrt(rms));
        Lms = (rms * rms + am2 - 2. * am * sqrt(rms)) / sqrt(rms) / tmp1;
        Ems = (rms * (rms - 2.) + am * sqrt(rms)) / rms / tmp1;
        Ut = (Ems * (r * (r2 + am2) + 2. * am2) - 2. * am * Lms) / r / delta;
//      U_phi = Lms;
        UrU_r = -1. + ((r2 + am2 + 2. * am2 / r) * Ems * Ems - 4. * am / 
                r * Ems * Lms - (1. - 2. / r) * Lms * Lms) / delta;
        if (UrU_r < 0.) UrU_r = 0.;
        U_r = -sqrt(UrU_r / delta) * r;
//      Ur = -sqrt(delta * UrU_r) / r;
      }
      tmp1 = Ut - pr_final * U_r;
// gfactor  from the axis to the disc
      gfac[i] = tmp1 / ULt;
// cosin at the disc
//    cosin[i] = q_final / r / tmp1;
// phip_i at the disc
//    phiph[i] = atan2(-U_phi, r * (pr_final - Ur *tmp1));
// dWadSd from the axis to the disc
      y1 = dWadSd[i + nradius * (ih0 - 1) + nradius * nh * (irh0 - 1)];
      y2 = dWadSd[i + nradius * (ih0 - 1) + nradius * nh * irh0];
      y3 = dWadSd[i + nradius * ih0 + nradius * nh * irh0];
      y4 = dWadSd[i + nradius * ih0 + nradius * nh * (irh0 - 1)];
      transf_d[i] = utmp1 * (ttmp1 * y1 + ttmp * y2) + 
                    utmp * (ttmp * y3 + ttmp1 * y4);
// delay from the axis to the disc
      y1 = delay[i + nradius * (ih0 - 1) + nradius * nh * (irh0 - 1)];
      y2 = delay[i + nradius * (ih0 - 1) + nradius * nh * irh0];
      y3 = delay[i + nradius * ih0 + nradius * nh * irh0];
      y4 = delay[i + nradius * ih0 + nradius * nh * (irh0 - 1)];
      del[i] = utmp1 * (ttmp1 * y1 + ttmp * y2) + 
               utmp * (ttmp * y3 + ttmp1 * y4);
    }
  }
//******************************************************************************
//    fprintf(stdout,"%f %f\n", thetaO, transf_o);
//    for(i = 0; i < nradius; i++) 
//      fprintf(stdout,"%d %f %f %f %f %f\n", i, radius[i], gfac[i], cosin[i], 
//              phiph[i], transf_d[i]);
//******************************************************************************
}
/******************************************************************************/
// Let's read the reflionx tables
if ((strcmp(kydir, FGMSTR(pname)) != 0) || (rix != rix_old) || first_rix) {
  sprintf(text, "kynrefrev: initializing %s tables...", reflion);
  xs_write(text, 5);
  xs_write("Ross & Fabian (2005), MNRAS, 358, 211",5);
// The status parameter must always be initialized.
  status = 0;
// Open the FITS file for readonly access
// - if set try KYDIR directory, otherwise look in the working directory
//   or in the xspec directory where tables are usually stored...
  sprintf(kydir, "%s", FGMSTR(pname));
  if (strlen(kydir) == 0) sprintf(tables_file, "./%s", reflion);
  else if (kydir[strlen(kydir) - 1] == '/') sprintf(tables_file, "%s%s",
                                                    kydir, reflion);
  else sprintf(tables_file, "%s/%s", kydir, reflion);
// Let's read the reflion(x) tables
// The status parameter must always be initialized.
  status = 0;
  ffopen(&fptr, tables_file, READONLY, &status);
  if (status) {
    sprintf(tables_file, "%s%s", FGMODF(), reflion);
    status = 0;
    ffopen(&fptr, tables_file, READONLY, &status);
  }
  if (status) {
    if (status) ffrprt(stderr, status);
    ffclos(fptr, &status);
    xs_write("\nkynrefionx: set the KYDIR to the directory with the KY tables.", 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
  if (((strcmp(kydir, FGMSTR(pname)) != 0) || (rix != rix_old)) && !first_rix) {
// Free memory from tmp arrays...
    if(abun != NULL); free((void *) abun); abun = NULL;
    if(gam != NULL); free((void *) gam); gam = NULL;
    if(xi != NULL); free((void *) xi); xi = NULL;
    if(logxi != NULL); free((void *) logxi); logxi = NULL;
    if(energy0 != NULL); free((void *) energy0); energy0 = NULL;
    if(emission != NULL); free((void *) emission); emission = NULL;
    if(flux0 != NULL); free((void *) flux0); flux0 = NULL;
    if(flux1 != NULL); free((void *) flux1); flux1 = NULL;
    if(energy1 != NULL); free((void *) energy1); energy1 = NULL;
    if(energy2 != NULL); free((void *) energy2); energy2 = NULL;
    if(flx != NULL); free((void *) flx); flx = NULL;
  }
// Let's read tables (binary tables => hdutype=2)
// Move to the first extension ('PARAMETERS')
  ffmrhd(fptr, 1, &hdutype, &status);
// Read the values of parameters --> abundance, photon index and ionisation
  nelems = 1;
  colnum = 9;
  frow = 1;
  ffgcv(fptr, TLONG, colnum, frow, felem, nelems, &float_nulval, &nabun,
        &anynul, &status);
  frow = 2;
  ffgcv(fptr, TLONG, colnum, frow, felem, nelems, &float_nulval, &ngam,
        &anynul, &status);
  frow = 3;
  ffgcv(fptr, TLONG, colnum, frow, felem, nelems, &float_nulval, &nxi,
        &anynul, &status);
    nspec = nabun * ngam * nxi;
//******************************************************************************
//  fprintf(stdout,"%d %d %d\n", nabun, ngam, nxi);
//******************************************************************************
// Allocate memory for abun, gam and xi...
  if( (abun  = (float *) malloc(nabun * sizeof(float))) == NULL ||
      (gam   = (float *) malloc( ngam * sizeof(float))) == NULL ||
      (xi    = (float *) malloc(  nxi * sizeof(float))) == NULL ||
      (logxi = (float *) malloc(  nxi * sizeof(float))) == NULL ) {
    xs_write("kynrefrev: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
  colnum = 10;
  frow = 1;
  nelems = nabun;
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, abun,
        &anynul, &status);
  frow = 2;
  nelems = ngam;
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, gam,
        &anynul, &status);
  frow = 3;
  nelems = nxi;
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, xi,
        &anynul, &status);
  for (i = 0; i < nxi; i++) logxi[i] = log10(xi[i]);
//******************************************************************************
//  for(i = 0; i < 20; i++)
//    fprintf(stdout,"%4d %12.1f %12.1f %12.1f %12.1f\n", i, abun[i], gam[i], 
//            xi[i], logxi[i]);
//******************************************************************************
// Move to the second extension 'ENERGIES' and read energy values
  ffmrhd(fptr, 1, &hdutype, &status);
  ffgnrw(fptr, &ne_loc, &status);
  ne_loc++;
//******************************************************************************
//  fprintf(stdout,"%d\n", ne_loc);
//******************************************************************************
// Allocate memory for energy...
  if ((energy0 = (float *) malloc(ne_loc * sizeof(float))) == NULL) {
    xs_write("kynrefrev: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
  nelems = ne_loc - 1;
  colnum = 1;
  frow = 1;
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, energy0,
        &anynul, &status);
  nelems = 1;
  colnum = 2;
  frow = ne_loc - 1;
  ffgcv(fptr, TFLOAT, colnum, frow, felem, nelems, &float_nulval, &energy0[ne_loc - 1],
        &anynul, &status);
//******************************************************************************
//  for(i = 0; i < ne_loc; i++) fprintf(stdout,"%d %f\n", i, energy0[i]);
//******************************************************************************
// Allocate memory for emission...
  if ((emission = (float *) malloc(nspec * (ne_loc - 1) * sizeof(float))) == NULL) {
    xs_write("kynrefrev: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
// Let's read the tables for emission
  ffmrhd(fptr, 1, &hdutype, &status);
// to read the file only once we have to read in blocks (all columns
// from the extension are put to buffer together)
// let's find out how many rows are going to be read into the buffer
  ffgrsz(fptr, &nrow, &status);
  nelements = nrow * (ne_loc - 1);
  for (irow = 0; irow < nspec; irow += nrow) {
// the last block to read may be smaller:
    if ((nspec - irow) < nrow) nelements = (nspec - irow) * (ne_loc - 1);
    iabun = irow / (ngam * nxi) + 1;
    igam = (irow - (iabun - 1) * ngam * nxi) / nxi + 1;
    ixi = irow + 1 - (iabun - 1) * ngam * nxi - (igam - 1) * nxi;
    ffgcv(fptr, TFLOAT, 2, irow + 1, 1, nelements, &float_nulval, 
          &emission[(ne_loc - 1) * (ixi - 1) + (ne_loc - 1) * nxi * (igam - 1) + 
                    (ne_loc - 1) * nxi * ngam * (iabun - 1)], &anynul, &status);
  }
// The FITS file must always be closed before exiting the program.
  ffclos(fptr, &status);
/*******************************************************************************
  iabun = 7;
  igam = 5;
  ixi = 3;
  fprintf(stdout, "%f %f %f %f\n", abun[iabun - 1], gam[igam - 1], 
                                   pow(10, xi[ixi - 1]), xi[ixi - 1]);
  for(ie = 0; ie < ne_loc - 1; ie++)
    fprintf(stdout, "%d %f %f\n", ie, energy0[ie], emission[ie + 
      (ne_loc - 1) * (ixi - 1) + (ne_loc - 1) * nxi * (igam - 1) + 
      (ne_loc - 1 ) * nxi * ngam * (iabun - 1)]);
*******************************************************************************/
// Allocate memory for local flux...
  if( (flux0   = (double *) malloc(ne_loc * nxi * sizeof(double))) == NULL ||
      (flux1   = (double *) malloc( nener * nxi * sizeof(double))) == NULL ||
      (energy1 = (double *) malloc( nener * sizeof(double))) == NULL ||
      (energy2 = (double *) malloc( (nener + 1) * sizeof(double))) == NULL ||
      (flx     = (double *) malloc( nener * sizeof(double))) == NULL ) {
    xs_write("kynrefrev: Failed to allocate memory for tmp arrays.", 5);
    for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
    goto error;
  }
  first_rix = 0;
  xs_write("------------------------------------------------", 5);
}
// end of reading reflionx fits file............................................
/******************************************************************************/
// Let's interpolate the emission to get rid of the non-radial dependences, i.e.
// interpolate in abundance and gamma0
if ((rix != rix_old) || (strcmp(kydir, FGMSTR(pname)) != 0) || ((abun_old == -1.) 
    || (abundance != abun_old)) || ((gam_old == -1) || (gamma0 != gam_old))) {
// given gamma0, find the corresponding index in gam():
  imin = 0;
  imax = ngam;
  igam0 = ngam / 2;
  while ((imax - imin) > 1) {
    if (gamma0 >= gam[igam0 - 1]) imin = igam0;
    else imax = igam0;
    igam0 = (imin + imax) / 2;
  }
  if (igam0 == 0) igam0 = 1;
  ttmp = (gamma0 - gam[igam0 - 1]) / (gam[igam0] - gam[igam0 - 1]);
  if (ttmp < 0.) ttmp = 0.;
  if (ttmp > 1.) ttmp = 1.;
  ttmp1 = 1. - ttmp;
// given abundance, find the corresponding index in abun[]:
  imin = 0;
  imax = nabun;
  iabun0 = nabun / 2;
  while ((imax - imin) > 1) {
    if (abundance >= abun[iabun0 - 1]) imin = iabun0;
    else imax = iabun0;
    iabun0 = (imin + imax) / 2;
  }
  if (iabun0 == 0) iabun0 = 1;
  utmp = (abundance - abun[iabun0 - 1]) / (abun[iabun0] - abun[iabun0 - 1]);
  if (utmp < 0.) utmp = 0.;
  if (utmp > 1.) utmp = 1.;
  utmp1 = 1. - utmp;
// In the following the last energy index is not used 
// (we have one less flux indices then in energy)
  for (ixi = 0; ixi < nxi; ixi++)
    for (ie = 0; ie < ne_loc - 1; ie++) {
      y1 = emission[ie + (ne_loc - 1) * ixi + (ne_loc - 1) * nxi * (igam0 - 1) + 
                    (ne_loc - 1) * nxi * ngam * (iabun0 - 1)];
      y2 = emission[ie + (ne_loc - 1) * ixi + (ne_loc - 1) * nxi * igam0 + 
                    (ne_loc - 1) * nxi * ngam * (iabun0 - 1)];
      y3 = emission[ie + (ne_loc - 1) * ixi + (ne_loc - 1) * nxi * igam0 + 
                    (ne_loc - 1) * nxi * ngam * iabun0];
      y4 = emission[ie + (ne_loc - 1) * ixi + (ne_loc - 1) * nxi * (igam0 - 1) + 
                    (ne_loc - 1) * nxi * ngam * iabun0];
      flux0[ie + ne_loc * ixi] = (utmp1 * (ttmp1 * y1 + ttmp * y2) + 
                                  utmp * (ttmp * y3 + ttmp1 * y4));
    }
/*******************************************************************************
  ixi = 8;
  fprintf(stdout, "%d %d %d\n", iabun0, igam0, ixi);
  fprintf(stdout, "%f %f %f %f\n", abun[iabun0 - 1], gam[igam0 - 1], 
          pow(10, logxi[ixi - 1]), xi[ixi - 1]);
  for(ie = 0; ie < ne_loc - 1; ie++) 
    fprintf(stdout, "%d %f %f\n", ie, energy0[ie], 
            flux0[ie + ne_loc * (ixi - 1)] / (energy0[i + 1] - energy0[i]));
*******************************************************************************/
// Let's rebin the spectra to the evenly log spaced energies
  if ((rix != rix_old) || (strcmp(kydir, FGMSTR(pname)) != 0)) {
    energy2[0] = elow / (1. + pow(ehigh / elow, 1. / (nener - 1.))) * 2.;
    for (ie = 1; ie <= nener; ie++) {
      energy2[ie] = energy2[0] * pow(ehigh / elow, ie / (nener - 1.));
      energy1[ie - 1] = elow * pow(ehigh / elow, (ie - 1.) / (nener - 1.));
      if(thermalisation != 0.)
//      the factor rg^2 is due to the fact that we integrate in
//      dS = r * dr * dphi which we do in radius in geometrical units!!!
//      the flux (photon number) is per keV, per s and per cm^2
//      we multiply by Dnorm *  RG2 * mass2 later on, after ide routine is done
        flx[ie-1] = 2. * pow( energy1[ie - 1] / H_KEVS / C_CMS, 2) / 
                    H_KEVS / pow(f_col, 4);
    }
  }
  for (ixi = 0; ixi < nxi; ixi++)
    for (ie = 0; ie < nener; ie++) flux1[ie + nener * ixi] = 0.;
  ie = 1;
  while (energy0[ie] <= energy2[0]) ie++;
  je = 1;
  quit = 0;
  while ((ie <= (ne_loc - 1)) && (energy0[ie - 1] < energy2[nener])) {
    ttmp = energy0[ie] - energy0[ie - 1];
    while (energy2[je - 1] < energy0[ie]) {
      if (energy2[je - 1] < energy0[ie - 1]) utmp = energy0[ie - 1];
      else utmp = energy2[je - 1];
      if (energy2[je] < energy0[ie]) utmp1 = energy2[je];
      else utmp1 = energy0[ie];
      if (utmp1 > utmp)
        for (ixi = 0; ixi < nxi; ixi++) 
          flux1[je - 1 + nener * ixi] +=
            flux0[ie - 1 + ne_loc * ixi] * (utmp1 - utmp) / ttmp;
      if (je < nener) je++;
      else {
        quit = 1;
        break;
      }
    }
    if (quit) break;
    if (energy2[je - 1] > energy0[ie]) je--;
    ie++;
  }
  for (ixi = 0; ixi < nxi; ixi++)
    for (ie = 0; ie < nener; ie++) 
      flux1[ie + nener * ixi] /= (energy2[ie + 1] - energy2[ie]);
/*******************************************************************************
  ixi = 8;
  fprintf(stdout, "%d %d %d\n", iabun0, igam0, ixi);
  fprintf(stdout, "%f %f %f %f\n", abun[iabun0 - 1], gam[igam0 - 1], 
          pow(10, logxi[ixi - 1]), xi[ixi - 1]);
  for(ie = 0; ie < nener; ie++) 
    fprintf(stdout, "%d %f %f\n", ie, energy1[ie], flux1[ie + nener * (ixi - 1)]);
*******************************************************************************/
  abun_old = abundance;
  gam_old = gamma0;
}
sprintf(kydir, "%s", FGMSTR(pname));
rix_old = rix;
am_old = am;
thetaO_old = thetaO;
h_rh_old = h_rh;
/******************************************************************************/
//check if energy ranges defined are reasonable
for(iband=0;iband<nbands;iband++){
  if(ener_low[iband] < energy1[0] || ener_high[iband] > energy1[nener-1]){
    sprintf(errortxt,"kynrefrev: defined energy band reaches beyond the energy of local flux definition!");
    xs_write(errortxt, 5);
  }
}
// delay_phi
ide_param[19]=del_a;
/******************************************************************************/
//Let's compute the Np according to its definition, i.e. transform from 
//intrinsic or observed in 2-10keV to total intrinsic photon flux
//Lx is intrinsic photon flux in 2-10keV
g_L = sqrt(1. - 2. * h / (am2 + h * h));
if( Np < 0. ){
  Lx = -Np;
  Np *= - incgamma(2. - gamma0, E0 / EC) / 
        ( incgamma(2. - gamma0, 2. / EC) - incgamma(2. - gamma0, 10. / EC));
}else{
  Lx = Np / g_L / g_L / transf_o /
       ( incgamma(2. - gamma0, 2. / g_L / EC) - incgamma(2. - gamma0, 10. / g_L / EC));
  Np = Lx * incgamma(2. - gamma0, E0 / EC);        
  Lx *= ( incgamma(2. - gamma0, 2. / EC) - incgamma(2. - gamma0, 10. / EC));
}
//Let's write the intrinsic photon flux in 2-10keV into the xset XSPEC variable
//KYLxLamp
sprintf(kyLxLamp, "%e", Lx);
FPMSTR(pkyLxLamp, kyLxLamp);
// Let's compute the ionisation at rin and rout
for (i = 0; i < 2; i++) {
  if (i == 0)
    if (rin <= (radius[0] + r_plus)) r = radius[1] + r_plus;
    else r = rin;
  else r = rout;
  imin = 0;
  imax = nradius;
  ir0 = nradius / 2;
  while ((imax - imin) > 1) {
    if (r >= (radius[ir0 - 1] + r_plus)) imin = ir0;
    else imax = ir0;
    ir0 = (imin + imax) / 2;
  }
  if ((ir0 == nradius) && (r == (radius[nradius - 1] + r_plus))) ir0 = ir0 - 1;
  ttmp = (r - radius[ir0 - 1] - r_plus) / (radius[ir0] - radius[ir0 - 1]);
  ttmp1 = 1. - ttmp;
// Let's interpolate gfactor between two radii
  gfactor0 = ttmp * gfac[ir0] + ttmp1 * gfac[ir0 - 1];
// Let's interpolate lensing between two radii
  lensing = (ttmp * transf_d[ir0] + ttmp1 * transf_d[ir0 - 1]) * 
            h * sqrt(1. - 2. * h / (h * h + am2)) / (r * r + h * h) / r;
  if (lensing != 0.) {
    if (qn != 0.) {
      if (sw == 1) ionisation = LOGXI_NORM0 + log10(pow(r, -qn) * lensing * 
                                gfactor0 * Np / mass / nH0);
      //sw==2      
      else ionisation = LOGXI_NORM0 + log10(pow(r, -qn) * lensing * 
                        pow(gfactor0, gamma0 - 1.) * Np / mass / nH0);
    }
    else {
      if (sw == 1) ionisation = LOGXI_NORM0 + 
                                log10(lensing * gfactor0 * Np / mass / nH0);
      //sw==2      
      else ionisation = LOGXI_NORM0 + 
                 log10(lensing * pow(gfactor0, gamma0 - 1.) * Np / mass / nH0);
    }
    if (i == 0) {
      sprintf(kyxiin, "%e", pow(10, ionisation));
      FPMSTR(pkyxiin, kyxiin);
    }
    else {
      sprintf(kyxiout, "%e", pow(10, ionisation));
      FPMSTR(pkyxiout, kyxiout);
    }
  }
}
//******************************************************************************
//******************************************************************************
//******************************************************************************
// Let's integrate local emission over the accretion disc
if (ide(ear, ne, nt, far, qar, uar, var, ide_param, emis_KYNrefrev, 
        nener)) {
  for (ie = 0; ie < ne_xspec; ie++) photar[ie] = 0.;
  goto error;
}
//******************************************************************************
//******************************************************************************
// Let's normalize the reflected flux properly
for (ie = 0; ie < ne; ie++)
  for (it = 0; it <nt; it++)
    far[ie+ne*it] *= Dnorm * mass2 * RG2;
// Let's add primary flux to the solution (note we multiply by dt later on)
refl_ratio=-1.;
if (NpNr != 0.) {
//*
//  Let's use our own incomplete gamma function computations
  Anorm = LEDD * mass * ERG * Np / pow(EC, 2. - gamma0) / PI2 / 2. / 
          incgamma( 2.-gamma0, E0 / EC );
  gf[0] = incgamma(1.-gamma0, ear[0]/(g_L*zzshift*EC));
  for ( ie = 1; ie <= ne; ie++){
    gf[ie] = incgamma(1.-gamma0, ear[ie]/(g_L*zzshift*EC));
    photar1[ie-1]= ( gf[ie-1]-gf[ie] );
  }
  Anorm *= Dnorm * fabs(NpNr) * transf_o * pow(g_L*zzshift, gamma0) * pow( g_L*zzshift*EC, 1.-gamma0);
//*/
/*
// Let's compute the incomplete gamma function with the XSPEC incgamma 
// function
  Anorm = LEDD * mass * ERG * Np / pow(EC, 2. - gamma0) / PI2 / 2. / 
          incgamma(2. - gamma0, E0 / EC);
// let's compute the cut-off powerlaw with the XSPEC routine cutoffPowerLaw
  for (ie = 0; ie <= ne; ie++) ear1[ie] = ear[ie];
  param1[0] = param[9];
  param1[1] = g_L * zzshift * EC;
  cutoffpl(ear1, ne, param1, photar1);
  Anorm *= Dnorm * fabs(NpNr) * transf_o * pow(g_L*zzshift, gamma0);
*/
// let's compute the primary flux
  flux_refl_tot = flux_prim_tot = 0.;
  for (ie = 0; ie < ne; ie++){
    for (it = 0; it < nt; it++) flux_refl_tot += far[ie+ne*it];
    far_prim[ie]=0.;
  }
  for ( ie=0; ie < ne; ie++){
    if(ear[ie+1] > g_L*zzshift*E0){
      far_prim[ie] = Anorm * photar1[ie];
      flux_prim_tot += far_prim[ie];
    }
  }
  flux_prim_tot *= (flare_duration_rg / dt);
  refl_ratio = flux_refl_tot / flux_prim_tot;
}
sprintf(kyRefl, "%e", refl_ratio);
FPMSTR(pkyRefl, kyRefl);
/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// let's write the input parameters to a file
sprintf(text,"kynrefrev_%s.txt",cparam);
fw = fopen(text, "w");
fprintf(fw, "a/m         %12.6f\n", param[0]);
fprintf(fw, "thetaO      %12.6f\n", param[1]);
fprintf(fw, "rin         %12.6f\n", param[2]);
fprintf(fw, "ms          %12d\n", (int) param[3]);
fprintf(fw, "rout        %12.6f\n", param[4]);
fprintf(fw, "phi         %12.6f\n", param[5]);
fprintf(fw, "dphi        %12.6f\n", param[6]);
fprintf(fw, "M/M8        %12.6f\n", param[7]);
fprintf(fw, "height      %12.6f\n", param[8]);
fprintf(fw, "PhoIndex    %12.6f\n", param[9]);
fprintf(fw, "L/Ledd      %12.6f\n", param[10]);
fprintf(fw, "Np:Nr       %12.6f\n", param[11]);
fprintf(fw, "nH0         %12.6f\n", param[12]);
fprintf(fw, "qn          %12.6f\n", param[13]);
fprintf(fw, "abun        %12.6f\n", param[14]);
fprintf(fw, "therm       %12.6f\n", param[15]);
fprintf(fw, "arate       %12.6f\n", param[16]);
fprintf(fw, "f_col       %12.6f\n", param[17]);
fprintf(fw, "alpha       %12.6f\n", ide_param[21]);
fprintf(fw, "beta        %12.6f\n", ide_param[22]);
fprintf(fw, "rcloud      %12.6f\n", ide_param[23]);
fprintf(fw, "zshift      %12.6f\n", param[21]);
fprintf(fw, "limb        %12d\n", (int) param[22]);
fprintf(fw, "tab         %12d\n", (int) param[23]);
fprintf(fw, "sw          %12d\n", (int) param[24]);
fprintf(fw, "duration_rg %12.6f\n", DURATION);
fprintf(fw, "deltaT      %12.6f\n", dt0);
fprintf(fw, "ntbin       %12d\n", nt0);
fprintf(fw, "ntable      %12d\n", (int) param[25]);
fprintf(fw, "nradius     %12d\n", (int) param[26]);
fprintf(fw, "division    %12.6f\n", ide_param[8]);
fprintf(fw, "nphi        %12d\n", (int) param[28]);
fprintf(fw, "nt          %12d\n", nt_ratio);
fprintf(fw, "t1/f1/E1    %12.4g\n", param[31]);
fprintf(fw, "t2/f2/E2    %12.4g\n", param[32]);
fprintf(fw, "E3          %12.4g\n", param[33]);
fprintf(fw, "E4          %12.4g\n", param[34]);
fprintf(fw, "Af          %12.4g\n", param[35]);
fprintf(fw, "qf          %12.4g\n", param[36]);
fprintf(fw, "photar_sw   %12d\n", (int) param[37]);
fprintf(fw, "smooth      %12d\n", SMOOTH);
fprintf(fw, "Stokes      %12d\n", stokes);
fprintf(fw, "polar       %12d\n", polar);
fprintf(fw, "r_horizon   %12.6f\n", r_plus);
fprintf(fw, "r_ms        %12.6f\n", rms);
fprintf(fw, "edivision   %12d\n", (int) ide_param[14]);
fprintf(fw, "nener       %12ld\n", nener);
fprintf(fw, "norm        %12.6f\n", param[39]);
fprintf(fw, "nthreads    %12d\n", (int) param[38]);
fprintf(fw, "Xi_in       %s\n", kyxiin);
fprintf(fw, "Xi_out      %s\n", kyxiout);
fprintf(fw, "refl_ratio  %12.4g\n", refl_ratio);
fclose(fw);
#endif
/******************************************************************************/

// let's increase the time bin due to finate time resolution of the detector
for(it=0; it<nt/nt_ratio; it++){
  for(ie=0; ie<ne; ie++){
    far_final[ie+ne*it]=0.;
    for(i=0; i<nt_ratio; i++) far_final[ie+ne*it] += far[ie+ne*(nt_ratio*it+i)];
// we finaly use flux per second, i.e. we have to divide by dt * nt_ratio,
// however, we have not multiplied the flux by dt, so we divide only by nt_ratio
    far_final[ie+ne*it] /= nt_ratio;
  }
}
//lets return to the final time bin and deltaT
nt /= nt_ratio;
dt *= nt_ratio;

// final spectrum output
for(it=0;it<nt;it++)
  for(ie=0;ie<ne;ie++) far[ie+ne*it]=far_final[ie+ne*it];

// light curve
//itmin=1;
for(it=0;it<nt;it++){
  flux[it]=0.;
  for(ie=0;ie<ne;ie++){
//  far[ie+ne*it] is flux per energy bin not per energy!!!
//  therefore it must not be multiplied by the energy bin
//  when integrated into the light curve
//  light curve in the units of "counts per second"
    flux[it] += far[ie+ne*it];
// light curve in the units of "energy per second"
//    flux[it] += far[ie+ne*it]*(ear[ie]+ear[ie+1])/2.;
  }
//  if(itmin==it && flux[it]==0.)itmin++;
}
if(NpNr != 0.){
  flux_prim=0.;  
  for(ie=0;ie<ne;ie++) flux_prim += far_prim[ie];
}
// time integrated spectra divided by duration of the flare
for(ie=0;ie<ne;ie++){
  spectrum[ie]=0.;
  for(it=0;it<nt;it++) spectrum[ie] += far[ie+ne*it];
// we have to divide by duration of the flare, flare_duration_rg, however,
// we did not multiply by deltaT, thus we have to divide by
// flare_duration_rg / (deltaT)
  spectrum[ie] *= dt / flare_duration_rg;
  if(NpNr != 0.)spectrum_prim[ie] = far_prim[ie];
}

//let's define the time array in seconds
deltaT = dt * SEC * mass;
t0 = (t0 + 0.5*(nt_ratio-1.)/nt_ratio*dt) * SEC * mass;
for(it=0;it<nt;it++) time[it] = (it - 1) * deltaT + t0;
//For the FFT we need t0 to be the first value of time[] array
t0 = time[0];

/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// far is flux per energy bin and per unit time, i.e not per unit time bin
if(abs(photar_sw)<=10){
  sprintf(text,"kynrefrev_%s_far.dat",cparam);
  fw = fopen(text, "w");
  for(ie=0;ie<ne;ie++){
//    fprintf(fw, "%e\t", 0.5*(ear[ie]+ear[ie+1]));
    for(it=0;it<nt;it++){
//    for(it=itmin;it<nt;it++){
      fprintf(fw, "%e\t", far[ie+ne*it]/(ear[ie+1]-ear[ie]));
    }
    fprintf(fw, "\n");
  }
  fclose(fw);
}
sprintf(text,"kynrefrev_%s_lc.dat",cparam);
fw = fopen(text, "w");
for(it=0;it<nt;it++){
//  for(it=itmin;it<nt;it++){
  fprintf(fw, "%e\t%e\n", time[it], flux[it]);
}
fclose(fw);
if(NpNr != 0.){
  sprintf(text,"kynrefrev_%s_flux_prim.dat",cparam);
  fw = fopen(text, "w");
  fprintf(fw, "%e", flux_prim);
  fclose(fw);
  if(abs(photar_sw)<=10){
    sprintf(text,"kynrefrev_%s_spectrum.dat",cparam);
    fw = fopen(text, "w");
    for(ie=0;ie<ne;ie++){
      fprintf(fw, "%e\t%e\t%e\n", 0.5*(ear[ie]+ear[ie+1]),
              spectrum[ie]/(ear[ie+1]-ear[ie]),
              spectrum_prim[ie]/(ear[ie+1]-ear[ie]));
    }
    fclose(fw);
  }else{
    sprintf(text,"kynrefrev_%s_spectrum.dat",cparam);
    fw = fopen(text, "w");
    for(ie=0;ie<ne;ie++){
      fprintf(fw, "%e\t%e\n", 0.5*(ear[ie]+ear[ie+1]),
              spectrum[ie]/(ear[ie+1]-ear[ie]));
    }
    fclose(fw);
  }
}
#endif
/******************************************************************************/
if(photar_sw){
// Let's add the tail to the light curve only if rout=1000
  for(it=nt;it<nn;it++){
    time[it] = time[nt-1] + deltaT*(it-nt+1);
//    flux[it] = flux[nt-1] * pow((time[it]/time[nt-1]), -2.);
//    for(ie=0;ie<ne;ie++)
//      far[ie+ne*it] = far[ie+ne*(nt-1)] * pow((time[it]/time[nt-1]), -2.);
  }
  if( param[4] != 1000. ){
    for(it=nt;it<nn;it++){
      flux[it] = 0.;
      for(ie=0;ie<ne;ie++)
        far[ie+ne*it] = 0.;
    }
  }else{
//  let's use least squares on log(flux) vs. log(time) from last n points
    npoints=30.;
    sumt = sumf = sumt2 = sumtf = 0.;
    for(it=0;it<npoints;it++){
      ttmp = log10(time[nt-npoints+it]);
      utmp = log10(flux[nt-npoints+it]);
      sumt += ttmp;
      sumf += utmp;
      sumt2 += ttmp * ttmp;
      sumtf += ttmp * utmp;
    }
    qtime = ( npoints * sumtf - sumt * sumf ) /
            ( npoints * sumt2 - sumt * sumt );
    Atime = pow(10., ( sumf - qtime * sumt ) / npoints );
    for(it=nt;it<nn;it++){
      flux[it] = Atime * pow( time[it], qtime );
      if( flux[it] < 0. || isnan(flux[it]) || isinf(flux[it]) ) flux[it] = 0.;
    }
    for(ie=0;ie<ne;ie++){
      sumt = sumf = sumt2 = sumtf = 0.;
      for(it=0;it<npoints;it++){
        ttmp = log10(time[nt-npoints+it]);
        utmp = log10(far[ie+ne*(nt-npoints+it)]);
        sumt += ttmp;
        sumf += utmp;
        sumt2 += ttmp * ttmp;
        sumtf += ttmp * utmp;
      }
      qtime = ( npoints * sumtf - sumt * sumf ) /
              ( npoints * sumt2 - sumt * sumt );
      Atime = pow(10., ( sumf - qtime * sumt ) / npoints );
      for(it=nt;it<nn;it++){
        far[ie+ne*it] = Atime * pow( time[it], qtime );
        if( far[ie+ne*it] < 0. || isnan(far[ie+ne*it]) || isinf(far[ie+ne*it]) )
          far[ie+ne*it] = 0.;
      }
    }
  }
// Compute the flux in different energy bands
  for(iband=0;iband<nbands;iband++){
    for(it=0;it<nn;it++) flux_bands[iband+nbands*it] = 0.;
    flux_bands_prim[iband]=0.;
  }
  for(ie=0;ie<ne;ie++){
    for(iband=0;
        iband < ( ( exclude_energy && abs(photar_sw) >= 15 ) ? nbands-1 : nbands ); 
        iband++){
//    whole bins
      if(ear[ie] >= ener_low[iband] && ear[ie+1] <= ener_high[iband])ttmp=1.;
//    first bin
      else if(ear[ie] < ener_low[iband] && ear[ie+1] > ener_low[iband]){
        if(ear[ie+1] > ener_high[iband])ttmp=ener_high[iband]-ener_low[iband];
        else ttmp=ear[ie+1]-ener_low[iband];
        ttmp /= (ear[ie+1]-ear[ie]);
//    last bin
      }else if(ear[ie+1] > ener_high[iband] && ear[ie] < ener_high[iband] && 
               ear[ie] >= ener_low[iband])
        ttmp = (ener_high[iband]-ear[ie]) / (ear[ie+1]-ear[ie]);
      else ttmp=0;
      if(ttmp){
        for(it=0;it<nn;it++)
          flux_bands[iband+nbands*it] += far[ie+ne*it] * ttmp;
        flux_bands_prim[iband] += far_prim[ie] * ttmp;
      }
    }
//  reference energy band (iband=nbands-1) if we exclude the energy band of interest
    if( exclude_energy && abs(photar_sw) >= 15 ){
     iband=nbands-1;
//    whole bins
      if( ( ear[ie] >= ener_low[iband] && ear[ie+1] <= ener_high[iband]) &&
          ( ear[ie] >= ener_high[0] || ear[ie+1] <= ener_low[0] ) ) ttmp=1.;
//    partial bins
      else{
        en1 = ( ear[ie] > ener_low[iband] ? ear[ie] : ener_low[iband]);
        en1 = ( en1 > ener_high[0] ? en1 : ( en1 >= ener_low[0] ? ener_high[0] : en1 ) );
        en2 = ( ear[ie+1] < ener_high[iband] ? ear[ie+1] : ener_high[iband]);
        en2 = ( en2 < ener_low[0] ? en2 : ( en2 <= ener_high[0] ? ener_low[0] : en2 ) );
        ttmp=en2-en1;
        if(ttmp<0.)ttmp=0.;
        else ttmp /= (ear[ie+1]-ear[ie]);
      }
      if(ttmp){
        for(it=0;it<nn;it++)
          flux_bands[iband+nbands*it] += far[ie+ne*it] * ttmp;
        flux_bands_prim[iband] += far_prim[ie] * ttmp;
      }
    }
  }
/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// write the flux in energy bands light curves
  sprintf(filename, "kynrefrev_%s_bands_lc.dat",cparam);
  fw = fopen(filename, "w");
  for(it=0;it<nn;it++){
    fprintf(fw, "%lg\t", time[it]);
    for(iband=0;iband<nbands;iband++){
      fprintf(fw, "%E\t", flux_bands[iband+nbands*it]);
    }
    fprintf(fw, "\n");
  }
  fclose(fw);
// write the primary flux in energy bands
  sprintf(filename, "kynrefrev_%s_bands_prim.dat",cparam);
  fw = fopen(filename, "w");
  for(iband=0;iband<nbands;iband++)
    fprintf(fw, "%E\t", flux_bands_prim[iband]);
  fclose(fw);
#endif
/******************************************************************************/
  tot_time_sec = nn*deltaT;
// Compute the Fourier transform
  fft_reverberation(ear_xspec, ne_xspec, photar, frequency1, frequency2, photar_sw, 
                    time, nn, nt, far, far_prim, flux_bands, flux_bands_prim, 
                    nbands, tot_time_sec, flare_duration_sec, cparam);
  sprintf(kyfwrap, "%e", fwrap);
  FPMSTR(pkyfwrap, kyfwrap);
  if(photar_sw == 5 || photar_sw == 7) for(ie=0;ie<ne;ie++)photar[ie] *= ampl_ampl;
  if(photar_sw == 6 || photar_sw == 8) for(ie=0;ie<ne;ie++)photar[ie] += lag_shift * (ear[ie+1]-ear[ie]);
  if(photar_sw == 16 || photar_sw == 18) for(ie=0;ie<ne_xspec;ie++){
    if(qf == 1)photar[ie] += Af*log(ear_xspec[ie+1]/ear_xspec[ie]);
    else photar[ie] += Af*(pow(ear_xspec[ie+1],1.-qf)-pow(ear_xspec[ie],1.-qf))/(1.-qf);
  }
}else{
// interface with XSPEC
  if(time1 == 0.){
    for(ie=0;ie<ne;ie++){
      photar[ie] = spectrum[ie];
      if(NpNr != 0.)photar[ie] += spectrum_prim[ie];
//      photar[ie] /= flare_duration_rg;  <-- we have already divided by flare duration!
    }
  }else if(time1 >= time2){
//  given time1, find the corresponding index in time[]:
    it0 = (int) ceil( (time1 - time[0]) / deltaT + 1 );
    if(it0 < 1) it0 = 1;
    else if(it0 > nt-1) it0 = nt-1;
    ttmp = (time1 - time[it0 - 1]) / (time[it0] - time[it0 - 1]);
    ttmp1 = 1. - ttmp;
    for(ie=0;ie<ne;ie++){
      y1 = far[ie+ne*(it0-1)];
      y2 = far[ie+ne*it0];
      photar[ie] = (ttmp1 * y1 + ttmp * y2);
    }
    if(NpNr != 0. && time1_rg <= flare_duration_rg ) photar[ie] += far_prim[ie];
  }else{
//  given time1 and time2, find the corresponding index in time[]:
    it0 = (int) ceil( (time1 - time[0]) / deltaT );
    if(it0 < 1) it0 = 1;
    else if(it0 > nt-1) it0 = nt;
    itn = (int) ceil( (time2 - time[0]) / deltaT );
    if(itn < 1) itn = 0;
    else if(itn > nt-1) itn = nt-1;
//  the first and last partial bin
    ttmp = (time1 - time[it0 - 1]) / (time[it0] - time[it0 - 1]);
    ttmp1 = 1. - ttmp;
    utmp = (time2 - time[itn - 1]) / (time[itn] - time[itn - 1]);
    utmp1 = 1. - utmp;
    for(ie=0;ie<ne;ie++){
      y1 = far[ie+ne*(it0-1)];
      y2 = far[ie+ne*it0];
      y3 = far[ie+ne*(itn-1)];
      y4 = far[ie+ne*itn];
      if( it0 < itn ){
        photar[ie] = ( (ttmp1 * y1 + ttmp * y2) + far[ie+ne*it0] ) / 2. *
                     ( time[it0] - time1 ) / deltaT;
        photar[ie] += ( far[ie+ne*(itn-1)] + (utmp1 * y3 + utmp * y4) ) / 2. *
                     ( time2 - time[itn-1] ) / deltaT;
      }else if( it0 == itn ){
        photar[ie] = ( (ttmp1 * y1 + ttmp * y2) + (utmp1 * y3 + utmp * y4) ) / 2.
                      * ( time2 - time1 ) / deltaT;
      }
//    all the whole bins
      for(it=it0+1;it<itn;it++) photar[ie] += (far[ie+ne*it]+far[ie+ne*(it-1)])/2.;
      photar[ie] *= dt / flare_duration_rg;
      if(NpNr != 0.)photar[ie] += far_prim[ie];
    }
  }
}
/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// photar output
sprintf(filename, "kynrefrev_photar_%s.dat",cparam);
fw = fopen(filename, "w");
for (ie = 0; ie < ne_xspec; ie++) 
  if(abs(photar_sw) <=10) fprintf(fw, "%E\t%E\n", 0.5*(ear[ie]+ear[ie+1]), 
                                  photar[ie] / (ear[ie+1] - ear[ie]));
  else fprintf(fw, "%E\t%E\n", 0.5*(ear_xspec[ie]+ear_xspec[ie+1]), 
                             photar[ie] / (ear_xspec[ie+1] - ear_xspec[ie]));
fclose(fw);
#endif
/******************************************************************************/

error:
//Let's free all allocated arrays
if (time != NULL){free((void *) time); time = NULL;}
if (far != NULL){free((void *) far); far = NULL;}
if (far_final != NULL){free((void *) far_final); far_final = NULL;}
if (flux != NULL){free((void *) flux); flux = NULL;}
if (flux_bands != NULL){free((void *) flux_bands); flux_bands = NULL;}
if (spectrum != NULL){free((void *) spectrum); spectrum = NULL;}
if (far_prim != NULL){free((void *) far_prim); far_prim = NULL;}
if (spectrum_prim != NULL){free((void *) spectrum_prim); spectrum_prim = NULL;}
if (photar1 != NULL){free((void *) photar1); photar1 = NULL;}
if (gf != NULL){free((void *) gf); gf = NULL;}
////if (ear1 != NULL){free((void *) ear1); ear1 = NULL;}

return;
}
/*******************************************************************************
*******************************************************************************/
void emis_KYNrefrev(double** ear_loc, const int ne_loc, const long nt, 
                    double *far_loc, double *qar_loc, double *uar_loc, 
                    double *var_loc, const double r, const double phi, 
                    const double cosmu, const double phiphoton, 
                    const double alpha_o, const double beta_o, 
                    const double delay, const double g) {
// local emissivity --> far_loc(:) array
// local energy array --> ear_loc()
// ne_loc --> number of points in local energy where the local photon flux
// density in keV is defined;
// disc surface in polar coords r, phi;
// cosine of local emission angle --> cosmu

double factor, factor1, gfactor, lensing, ionisation, fluxe[ne_loc];
double temp, x, Ccal, Lcal, Fbb, Finc, Frefl, Ftherm, temp_new, flux_thermal;
double time, delay0;
double ttmp, ttmp1, y1, y2;
long   ixi0, it;
int    ie, imin, imax, ir0;
char   errortext[255];

*ear_loc = energy1;
// Normalization due to an imposed emissivity law
if (limb == 0) factor = 1. / PI;
if (limb == 1) factor = 1. / PI / (1. + 2.06 * 2. / 3.) * (1. + 2.06 * cosmu);
if (limb == 2) factor = 1. / PI * log(1. + 1. / cosmu);
factor *= nH0 * REFLIONX_NORM;
// given r, find corresponding indices in radius:
imin = 0;
imax = nradius;
ir0 = nradius / 2;
while ((imax - imin) > 1) {
  if (r >= (radius[ir0 - 1] + r_plus)) imin = ir0;
  else imax = ir0;
  ir0 = (imin + imax) / 2;
}
//if (ir0 == 0) ir0 = 1;
//if ((imax == nradius) && (r > (radius[nradius - 1] + r_plus)) ir0 = nradius;
if ((ir0 == nradius) && (r == (radius[nradius - 1] + r_plus))) ir0 -= 1;
if ((ir0 == 0) || (ir0 >= nradius)) {
  for(it=0;it<nt;it++) far_loc[ne_loc+(ne_loc+1)*it]=0.;
}else  {
  ttmp = (r - radius[ir0 - 1] - r_plus) / (radius[ir0] - radius[ir0 - 1]);
  ttmp1 = 1. - ttmp;
// Let's interpolate delay between two radii
  delay0 = ttmp * del[ir0] + ttmp1 * del[ir0 - 1];
    if((delay0+delay-t0 > nt*dt) || (t0-delay0-delay > flare_duration_rg)){
    for(it=0; it < nt; it++) far_loc[ne_loc+(ne_loc+1)*it]=0.;
  }else{
//  Let's interpolate gfactor between two radii
    gfactor = ttmp * gfac[ir0] + ttmp1 * gfac[ir0 - 1];
//  Let's interpolate cosmu0 between two radii
//    cosmu0 = ttmp * cosin[ir0] + ttmp1 * cosin[ir0 - 1];
//  Let's interpolate lensing between two radii
    lensing = (ttmp * transf_d[ir0] + ttmp1 * transf_d[ir0 - 1]) * 
              h * sqrt(1. - 2. * h / (h * h + am2)) / (r * r + h * h) / r;
// Let's compute the emitted flux at the particular radius
    if (lensing != 0.) {
      if(sw == 1) Finc = Np * LEDD / RG2 / mass * lensing * gfactor;
      else Finc = Np * LEDD / RG2 / mass * lensing * pow(gfactor, gamma0 - 1.);
      if (qn != 0.) ionisation = log10( Finc / ( nH0 * 1e15 * pow(r, qn) ) );
      else ionisation = log10( Finc / ( nH0 * 1e15 ) );
//    given ionisation, find the corresponding index in logXi():
      imin = 0;
      imax = nxi;
      ixi0 = nxi / 2;
      while ((imax - imin) > 1) {
        if (ionisation >= logxi[ixi0 - 1]) imin = ixi0;
        else imax = ixi0;
        ixi0 = (imin + imax) / 2;
      }
      if (ixi0 == 0) ixi0 = 1;
      ttmp = (ionisation - logxi[ixi0 - 1]) / (logxi[ixi0] - logxi[ixi0 - 1]);
      factor1 = 1.;
      if (ttmp < 0.) {
        ttmp = 0.;
        factor1 = pow(10, ionisation - logxi[0]);
      }
      if (ttmp > 1.) {
        ttmp = 1.;
        factor1 = pow(10, ionisation - logxi[nxi - 1]);
      }
      ttmp1 = 1. - ttmp;
      for (ie = 0; ie < ne_loc; ie++) {
        y1 = flux1[ie + ne_loc * (ixi0 - 1)];
        y2 = flux1[ie + ne_loc * ixi0];
        fluxe[ie] = (ttmp1 * y1 + ttmp * y2) * factor * factor1;
        if (sw == 1) fluxe[ie] *= pow(gfactor, gamma0 - 2.);
        if (qn != 0.) fluxe[ie] *= pow(r, qn);
      }
// add the thermal component
      if(thermalisation != 0.){
//      compute thermalised total flux (Ftherm = Finc - Frefl 
//      or Ftherm = thermalisation * Finc)...
//      incident flux in SI units ( i.e. J / s / m^2 )
        Finc *= ( 1e-3 / 4. / PI );
        if( fabs(thermalisation) <= 1. ) Ftherm = fabs(thermalisation) * Finc;
        else{
//        reflected flux in SI units
          Frefl=0.;
          ttmp = fluxe[0] * (*(*ear_loc));
          for (ie = 1; ie < ne_loc; ie++){
            ttmp1 = fluxe[ie] * (*(*ear_loc+ie));
            Frefl += (ttmp + ttmp1) / 2. * ((*(*ear_loc+ie))-(*(*ear_loc+ie-1)));
            ttmp = ttmp1;
          }
          Frefl *= KEV * 1e4;
//        fprintf(stdout,"%lg\n", Frefl / Finc);
          if(Finc < Frefl){
            sprintf(errortext, 
              "kynrefrev: reflection is larger than illumination at radius %lg GM/c^2!",
              r);
            xs_write(errortext, 5);
            Frefl=Finc;
          }
          Ftherm = Finc - Frefl;
        }
//      compute the temperature of the accretion disc
        x = sqrt(r);
        Ccal = 1. - 3. / (x*x) + 2. * am / (x * x * x);
        if(am < 1.)
          Lcal = 1. / x * (x - x0 - 1.5 * am * log(x / x0) - 
                 3. * pow(x1 - am,2.)/x1/(x1 - x2)/(x1 - x3) * log((x - x1)/(x0 - x1))-
                 3. * pow(x2 - am,2.)/x2/(x2 - x1)/(x2 - x3) * log((x - x2)/(x0 - x2))-
                 3. * pow(x3 - am,2.)/x3/(x3 - x1)/(x3 - x2) * log((x - x3)/(x0 - x3)));
        else Lcal = 1. / x * (x - 1 + 1.5 * log((x + 2.) / 3. / x));
        temp = Tnorm * pow(x, -1.5) * pow( arate / mass2 * Lcal / Ccal, 0.25 );
//      compute new temperature
        Fbb = SIGMA * pow( temp / K_KEVK / f_col, 4. );
        temp_new = K_KEVK * pow( ( Fbb + Ftherm ) / SIGMA, 0.25 ) * f_col;
//      compute the additional flux due to thermalisation (i.e. "above"
//      the stationary thermal flux)
        for(ie = 0; ie < ne_loc; ie++){
          flux_thermal = flx[ie] * ( 1. / (exp((*(*ear_loc+ie)) / temp_new) - 1.) 
                                   - 1. / (exp((*(*ear_loc+ie)) / temp) - 1.) );
          if( thermalisation < 0. ) fluxe[ie] = flux_thermal;
          else fluxe[ie] += flux_thermal;
        }
      }
    }else{
      for(it=0;it<nt;it++) far_loc[ne_loc+(ne_loc+1)*it]=0.;
      return;
    }
    for(it=0;it<nt;it++){
      time = it * dt - delay0 - delay + t0 - dt*nt_ratio;
      if(time>=0. && time<=flare_duration_rg){
        far_loc[ne_loc+(ne_loc+1)*it]=1.;
        for(ie=0;ie<ne_loc;ie++) far_loc[ie+(ne_loc+1)*it]=fluxe[ie];
      }else far_loc[ne_loc+(ne_loc+1)*it]=0.;
    }
  }
}
/*******************************************************************************
// local spectrum output -- write energy1[] and far_local[] into file
  fw = fopen("kynrefrev_photar_loc.dat", "w");
  for (ie = 0; ie < ne_loc; ie++)
    fprintf(fw, "%14.6f\t%E\n", energy1[ie], fluxe[ie]);
  fclose(fw)
*******************************************************************************/
return;
}
/*******************************************************************************
*******************************************************************************/
double incgamma(double a, double x){

double gfunc;
long   i;

gfunc=0.;
for(i=100000;i>=1;i--) gfunc=(a-i)*i/(2.*i+1.+x-a+gfunc);
gfunc = pow(x,a) * exp(-x) / (1.+x-a+gfunc);
return(gfunc);
}
/*******************************************************************************
*******************************************************************************/
