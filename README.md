Table of contents
=================

  * [Model description](#model-description)
  * [Installation](#installation)
    * [Required files](#required-files)
    * [Usage in XSPEC](#usage-in-xspec)
    * [Usage outside of XSPEC](#usage-outside-of-xspec)
  * [Parameters of the model](#parameters-of-the-model)
    * [Definition in XSPEC](#definition-in-xspec)
    * [Definition outside XSPEC](#definition-outside-xspec)
  * [Output files created](#output-files-created)

Model description
=================

The KYNrefrev model computes the time dependent reflection spectra of the disc 
as a response to a flash of primary power-law radiation from a point source 
located on the axis of the black-hole accretion disc.

_Assumptions of the model:_

* central Kerr black hole,

* Keplerian, geometrically thin, optically thick, ionised disc with different 
  radial density profiles,

* stationary hot point-like patch of plasma located on the system rotation axis 
  and emitting isotropic power-law radiation,

* full relativistic ray-tracing code in vacuum is used for photon paths from the 
  corona to the disc and to the observer and from the disc to the observer,

* re-processing in the ionised accretion disc is computed for each radius from 
  REFLIONX tables for constant density slab illuminated by power-law radiation,

* the ionisation of the disc is set for each radius according to the amount of 
  the incident primary flux and the density of the accretion disc,

* several limb brightening/darkening prescriptions for directionality of the 
  re-processed emission are used.

_Output of the code:_

* time dependent spectra (only when used outside of XSPEC) of the disc response 
  and observed primary flash,

* integrated spectrum,

* light curve for a given energy band,

* lag as a function of frequency between given energy bands,

* lag as a function of energy for different frequencies.


Installation
============

Required files
--------------

* Source files in the main repository directory.

* KY tables: [KBHlamp_qt.fits](https://owncloud.asu.cas.cz/index.php/s/xg64GRMSRGiWOPR) 
  (also [here](http://www.astro.cas.cz/dovciak/pub/KY/KBHlamp_qt.fits)) 
  and [KBHtables80.fits](https://owncloud.asu.cas.cz/index.php/s/WP8aLN168MJgcB9) 
  (also [here](http://www.astro.cas.cz/dovciak/pub/KY/KBHtables80.fits)).

* [REFLION(X)](https://heasarc.gsfc.nasa.gov/xanadu/xspec/models/reflion.html) 
  tables (unpack gzipped files): 

   - [reflion.mod](https://heasarc.gsfc.nasa.gov/xanadu/xspec/models/reflion.mod.gz) (old),
   - [reflionx.mod](https://heasarc.gsfc.nasa.gov/xanadu/xspec/models/reflionx.mod.gz),

   or in case the links are not available or if the tables there are updated and 
   their format/structure has changed:

   - [reflion.mod](https://owncloud.asu.cas.cz/index.php/s/6CWcb0o5Ssjehju) 
     (or [here](http://www.astro.cas.cz/dovciak/pub/KY-external/reflion.mod))
     (old),
   - [reflionx.mod](https://owncloud.asu.cas.cz/index.php/s/Q6biiTPM1QBMtiT)
     (or [here](http://www.astro.cas.cz/dovciak/pub/KY-external/reflionx.mod)).

Usage in XSPEC
--------------

The code is compiled inside XSPEC with the following command (assuming all the 
source files and FITS tables are in the directory /path/to/KYNrefrev):

* `initpackage kynrefrev lmodel.dat /path/to/KYNrefrev`.

To use the KYNrefrev model inside XSPEC, first the package needs to be loaded 
and directory with KYNrefrev set:

* `lmod kynrefrev /path/to/KYNrefrev`,

* `xset KYDIR /path/to/KYNrefrev`.

Then the model may be used:

* `mo kynrefrev`.

_Note_: 
In case of segmentation fault, one may need to increase the stack size, e.g. 
with the command `ulimit -s unlimited` or `ulimit -s 65532`.


Usage outside of XSPEC
----------------------

* One also needs the Makefile and libxspec library included in the directory 
  'other'.

* The library to work with FITS files (libcfitsio.so) is needed, thus one needs 
  to define the name of the library and path to it in the provided Makefile.

* The model parameters have to be changed inside the source file.

* Compile with the make command:

    * `make kynrefrev`

* Run the code:

    * `./kynrefrev`.

* The model creates various files described below.

_Note_:
In case of segmentation fault, one may need to increase the stack size, e.g. 
with the command `ulimit -s unlimited` or `ulimit -s 65532`.


Parameters of the model
=======================

Definition in XSPEC
-------------------

   * The meaning of the input parameters are also explained at the beginning of 
     the kynrefrev.c file. The parameters when the code runs under XSPEC are 
     defined in the usual way as for other XSPEC models. The parameter 
     definitions when run outside of XSPEC must be changed directly inside the 
     source code. Summary of the parameters:

     * **par1  ... a/M**
       - black hole angular momentum (-1 &le; a/M &le; 1)
     * **par2  ... theta_o**
       - observer inclination in degrees (0&deg;-pole, 90&deg;-disc)
     * **par3  ... rin**
       - inner edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
     * **par4  ... ms**
       - switch for inner edge
       - 0: we integrate from inner edge = par3 
       - 1: if the inner edge of the disc is below marginally stable orbit (MSO) 
            then we integrate emission above MSO only
       - 2: we integrate from inner edge given in units of MSO, i.e. inner 
            edge = par3 &times; r~mso~ (the same applies for outer edge)
     * **par5  ... rout**
       - outer edge of non-zero disc emissivity (in GM/c^2 or in r~mso~)
     * **par6  ... phi**
       - lower azimuth of non-zero disc emissivity (degrees)
     * **par7  ... dphi**
       - (phi + dphi) is upper azimuth of non-zero disc emissivity 0&deg; &le; 
         dphi &le; 360&deg;
     * **par8 ... M/M8**
       - black hole mass in units of 10^8 solar masses
     * **par9 ... height**
       - height on the axis (measured from the center) at which the primary 
         source is located (GM/c^(2))
     * **par10 ... PhoIndex**
       - power-law energy index of the primary flux
     * **par11 ... Np**
       - dN/dt/d&Omega;, the intrinsic local (if negative) or the observed 
         (if positive) primary isotropic flux in the X-ray energy range 2-10keV 
         in units of L~Edd~
     * **par12 ... NpNr**
       - ratio of the primary to the reflected normalization
       - 1: self-consistent model for isotropic primary source
       - 0: only reflection, primary source is hidden
       - if positive then Np (par11) means the luminosity towards the observer
       - if negative then Np (par11) means the luminosity towards the disc
     * **par13 ... nH0**
       - density profile normalization in 10^15 cm^(-3)
     * **par14 ... q_n**
       - radial power-law density profile
     * **par15 ... abun**
       - Fe abundance (in solar abundance)
     * **par16 ... alpha**
       - position of the cloud centre in GM/c^2 in alpha coordinate (alpha being 
         the impact parameter in &phi;-direction, positive for approaching side 
         of the disc)
     * **par17 ... beta**
       - position of the cloud centre in GM/c^2 in beta coordinate (beta being 
         the impact parameter in &theta;-direction, positive in up direction, 
         i.e. away from the disc)
     * **par18 ... rcloud**
       - radius of the obscuring cloud
       - the meaning of cloud is inverted for negative values of rcloud, i.e. 
         only the radiation behind the cloud is computed
     * **par19 ... zshift**
       - overall Doppler shift
     * **par20 ... limb**
       - 0: for isotropic emission (flux ~ 1)
       - 1: for Laor's limb darkening (flux ~ 1+2.06&mu;)
       - 2: for Haardt's limb brightening (flux ~ ln (1+1/&mu;))
     * **par21 ... tab**
       - which reflion table to use
       - 1: reflion (the old one, lower cut-off energy at 1eV, not good for 
            PhoIndex > 2)
       - 2: reflionx (the newer one, lower cut-off energy at 100eV)
     * **par22 ... sw**
       - switch for the way how to compute the refl. spectra
       - 1: use the computed ionisation parameter, &xi;, for the interpolation 
            in reflion, i.e. use proper total incident intensity with the 
            shifted cut-offs
       - 2: use the ionisation parameter, &xi;, correspondent to the computed 
            normalization of the incident flux, i.e. do not shift the cut-offs 
            when computing the total incident intensity
     * **par23 ... ntable**
       - defines fits file with tables (0 &le; ntable &le; 99), currently the 
         tables with ntable=80 are correct for this model
     * **par24 ... nradius**
       - number of grid points in radius
       - if negative than the number of radial grid points is dependent on 
         height as 
       -nradius&nbsp;/&nbsp;height^(&nbsp;0.66) 
     * **par25 ... division**
       - type of division in radial integration
       - 0: equidistant radial grid (constant linear step)
       - 1: exponential radial grid (constant logarithmic step)
       - >1: mixed radial grid with a constant logarithmic step in the inner 
           region and with a constant linear step in the outer region; the 
           total nradius (par24) number of points is divided in the 3:2 ratio 
           in these regions; the value of par25 gives the transition radius 
           between these regions (in GM/c^(2))
       - -1: mixed radial grid with the transition radius at 2&times;height
     * **par26 ... nphi**
       - number of grid points in azimuth
     * **par27 ... deltaT**
       - length of the time bin (GM/c^(3))
     * **par28 ... nt**
       - number of time subbins per one time bin
     * **par29 ... t1/f1/E1**
       - the time to be used in XSPEC for the spectrum (0 means average 
         spectrum, i.e. divided by the flare duration)
       - the frequency to be used in XSPEC for the energy dependent Fourier 
         transform (0 means average values in the range of 0 to the first 
         wrapping frequency)
       - positive values are in sec or Hz
       - negative values are in GM/c^3 or (GM/c^(3))^(-1)
       - if different than par30, the value gives the lower end of the 
         time/frequency interval of interest
       - if same as par30, then the functions are computed for this value of 
         the time/frequency of interest
       - in case of frequency dependent lags it defines the lower value of the 
         energy band of interest in keV
     * **par30 ... t2/f2/E2**
       - used only if different than par29 and if par29 is nonzero
       - its value gives the upper end of the time/frequency interval of 
         interest
       - positive values are in sec or Hz
       - negative values are in GM/c^3 or (GM/c^(3))^(-1)
       - in case of frequency dependent lags it defines the upper value of the 
         energy band of interest in keV
     * **par31 ... E3**
       - it defines the lower value of the reference energy band for lag or 
         amplitude energy dependence as well as in case of frequency dependent 
         lags and amplitudes
       - if zero no reference band is used
       - if negative, the whole energy band is used as a reference band for 
         lag-energy spectra, always excluding the current energy bin; it must be
         non-negative in case of lag-frequency dependence
     * **par32 ... E4**
       - it defines the upper value of the reference energy band for lag-energy
         dependence as well as in case of frequency dependent lags
     * **par33 ... tshift/Af**
       - lag shift for lag-energy dependence in case of par35=+6
       - multiplicative factor in case of adding empirical hard lags 
         Af&times;f^(qf), used for par35=+16
     * **par34 ... Amp/qf**
       - multiplicative factor for the amplitude-energy dependence in case of 
         par35=+5
       - powerlaw index in case of adding empirical hard lags Af&times;f^(qf), 
         used for par35=+16
     * **par35 ... photar_sw**
       - defines output in the XSPEC (photar array)
       - 0: spectrum for time interval defined by par29 and par30
       - _the following values correspond to energy dependent Fourier transform 
         at the frequency band defined by par29 and par30:_
         - -1: real part of FT of the relative reflection
         - -2: imaginary part of FT of the relative reflection
         - -3: amplitude of FT of the relative reflection
         - -4: phase of FT of the relative reflection
         - -5: amplitude for the relative reflection divided by amplitude in the 
              reference energy band defined by par31 and par32
         - -6: lag for the relative reflection with respect to reference energy 
              band defined by par31 and par32
         - 1: real part of FT including primary radiation
         - 2: imaginary part of FT including primary radiation
         - 3: amplitude of FT including primary radiation
         - 4: phase of FT including primary radiation
         - 5: amplitude including the primary radiation divided by amplitude in 
              the reference energy band defined by par31 and par32
         - 6: lag diluted by primary radiation with respect to reference energy 
              band defined by par31 and par32
       - _the following values correspond to frequency dependent Fourier 
         transform for the energy band of interest defined by par29 and par30:_
         - -11: real part of FT of the relative reflection
         - -12: imaginary part of FT of the relative reflection
         - -13: amplitude of FT of the relative reflection
         - -14: phase of FT of the relative reflection
         - -15: amplitude  for the relative reflection divided by amplitude in 
              the reference energy band defined by par31 and par32
         - -16: lag for the relative reflection with respect to reference energy 
              band defined by par31 and par32
         - 11: real part of FT including primary radiation
         - 12: imaginary part of FT including primary radiation
         - 13: amplitude of FT including primary radiation
         - 14: phase of FT including primary radiation
         - 15: amplitude including the primary radiation divided by amplitude in 
              the reference energy band defined by par31 and par32
         - 16: lag diluted by primary radiation with respect to reference energy 
              band defined by par31 and par32
     * **par36 ... nthreads**
       - how many threads should be used for computations
     * **par37 ... norm**
       - **has to be set to unity!**

Definition outside XSPEC
------------------------

   * The model parameters need to be defined inside the **kynrefrev.c** code 
     when run outside of XSPEC (the code needs to be recompiled after changing 
     them):

       - _energy_ in the following lines:

                  #define NE     30
                  #define E_MIN  0.3
                  #define E_MAX  80.

       - choose the _energy bands of interest_ in the following lines:

                  #define NBANDS 5
                  .
                  .
                  .
                  //definition of energy band of interest, reference band is defined as the last 
                  //one, usually the whole energy range
                  ener_low[0] = 0.3;
                  ener_high[0] = 0.8;
                  ener_low[1] = 1.;
                  ener_high[1] = 3.;
                  ener_low[2] = 3.;
                  ener_high[2] = 9.;
                  ener_low[3] = 12.;
                  ener_high[3] = 40.;
                  ener_low[4] = E_MIN;
                  ener_high[4] = E_MAX;

       - _all basic parameters_  of the model (physical ones as well as those 
         defining resolution grid for computations) are defined in the following 
         lines:

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
                  param[10] = 0.001;    // Np
                  param[11] = 1.;       // NpNr
                  param[12] = 1.;       // nH0
                  param[13] = 0.;       // q_n
                  param[14] = 1.;       // abun
                  param[15] = -6.;      // alpha
                  param[16] = 0.;       // beta
                  param[17] = 0.;       // rcloud
                  param[18] = 0.;       // zshift
                  param[19] = 0.;       // limb
                  param[20] = 2.;       // tab
                  param[21] = 2.;       // sw
                  param[22] = 80.;      // ntable
                  param[23] = -4488.;   // nradius
                  param[24] = -1.;      // division
                  param[25] = 180.;     // nphi
                  param[26] = 1.;       // deltaT
                  param[27] = 1.;       // nt
                  param[28] = 2.e-4;    // time/frequency/energy-lower
                  param[29] = 8.e-4;    // time/frequency/energy-upper
                  param[30] = -1.;      // reference energy band-lower
                  param[31] = 3.;       // reference energy band-upper
                  param[32] = 0.;       // lag shift or multiplicative factor for hard lags
                  param[33] = 1.;       // amplitude multiplicative factor or power-law index for hard lags
                  param[34] = 6.;       // photar_sw
                  param[35] = 4.;       // nthreads
                  param[36] = 1.;       // norm

       - some parameters are later changed in the loops for convenience (to 
         create files for grid of parameters), see lines as:

                  for (ia=0;ia<=1;ia++){
                    param[0] = (double) ia;
                    for (iinc=20;iinc<=80;iinc+=20){
                      param[1] = (double) iinc;
                      for (ih=1;ih<=20;ih++){
                        param[8] = 1.5 * (100./1.5)**((ih-1.)/19.);


Output files created
====================

   * The output files are created only when the code is run outside XSPEC.
   * The following naming scheme is used for the output files:
     - **AAA** is 100&times; the horizon value (thus 100 means a=1 and 200 means 
       a=0),
     - **BB** is the inclination in degrees,
     - **CCCC** is 10&times; the height (e.g. 0030 means h=3),
     - **u1** is used for phase unwrapped in frequency dependence,
     - **u2** is used for phase unwrapped in energy dependence.
   * All spectra and light curves are always computed in photon numbers and per 
     keV and per second, time is in seconds and frequency in Hz.
   * The output files created by **kynrefrev.c** code:
     - below by relative reflection it is meant the disc response divided by the 
       total primary flux in the flash.
     - **kynrefrev_photar.dat** &rarr; the values as would be given inside 
       XSPEC,
     - **kynrefrev_AAA_BB_CCCC.txt** &rarr; the summary of parameter values,
     - **kynrefrev_AAA_BB_CCCC_far.dat** &rarr; the time evolving observed 
       reflection spectrum where the energy changes with rows and the time 
       changes with columns,
     - **kynrefrev_AAA_BB_CCCC_flux_prim.dat** &rarr; the total observed primary 
       flux (i.e. integrated in energy) per second, it is constant during the 
       duration of the flare,
     - **kynrefrev_AAA_BB_CCCC_lc.dat** &rarr; the light curve of the observed 
       reflection (2^nd column) where we integrated over the whole energy range, 
       the 1^st column contains the time,
     - **kynrefrev_AAA_BB_CCCC_spectrum.dat** &rarr; the time integrated 
       spectrum of the observed reflection (2^nd column) and the observed 
       primary (3^rd column), both are divided by the flare duration, the 1^st 
       column contains the central value of the energy bins in keV.
     - **kynrefrev_AAA_BB_CCCC_bands_lc.dat** &rarr; the light curves for the 
       observed reflection (2^nd and higher columns) where in each column the 
       light curve is integrated over different energy band (defined in the 
       code), the 1^st column contains the time,
     - **kynrefrev_AAA_BB_CCCC_bands_prim.dat** &rarr; the observed primary flux 
       per second, it is constant during the duration of the flare where in each 
       column the flux is integrated over different energy band (defined in the 
       code).
     - **kynrefrev_AAA_BB_CCCC_real.dat** 
       &rarr; the real part of the FFT of the relative reflection with frequency 
       changing with rows and energy with columns,
     - **kynrefrev_AAA_BB_CCCC_imag.dat** 
       &rarr; the imaginary part of the FFT of the relative reflection with 
       frequency changing with rows and energy with columns,
     - **kynrefrev_AAA_BB_CCCC_ampl.dat** 
       &rarr; the amplitude of the FFT of the relative reflection with frequency 
       changing with rows and energy with columns,
     - **kynrefrev_AAA_BB_CCCC_phase.dat**, 
       **kynrefrev_AAA_BB_CCCC_phase_u1.dat**,
       **kynrefrev_AAA_BB_CCCC_phase_u2.dat** &rarr; the phase of the FFT of the 
       relative reflection with frequency changing with rows and energy with 
       columns,
     - **kynrefrev_AAA_BB_CCCC_real_tot.dat** &rarr; the real part of the FFT of 
       the total signal (reflection response plus primary flash) with frequency 
       changing with rows and energy with columns,
     - **kynrefrev_AAA_BB_CCCC_imag_tot.dat** &rarr; the imaginary part of the 
       FFT of the total signal (reflection response plus primary flash) with 
       frequency changing with rows and energy with columns,
     - **kynrefrev_AAA_BB_CCCC_ampl_tot.dat** &rarr; the amplitude of the FFT of 
       the total signal (reflection response plus primary flash) with frequency 
       changing with rows and energy with columns,
     - **kynrefrev_AAA_BB_CCCC_phase_tot.dat**, 
       **kynrefrev_AAA_BB_CCCC_phase_tot_u1.dat**,
       **kynrefrev_AAA_BB_CCCC_phase_tot_u2.dat** &rarr; the phase of the FFT of 
       the total signal (reflection response plus primary flash) with frequency 
       changing with rows and energy with columns,
     - **kynrefrev_AAA_BB_CCCC_bands_real.dat** &rarr; the real part, as a 
       function of frequency, of the FFT of the relative reflection integrated 
       in different energy bands as defined in the code (2^nd and higher 
       columns), the 1^st column contains the frequency,
     - **kynrefrev_AAA_BB_CCCC_bands_imag.dat** &rarr; the imaginary part, as a 
       function of frequency, of the FFT of the relative reflection integrated 
       in different energy bands as defined in the code (2^nd and higher 
       columns), the 1^st column contains the frequency,
     - **kynrefrev_AAA_BB_CCCC_bands_ampl.dat** &rarr; the amplitude, as a 
       function of frequency, of the FFT of the relative reflection integrated 
       in different energy bands as defined in the code (2^nd and higher 
       columns), the 1^st column contains the frequency,
     - **kynrefrev_AAA_BB_CCCC_bands_phase.dat**, 
       **kynrefrev_AAA_BB_CCCC_bands_phase_u1.dat** &rarr; the phase, as a 
       function of frequency, of the FFT of the relative reflection integrated 
       in different energy bands as defined in the code (2^nd and higher 
       columns), the 1^st column contains the frequency,
     - **kynrefrev_AAA_BB_CCCC_bands_real_tot.dat** &rarr; the real part, as a 
       function of frequency, of the FFT of the total signal (reflection 
       response plus primary flash) integrated in different energy bands as 
       defined in the code (2^nd and higher columns), the 1^st column contains 
       the frequency,
     - **kynrefrev_AAA_BB_CCCC_bands_imag_tot.dat** &rarr; the imaginary part, 
       as a function of frequency, of the FFT of the total signal (reflection 
       response plus primary flash) integrated in different energy bands as 
       defined in the code (2^nd and higher columns), the 1^st column contains 
       the frequency,
     - **kynrefrev_AAA_BB_CCCC_bands_ampl_tot.dat** &rarr; the amplitude, as a 
       function of frequency, of the FFT of the total signal (reflection 
       response plus primary flash) integrated in different energy bands as 
       defined in the code (2^nd and higher columns), the 1^st column contains 
       the frequency,
     - **kynrefrev_AAA_BB_CCCC_bands_phase_tot.dat**, 
       **kynrefrev_AAA_BB_CCCC_bands_phase_tot_u1.dat** &rarr; the phase, as a 
       function of frequency, of the FFT of the total signal (reflection 
       response plus primary flash) integrated in different energy bands as 
       defined in the code (2^nd and higher columns), the 1^st column contains 
       the frequency,
     - **kynrefrev_AAA_BB_CCCC_freq_wrap.dat** &rarr; the first wrapping 
       frequency for the phase computed for the relative reflection (the lowest 
       one from all energy bins)
     - **kynrefrev_AAA_BB_CCCC_fft_tot_int.dat** &rarr; energy dependent average 
       values of the FFT of the total signal (real part, imaginary part, 
       amplitude, phase, unwrapped phase, delay, ratio of the amplitudes for the 
       energy band of interest and reference energy band, delay between the 
       energy band of interest and reference energy band computed from wrapped 
       and unwrapped phases and directly averaged delay between the two energy 
       bands as well as the ratio of the amplitudes and delay difference between 
       the energy band of interest and reference energy band where reference 
       energy band always excludes the current energy bin), the average is 
       computed in the range of 0 to the first wrapping frequency, the 1^st 
       column contains the central value of energy bins in keV; the FFT here is 
       averaged over frequencies for real and imaginary parts first and then, 
       from the result, all the rest quantities are computed except for the 
       directly averaged delay, where the delay is computed first from real and 
       imaginary parts of the FFT for each frequency and only then it is 
       averaged (just for comparison).
     - **kynrefrev_AAA_BB_CCCC_fft_tot_fband.dat** &rarr; energy dependent 
       average values of the FFT of the total signal (real part, imaginary part, 
       amplitude, phase, unwrapped phase, delay, ratio of the amplitudes for the 
       energy band of interest and reference energy band, delay between the 
       energy band of interest and reference energy band computed from wrapped 
       and unwrapped phases and directly averaged delay between the two energy 
       bands as well as the ratio of the amplitudes and delay difference between 
       the energy band of interest and reference energy band where reference 
       energy band always excludes the current energy bin), the average is 
       computed in the frequency range given by param[28] and param[29], the 
       1^st column contains the central value of energy bins in keV; the FFT 
       here is averaged over frequencies for real and imaginary parts first and 
       then, from the result, all the rest quantities are computed except for 
       the directly averaged delay, where the delay is computed first from real 
       and imaginary parts of the FFT for each frequency and only then it is 
       averaged (just for comparison).
