Table of contents
=================

  * [Model description](#model-description)
  * [Installation](#installation)
    * [Required files](#required-files)
    * [Usage in XSPEC](#usage-in-xspec)
    * [Usage outside of XSPEC](#usage-outside-of-xspec)
  * [Parameters of the KYNrefrev model](#parameters-of-the-kynrefrev-model)
    * [Definition in XSPEC](#definition-kynrefrev-in-xspec)
    * [Definition outside XSPEC](#definition-kynrefrev-outside-xspec)
  * [Parameters of the KYNxilrev model](#parameters-of-the-kynxilrev-model)
    * [Definition in XSPEC](#definition-kynxilrev-in-xspec)
    * [Definition outside XSPEC](#definition-kynxilrev-outside-xspec)
  * [Output files created](#output-files-created)


Model description
=================

The KYNrefrev and KYNxilrev models compute the time dependent reflection spectra 
of the disc as a response to a flash of primary power-law radiation from a point 
source located on the axis of the black-hole accretion disc.

_Assumptions of the model:_

* central Kerr black hole,

* Keplerian, geometrically thin, optically thick, ionised disc with different 
  radial density profiles,

* stationary hot point-like patch of plasma located on the system rotation axis 
  and emitting isotropic power-law radiation,

* full relativistic ray-tracing code in vacuum is used for photon paths from the 
  corona to the disc and to the observer and from the disc to the observer,

* re-processing in the ionised accretion disc is computed for each radius from 
  REFLIONX tables (KYNrefrev) or XILLVER tables (KYNxilrev) for constant density 
  slab illuminated by power-law radiation,

* increase in the disc temperature due to partial thermalisation of the 
  illuminating flux

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

* KY tables: [KBHlamp80.fits](https://owncloud.asu.cas.cz/index.php/s/abuFcygHKEKFiSa) 
  (also [here](http://www.astro.cas.cz/dovciak/pub/KY/KBHlamp80.fits)) 
  and [KBHtables80.fits](https://owncloud.asu.cas.cz/index.php/s/WP8aLN168MJgcB9) 
  (also [here](http://www.astro.cas.cz/dovciak/pub/KY/KBHtables80.fits)).

* [REFLION(X)](https://heasarc.gsfc.nasa.gov/xanadu/xspec/models/reflion.html) 
  tables for KYNrefrev (Ross & Fabian 2005, MNRAS, 358, 211) - unpack gzipped files:

    - [reflion.mod](https://heasarc.gsfc.nasa.gov/xanadu/xspec/models/reflion.mod.gz) (old),
    - [reflionx.mod](https://heasarc.gsfc.nasa.gov/xanadu/xspec/models/reflionx.mod.gz),
 
        or in case the links are not available or if the tables there are updated and 
        their format/structure has changed:
 
    - [reflion.mod](https://owncloud.asu.cas.cz/index.php/s/6CWcb0o5Ssjehju) 
      (or [here](http://www.astro.cas.cz/dovciak/pub/KY-external/reflion.mod))
      (old),
    - [reflionx.mod](https://owncloud.asu.cas.cz/index.php/s/Q6biiTPM1QBMtiT)
      (or [here](http://www.astro.cas.cz/dovciak/pub/KY-external/reflionx.mod)).
 
* [XILLVER](https://hea-www.cfa.harvard.edu/%7Ejavier/xillver/)
  tables for KYNxilrev (Garcia & Kallman 2010, ApJ, 718, 695, 
  Garcia et al. 2013, ApJ, 768, 2 and Garcia et al. 2016, MNRAS, 462, 751) - 
  unpack gzipped files:

    - [xillver-a-Ec5.fits](https://hea-www.cfa.harvard.edu/%7Ejavier/xillver/tables/xillver-a-Ec5.fits),
    - [xillverD-4.fits](http://www.sternwarte.uni-erlangen.de/~dauser/research/relxill/xillverD-4.fits.gz).

Usage in XSPEC
--------------

The code is compiled inside XSPEC with the following command (assuming all the 
source files and FITS tables are in the directory /path/to/KYNreverb):

* `initpackage kynreverb lmodel-kynreverb.dat /path/to/KYNreverb`.

To use the KYNrefrev or KYNxilrev model inside XSPEC, first the package needs to 
be loaded and directory with KYNreverb set:

* `lmod kynreverb /path/to/KYNreverb`,

* `xset KYDIR /path/to/KYNreverb`.

Then the models may be used:

* `mo kynrefrev`

or

* `mo kynxilrev`.

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

    or

    * `make kynxilrev`

* Run the code:

    * `./kynrefrev`

    or

    * `./kynxilrev`.

* The models create various files described below.

_Note_:
In case of segmentation fault, one may need to increase the stack size, e.g. 
with the command `ulimit -s unlimited` or `ulimit -s 65532`.


Parameters of the KYNrefrev model
=================================

Definition in XSPEC
-------------------

The meaning of the input parameters are also explained at the beginning of 
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
  * **par11 ... L/L~Edd~**
    - dE/dt, the intrinsic local (if negative) or the observed 
      (if positive) primary isotropic flux in the X-ray energy range 2-10keV 
      in units of L~Edd~
  * **par12 ... Np:Nr**
    - ratio of the primary to the reflected normalization
    - 1: self-consistent model for isotropic primary source
    - 0: only reflection, primary source is hidden
    - if positive then L/L~Edd~ (par11) means the luminosity towards the 
      observer
    - if negative then L/L~Edd~ (par11) means the luminosity towards the disc
  * **par13 ... density/ionisation**
    - density profile normalization in 10^15 cm^(-3) if positive
    - ionisation profile normalisation if it is negative
    - this parameter cannot be zero
  * **par14 ... den_prof/ion_prof**
    - radial power-law density profile if par13 is positive
    - radial ionisation profile if par13 is negative
    - the radial profiles in both cases are given by abs(par13) &times; 
      r^(par14)
  * **par15 ... abun**
    - Fe abundance (in solar abundance)
  * **par16 ... therm**
    - fraction of thermalised flux from the overal incident flux illuminating 
      the disc
    - = 0: only the reverberation of reflected radiation is computed
    - < 0: only the reverberation of thermal radiation is  computed
    - > 0: both the thermal and reflection reverberation is included
    - abs(par16) > 1: the fraction of thermalisation is computed from 
                      difference between the incident and reflected fluxes
  * **par17 ... arate**
    - accretion rate in units of L~Edd~ if positive or in Solar mass per  
      Julian year (365.25 days) if negative
  * **par18 ... f_col**  
    - spectral hardening factor
  * **par19 ... alpha**
    - position of the cloud centre in GM/c^2 in alpha coordinate (alpha being 
      the impact parameter in &phi;-direction, positive for approaching side 
      of the disc)
  * **par20 ... beta**
    - position of the cloud centre in GM/c^2 in beta coordinate (beta being 
      the impact parameter in &theta;-direction, positive in up direction, 
      i.e. above the disc)
  * **par21 ... rcloud**
    - radius of the obscuring cloud
    - the meaning of cloud is inverted for negative values of rcloud, i.e. 
      only the radiation transmitted through the cloud is computed
  * **par22 ... zshift**
    - overall Doppler shift
  * **par23 ... limb**
    - 0: for isotropic emission (flux ~ 1)
    - 1: for Laor's limb darkening (flux ~ 1+2.06&mu;)
    - 2: for Haardt's limb brightening (flux ~ ln (1+1/&mu;))
  * **par24 ... tab**
    - which reflion table to use
    - 1: reflion (the old one, lower cut-off energy at 1eV, not good for 
         PhoIndex > 2)
    - 2: reflionx (the newer one, lower cut-off energy at 100eV)
  * **par25 ... sw**
    - switch for the way how to compute the refl. spectra
    - 1: use the computed ionisation parameter, &xi;, for the interpolation 
         in reflion, i.e. use proper total incident intensity with the 
         shifted cut-offs
    - 2: use the ionisation parameter, &xi;, correspondent to the computed 
         normalization of the incident flux, i.e. do not shift the cut-offs 
         when computing the total incident intensity
  * **par26 ... ntable**
    - defines fits file with tables (0 &le; ntable &le; 99), currently the 
      tables with ntable=80 are correct for this model
  * **par27 ... nrad**
    - number of grid points in radius
    - if negative than the number of radial grid points is dependent on 
      height as -nrad&nbsp;/&nbsp;height^(&nbsp;0.66) 
  * **par28 ... division**
    - type of division in radial integration
    - 0: equidistant radial grid (constant linear step)
    - 1: exponential radial grid (constant logarithmic step)
    - >1: mixed radial grid with a constant logarithmic step in the inner 
        region and with a constant linear step in the outer region; the 
        total nradius (par27) number of points is divided in the 3:2 ratio 
        in these regions; the value of par28 gives the transition radius 
        between these regions (in GM/c^(2))
    - -1: mixed radial grid with the transition radius at 2&times;height
  * **par29 ... nphi**
    - number of grid points in azimuth
  * **par30 ... deltaT**
    - length of the time bin (GM/c^(3))
  * **par31 ... nt**
    - number of time subbins per one time bin
  * **par32 ... t1/f1/E1**
    - the time to be used in XSPEC for the spectrum (0 means average 
      spectrum, i.e. divided by the flare duration)
    - the frequency to be used in XSPEC for the energy dependent Fourier 
      transform (0 means average values in the range of 0 to the first 
      wrapping frequency)
    - positive values are in sec or Hz
    - negative values are in GM/c^3 or (GM/c^(3))^(-1)
    - if different than par33, the value gives the lower end of the 
      time/frequency interval of interest
    - if same as par33, then the functions are computed for this value of 
      the time/frequency of interest
    - in case of frequency dependent lags it defines the lower value of the 
      energy band of interest in keV
  * **par33 ... t2/f2/E2**
    - used only if different than par32 and if par32 is nonzero
    - its value gives the upper end of the time/frequency interval of 
      interest
    - positive values are in sec or Hz
    - negative values are in GM/c^3 or (GM/c^(3))^(-1)
    - in case of frequency dependent lags it defines the upper value of the 
      energy band of interest in keV
  * **par34 ... Eref1**
    - it defines the lower value of the reference energy band for lag or 
      amplitude energy dependence as well as in case of frequency dependent 
      lags and amplitudes
    - if zero no reference band is used
    - if negative:
        * for lag-energy spectra, the whole energy band is used as a reference 
          band, always excluding the current energy bin
        * for lag-frequency dependence, the energy reference band is
          abs(par34) to abs(par35) excluding overlaping part with energy band 
          of interest abs(par32) to abs(par33)
  * **par35 ... Eref2**
    - it defines the upper value of the reference energy band for lag-energy
      dependence as well as in case of frequency dependent lags
  * **par36 ... dt/Af**
    - lag shift for lag-energy dependence in case of par38=+6
    - multiplicative factor in case of adding empirical hard lags 
      Af&times;f^(qf), used for par38=+16 and par38=+18; 
      if par36=-1 then the following hard lags prescription is used (see 
      Epitropakis & Papadakis, 2017):
      100 * log10(E~ref~/E) * (f/1e-4)^(-1) s
      with E~ref~ being middle of the reference energy band and E middle of 
      the energy band of interest
  * **par37 ... Amp/qf**
    - multiplicative factor for the amplitude-energy dependence in case of 
      par38=+5
    - powerlaw index in case of adding empirical hard lags Af&times;f^(qf), 
      used for par38=+16 and par38=+18
  * **par38 ... xsw**
    - defines output in the XSPEC (photar array)
        - 0: spectrum for time interval defined by par32 and par33
    - _the following values correspond to energy dependent Fourier transform 
       at the frequency band defined by par32 and par33:_
        - -1: real part of FT of the relative reflection
        - -2: imaginary part of FT of the relative reflection
        - -3: amplitude of FT of the relative reflection
        - -4: phase of FT of the relative reflection
        - -5: amplitude for the relative reflection divided by amplitude in the 
              reference energy band defined by par34 and par35  (integration in 
              frequencies is done in real and imaginary parts first and then 
              the amplitudes are computed)
        - -6: lag for the relative reflection with respect to reference energy 
              band defined by par34 and par35 (integration in frequencies is 
              done in real and imaginary parts first and then the lags are 
              computed with frequency at half of the wrapping frequency or 
              middle of the frequency band)
        - -7: amplitude  for the relative reflection divided by amplitude in 
              the reference energy band defined by par34 and par35 (integration 
              in frequencies here is done in amplitudes directly)
        - -8: lag for the relative reflection with respect to reference energy 
              band defined by par34 and par35 (integration in frequencies here 
              is done in lags directly)
        - 1: real part of FT including primary radiation
        - 2: imaginary part of FT including primary radiation
        - 3: amplitude of FT including primary radiation
        - 4: phase of FT including primary radiation
        - 5: amplitude including the primary radiation divided by amplitude in 
             the reference energy band defined by par34 and par35 (integration 
             in frequencies is done in real and imaginary parts first and then 
             the amplitudes are computed)
        - 6: lag diluted by primary radiation with respect to reference energy 
             band defined by par34 and par35 (integration in frequencies is 
             done in real and imaginary parts first and then the lags are 
             computed with frequency at half of the wrapping frequency or 
             middle of the frequency band)
        - 7: amplitude including the primary radiation divided by amplitude in 
             the reference energy band defined by par34 and par35 (integration 
             in frequencies here is done in amplitudes directly)
        - 8: lag diluted by primary radiation with respect to reference energy 
             band defined by par34 and par35 (integration in frequencies here 
             is done in lags directly)
    - _the following values correspond to frequency dependent Fourier 
       transform for the energy band of interest defined by par32 and par33:_
        - -11: real part of FT of the relative reflection
        - -12: imaginary part of FT of the relative reflection
        - -13: amplitude of FT of the relative reflection
        - -14: phase of FT of the relative reflection
        - -15: amplitude  for the relative reflection divided by amplitude in 
               the reference energy band defined by par34 and par35 (rebinning 
               here is done in real and imaginary parts first and then the 
               amplitudes are computed)
        - -16: lag for the relative reflection with respect to reference energy 
               band defined by par34 and par35 (rebinning here is done in real 
               and imaginary parts first and then the lags are computed)
        - -17: amplitude  for the relative reflection divided by amplitude in 
               the reference energy band defined by par34 and par35 (rebinning 
               here is done in amplitudes directly)
        - -18: lag for the relative reflection with respect to reference energy 
               band defined by par34 and par35 (rebinning here is done in lags 
               directly)
        - 11: real part of FT including primary radiation
        - 12: imaginary part of FT including primary radiation
        - 13: amplitude of FT including primary radiation
        - 14: phase of FT including primary radiation
        - 15: amplitude including the primary radiation divided by amplitude in 
              the reference energy band defined by par34 and par35 (rebinning 
              here is done in real and imaginary parts first and then the 
              amplitudes are computed)
        - 16: lag diluted by primary radiation with respect to reference energy 
              band defined by par34 and par35 (rebinning here is done in real 
              and imaginary parts first and then the lags are computed)
        - 17: amplitude including the primary radiation divided by amplitude in 
              the reference energy band defined by par34 and par35 (rebinning 
              here is done in amplitudes directly)
        - 18: lag diluted by primary radiation with respect to reference energy 
              band defined by par34 and par35 (rebinning here is done in lags 
              directly)
  * **par39 ... nthreads**
    - how many threads should be used for computations
  * **par40 ... norm**
    - **has to be set to unity!**

Definition outside XSPEC
------------------------

The model parameters need to be defined inside the **kynrefrev.c** code 
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
             param[31] = 2.e-4;    // t1/f1/E1
             param[32] = 8.e-4;    // t2/f2/E2
             param[33] = -1.;      // Eref1
             param[34] = 3.;       // Eref2
             param[35] = 0.;       // dt/Af
             param[36] = 1.;       // Amp/qf
             param[37] = 6.;       // xsw
             param[38] = 4.;       // nthreads
             param[39] = 1.;       // norm

  - some parameters are later changed in the loops for convenience (to 
    create files for grid of parameters), see lines as:

             for (ia=0;ia<=1;ia++){
               param[0] = (double) ia;
               for (iinc=20;iinc<=80;iinc+=20){
                 param[1] = (double) iinc;
                 for (ih=1;ih<=20;ih++){
                   param[8] = 1.5 * (100./1.5)**((ih-1.)/19.);


Parameters of the KYNxilrev model
=================================

Definition in XSPEC
-------------------

The meaning of the input parameters are also explained at the beginning of 
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
  * **par11 ... L/L~Edd~**
    - dE/dt, the intrinsic local (if negative) or the observed 
      (if positive) primary isotropic flux in the X-ray energy range 2-10keV 
      in units of L~Edd~
  * **par12 ... Np:Nr**
    - ratio of the primary to the reflected normalization
    - 1: self-consistent model for isotropic primary source
    - 0: only reflection, primary source is hidden
    - if positive then L/L~Edd~ (par11) means the luminosity towards the 
      observer
    - if negative then L/L~Edd~ (par11) means the luminosity towards the disc
  * **par13 ... density/ionisation**
    - density profile normalization in 10^15 cm^(-3) if positive, 
      i.e. n = par13 &times; r^(par14)
    - ionisation profile normalisation if it is negative and  constant density 
      xillver tables are used, i.e. &xi; = -par13 &times; r^(par14)
    - ionisation parameter if it is negative and xillver tables dependent on 
      density are used, i.e. &xi; = -par13
    - this parameter cannot be zero
  * **par14 ... den_prof/ion_prof**
    - radial power-law density profile if par13 is positive
    - radial ionisation profile if par13 is negative and constant density 
      xillver tables are used
    - density in 10^15 cm^(-3) if par13 is negative and xillver tables
      dependent on density are used
  * **par15 ... abun**
    - Fe abundance (in solar abundance)
  * **par16 ... E_cut**
    - the observed (if positive) or intrinsic local at the source 
      (if negative) cut-off energy of the primary X-ray radiation
  * **par17 ... therm**
    - fraction of thermalised flux from the overal incident flux illuminating 
      the disc
    - = 0: only the reverberation of reflected radiation is computed
    - < 0: only the reverberation of thermal radiation is  computed
    - > 0: both the thermal and reflection reverberation is included
    - abs(par17) > 1: the fraction of thermalisation is computed from 
                      difference between the incident and reflected fluxes
  * **par18 ... arate**
    - accretion rate in units of L~Edd~ if positive or in Solar mass per  
      Julian year (365.25 days) if negative
  * **par19 ... f_col**
    - spectral hardening factor
  * **par20 ... alpha**
    - position of the cloud centre in GM/c^2 in alpha coordinate (alpha being 
      the impact parameter in &phi;-direction, positive for approaching side 
      of the disc)
  * **par21 ... beta**
    - position of the cloud centre in GM/c^2 in beta coordinate (beta being 
      the impact parameter in &theta;-direction, positive in up direction, 
      i.e. above the disc)
  * **par22 ... rcloud**
    - radius of the obscuring cloud
    - the meaning of cloud is inverted for negative values of rcloud, i.e. 
      only the radiation transmitted through the cloud is computed
  * **par23 ... zshift**
    - overall Doppler shift
  * **par24 ... limb**
    - only used for angle averaged XILLVER tables
    - 0: for isotropic emission (flux ~ 1)
    - 1: for Laor's limb darkening (flux ~ 1+2.06&mu;)
    - 2: for Haardt's limb brightening (flux ~ ln (1+1/&mu;))
  * **par25 ... tab**
    - which XILLVER table to use
    - 1: xillver.fits, angle averaged with cut-off energy at 
         300 keV
    - 2: xillver-a.fits, angle dependent with cut-off energy at
         300 keV
    - 3: xillver-Ec.fits, angle averaged with free cut-off energy 
    - 4: xillver-a-Ec.fits, angle dependent with free cut-off 
         energy 
    - 5: xillver-a-Ec2.fits, angle dependent with free cut-off
         energy 
    - 6: xillver-a-Ec3.fits, angle dependent with free cut-off
         energy 
    - 7: xillver-a-Ec4.fits, angle dependent with free cut-off
         energy 
    - 8: xillver-a-Ec5.fits, angle dependent with free cut-off
         energy 
    - 11: xillverD-4.fits, angle dependent with cut-off energy at 300 keV for
          disc density 10^(15)-10^(19) cm^(-3)
  * **par26 ... ntable**
    - defines fits file with tables (0 &le; ntable &le; 99), currently the 
      tables with ntable=80 are correct for this model
  * **par27 ... nrad**
    - number of grid points in radius
    - if negative than the number of radial grid points is dependent on 
      height as -nrad&nbsp;/&nbsp;height^(&nbsp;0.66) 
  * **par28 ... division**
    - type of division in radial integration
    - 0: equidistant radial grid (constant linear step)
    - 1: exponential radial grid (constant logarithmic step)
    - >1: mixed radial grid with a constant logarithmic step in the inner 
        region and with a constant linear step in the outer region; the 
        total nradius (par27) number of points is divided in the 3:2 ratio 
        in these regions; the value of par28 gives the transition radius 
        between these regions (in GM/c^(2))
    - -1: mixed radial grid with the transition radius at 2&times;height
  * **par29 ... nphi**
    - number of grid points in azimuth
  * **par30 ... deltaT**
    - length of the time bin (GM/c^(3))
  * **par31 ... nt**
    - number of time subbins per one time bin
  * **par32 ... t1/f1/E1**
    - the time to be used in XSPEC for the spectrum (0 means average 
      spectrum, i.e. divided by the flare duration)
    - the frequency to be used in XSPEC for the energy dependent Fourier 
      transform (0 means average values in the range of 0 to the first 
      wrapping frequency)
    - positive values are in sec or Hz
    - negative values are in GM/c^3 or (GM/c^(3))^(-1)
    - if different than par33, the value gives the lower end of the 
      time/frequency interval of interest
    - if same as par33, then the functions are computed for this value of 
      the time/frequency of interest
    - in case of frequency dependent lags it defines the lower value of the 
      energy band of interest in keV
  * **par33 ... t2/f2/E2**
    - used only if different than par32 and if par32 is nonzero
    - its value gives the upper end of the time/frequency interval of 
      interest
    - positive values are in sec or Hz
    - negative values are in GM/c^3 or (GM/c^(3))^(-1)
    - in case of frequency dependent lags it defines the upper value of the 
      energy band of interest in keV
  * **par34 ... Eref1**
    - it defines the lower value of the reference energy band for lag or 
      amplitude energy dependence as well as in case of frequency dependent 
      lags and amplitudes
    - if zero no reference band is used
    - if negative:
        * for lag-energy spectra, the whole energy band is used as a reference 
          band, always excluding the current energy bin
        * for lag-frequency dependence, the energy reference band is
          abs(par34) to abs(par35) excluding overlaping part with energy band 
          of interest abs(par32) to abs(par33)
  * **par35 ... Eref2**
    - it defines the upper value of the reference energy band for lag-energy
      dependence as well as in case of frequency dependent lags
  * **par36 ... dt/Af**
    - lag shift for lag-energy dependence in case of par38=+6
    - multiplicative factor in case of adding empirical hard lags 
      Af&times;f^(qf), used for par38=+16 and par38=+18; 
      if par36=-1 then the following hard lags prescription is used (see 
      Epitropakis & Papadakis, 2017):
      100 * log10(E~ref~/E) * (f/1e-4)^(-1) s
      with E~ref~ being middle of the reference energy band and E middle of 
      the energy band of interest
  * **par37 ... Amp/qf**
    - multiplicative factor for the amplitude-energy dependence in case of 
      par38=+5
    - powerlaw index in case of adding empirical hard lags Af&times;f^(qf), 
      used for par38=+16 and par38=+18
  * **par38 ... xsw**
    - defines output in the XSPEC (photar array)
        - 0: spectrum for time interval defined by par32 and par33
    - _the following values correspond to energy dependent Fourier transform 
      at the frequency band defined by par32 and par33:_
        - -1: real part of FT of the relative reflection
        - -2: imaginary part of FT of the relative reflection
        - -3: amplitude of FT of the relative reflection
        - -4: phase of FT of the relative reflection
        - -5: amplitude for the relative reflection divided by amplitude in the 
              reference energy band defined by par34 and par35  (integration in 
              frequencies is done in real and imaginary parts first and then 
              the amplitudes are computed)
        - -6: lag for the relative reflection with respect to reference energy 
              band defined by par34 and par35 (integration in frequencies is 
              done in real and imaginary parts first and then the lags are 
              computed with frequency at half of the wrapping frequency or 
              middle of the frequency band)
        - -7: amplitude  for the relative reflection divided by amplitude in 
              the reference energy band defined by par34 and par35 (integration 
              in frequencies here is done in amplitudes directly)
        - -8: lag for the relative reflection with respect to reference energy 
              band defined by par34 and par35 (integration in frequencies here 
              is done in lags directly)
        - 1: real part of FT including primary radiation
        - 2: imaginary part of FT including primary radiation
        - 3: amplitude of FT including primary radiation
        - 4: phase of FT including primary radiation
        - 5: amplitude including the primary radiation divided by amplitude in 
             the reference energy band defined by par34 and par35 (integration 
             in frequencies is done in real and imaginary parts first and then 
             the amplitudes are computed)
        - 6: lag diluted by primary radiation with respect to reference energy 
             band defined by par34 and par35 (integration in frequencies is 
             done in real and imaginary parts first and then the lags are 
             computed with frequency at half of the wrapping frequency or 
             middle of the frequency band)
        - 7: amplitude including the primary radiation divided by amplitude in 
             the reference energy band defined by par34 and par35 (integration 
             in frequencies here is done in amplitudes directly)
        - 8: lag diluted by primary radiation with respect to reference energy 
             band defined by par34 and par35 (integration in frequencies here 
             is done in lags directly)
    - _the following values correspond to frequency dependent Fourier 
      transform for the energy band of interest defined by par32 and par33:_
        - -11: real part of FT of the relative reflection
        - -12: imaginary part of FT of the relative reflection
        - -13: amplitude of FT of the relative reflection
        - -14: phase of FT of the relative reflection
        - -15: amplitude  for the relative reflection divided by amplitude in 
               the reference energy band defined by par34 and par35 (rebinning 
               here is done in real and imaginary parts first and then the 
               amplitudes are computed)
        - -16: lag for the relative reflection with respect to reference energy 
               band defined by par34 and par35 (rebinning here is done in real 
               and imaginary parts first and then the lags are computed)
        - -17: amplitude  for the relative reflection divided by amplitude in 
               the reference energy band defined by par34 and par35 (rebinning 
               here is done in amplitudes directly)
        - -18: lag for the relative reflection with respect to reference energy 
               band defined by par34 and par35 (rebinning here is done in lags 
               directly)
        - 11: real part of FT including primary radiation
        - 12: imaginary part of FT including primary radiation
        - 13: amplitude of FT including primary radiation
        - 14: phase of FT including primary radiation
        - 15: amplitude including the primary radiation divided by amplitude in 
              the reference energy band defined by par34 and par35 (rebinning 
              here is done in real and imaginary parts first and then the 
              amplitudes are computed)
        - 16: lag diluted by primary radiation with respect to reference energy 
              band defined by par34 and par35 (rebinning here is done in real 
              and imaginary parts first and then the lags are computed)
        - 17: amplitude including the primary radiation divided by amplitude in 
              the reference energy band defined by par34 and par35 (rebinning 
              here is done in amplitudes directly)
        - 18: lag diluted by primary radiation with respect to reference energy 
              band defined by par34 and par35 (rebinning here is done in lags 
              directly)
  * **par39 ... nthreads**
    - how many threads should be used for computations
  * **par40 ... norm**
    - **has to be set to unity!**

Definition outside XSPEC
------------------------

The model parameters need to be defined inside the **kynxilrev.c** code 
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
             param[10] = 0.001;    // L/Ledd
             param[11] = 1.;       // Np:Nr
             param[12] = 1.;       // density
             param[13] = 0.;       // den_prof
             param[14] = 1.;       // abun
             param[15] = 300.;     // E_cut
             param[16] = 0.;       // thermalisation
             param[17] = 0.1;      // arate
             param[18] = 2.4;      // f_col
             param[19] = -6.;      // alpha
             param[20] = 0.;       // beta
             param[21] = 0.;       // rcloud
             param[22] = 0.;       // zshift
             param[23] = 0.;       // limb
             param[24] = 11.;      // tab
             param[25] = 80.;      // ntable
             param[26] = -4488.;   // nrad
             param[27] = -1.;      // division
             param[28] = 180.;     // nphi
             param[29] = 1.;       // deltaT
             param[30] = 1.;       // nt
             param[31] = 2.e-4;    // t1/f1/E1
             param[32] = 8.e-4;    // t2/f2/E2
             param[33] = -1.;      // Eref1
             param[34] = 3.;       // Eref2
             param[35] = 0.;       // dt/Af
             param[36] = 1.;       // Amp/qf
             param[37] = 6.;       // xsw
             param[38] = 4.;       // nthreads
             param[39] = 1.;       // norm

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
