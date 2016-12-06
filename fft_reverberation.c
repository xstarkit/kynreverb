#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ISGN -1
#define PI   3.14159265358979
#define PI2  6.28318530717959

extern int    xs_write(char* wrtstr, int idest);
extern double t0, fwrap;
extern int    exclude_energy;

void fft_reverberation(const double *ear, int ne, double *photar, 
                       double frequency1, double frequency2, int photar_sw, 
                       double *time, long nn, long nt, 
                       double *far, double *far_prim,
                       double *flux_bands, double *flux_bands_prim,
                       int nbands, double tot_time_sec, 
                       double flare_duration_sec, char *cparam){

void four(double *data_r, double *data_im, long n1, int isign);

double freq_wrap=0., freq_wrap0=0.;
double Re_sinc, Im_sinc;
double *data_r=NULL, *data_im=NULL, *freq=NULL;
double *r_part=NULL, *im_part=NULL, *ampl=NULL, *phase=NULL, *phase_u1=NULL, 
       *phase_u2=NULL, *ampl_etot2=NULL, *phase_etot2=NULL;
double *r_part_tot=NULL, *im_part_tot=NULL, *ampl_tot=NULL, *phase_tot=NULL, 
       *phase_tot_u1=NULL, *phase_tot_u2=NULL, *ampl_tot_etot2=NULL, *phase_tot_etot2=NULL;
double *r_part_bands=NULL, *im_part_bands=NULL, *ampl_bands=NULL, 
       *phase_bands=NULL, *phase_bands_u1=NULL;
double *r_part_bands_tot=NULL, *im_part_bands_tot=NULL, *ampl_bands_tot=NULL,
       *phase_bands_tot=NULL, *phase_bands_tot_u1=NULL;
double r_part_int[ne], im_part_int[ne], ampl_int[ne], ampl2_int[ne],
       phase_int[ne], phase_int_u2[ne], delay2_int[ne];
double r_part_tot_int[ne], im_part_tot_int[ne], ampl_tot_int[ne], ampl2_tot_int[ne],
       phase_tot_int[ne], phase_tot_int_u2[ne], delay2_tot_int[ne];
double r_part_fband[ne], im_part_fband[ne], ampl_fband[ne], ampl2_fband[ne],
       phase_fband[ne], phase_fband_u2[ne], delay2_fband[ne];
double r_part_tot_fband[ne], im_part_tot_fband[ne], ampl_tot_fband[ne], ampl2_tot_fband[ne],
       phase_tot_fband[ne], phase_tot_fband_u2[ne], delay2_tot_fband[ne];
double r_part_int_etot=0., im_part_int_etot=0., ampl_int_etot=0., 
       ampl2_int_etot=0., ampl2_int_etot2[ne],
       phase_int_etot=0., delay2_int_etot=0., delay2_int_etot2[ne];
double r_part_tot_int_etot=0., im_part_tot_int_etot=0., ampl_tot_int_etot=0., 
       ampl2_tot_int_etot=0., ampl2_tot_int_etot2[ne], 
       phase_tot_int_etot=0., delay2_tot_int_etot=0., delay2_tot_int_etot2[ne];
double r_part_fband_etot=0., im_part_fband_etot=0., ampl_fband_etot=0., 
       ampl2_fband_etot=0., ampl2_fband_etot2[ne], 
       phase_fband_etot=0., delay2_fband_etot=0., delay2_fband_etot2[ne];
double r_part_tot_fband_etot=0., im_part_tot_fband_etot=0., ampl_tot_fband_etot=0., 
       ampl2_tot_fband_etot=0.,ampl2_tot_fband_etot2[ne],
       phase_tot_fband_etot=0., delay2_tot_fband_etot=0., delay2_tot_fband_etot2[ne];
double r_part_etot=0., im_part_etot=0., ampl_etot=0., phase_etot=0.;
double r_part_tot_etot=0., im_part_tot_etot=0., ampl_tot_etot=0., 
       phase_tot_etot=0.;
double r_part_bands_rebin[2*ne], im_part_bands_rebin[2*ne];
double ampl_bands_rebin[2*ne], phase_bands_rebin[2*ne];
double ampl2_bands_rebin[2*ne], delay2_bands_rebin[2*ne];
double r_part_bands_tot_rebin[2*ne], im_part_bands_tot_rebin[2*ne];
double ampl_bands_tot_rebin[2*ne], phase_bands_tot_rebin[2*ne];
double ampl2_bands_tot_rebin[2*ne], delay2_bands_tot_rebin[2*ne];
double ttmp, ttmp1, utmp, utmp1, y1=0., y2=0., y3=0., y4=0.;
int    ie, iband, iiband, k1, k2, nonwrap;
long   it, ifr, ifr0, ifrn, ifreq;
char   filename[255], errortxt[255];
FILE   *fw1, *fw2, *fw3, *fw4, *fw5, *fw6;

//Let's create all large arrays dynamically so that they are allocated on the 
//heap instead of on an available size-limited stack!
if ( ( data_r  = (double *) malloc( nn * sizeof(double) ) ) == NULL ||
     ( data_im = (double *) malloc( nn * sizeof(double) ) ) == NULL ) {
  xs_write("kynrefrev: Failed to allocate memory for tmp arrays.", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  goto error;
}
if ( ( freq  = (double *) malloc( (nn/2+1) * sizeof(double) ) ) == NULL ) {
  xs_write("kynrefrev: Failed to allocate memory for tmp arrays.", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  goto error;
}

if( ( r_part   = (double *) malloc( ne * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( im_part  = (double *) malloc( ne * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( ampl     = (double *) malloc( ne * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( phase    = (double *) malloc( ne * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( phase_u1 = (double *) malloc( ne * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( phase_u2 = (double *) malloc( ne * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( ampl_etot2   = (double *) malloc( ne * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( phase_etot2  = (double *) malloc( ne * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( r_part_tot   = (double *) malloc( ne * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( im_part_tot  = (double *) malloc( ne * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( ampl_tot     = (double *) malloc( ne * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( phase_tot    = (double *) malloc( ne * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( phase_tot_u1 = (double *) malloc( ne * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( phase_tot_u2 = (double *) malloc( ne * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( ampl_tot_etot2  = (double *) malloc( ne * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( phase_tot_etot2 = (double *) malloc( ne * (nn/2+1) * sizeof(double) ) ) == NULL ) {
  xs_write("kynrefrev: Failed to allocate memory for tmp arrays.", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  goto error;
}

if( nbands > 0 && (
    ( r_part_bands   = (double *) malloc( nbands * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( im_part_bands  = (double *) malloc( nbands * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( ampl_bands     = (double *) malloc( nbands * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( phase_bands    = (double *) malloc( nbands * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( phase_bands_u1 = (double *) malloc( nbands * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( r_part_bands_tot   = (double *) malloc( nbands * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( im_part_bands_tot  = (double *) malloc( nbands * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( ampl_bands_tot     = (double *) malloc( nbands * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( phase_bands_tot    = (double *) malloc( nbands * (nn/2+1) * sizeof(double) ) ) == NULL ||
    ( phase_bands_tot_u1 = (double *) malloc( nbands * (nn/2+1) * sizeof(double) ) ) == NULL ) ){
  xs_write("kynrefrev: Failed to allocate memory for tmp arrays.", 5);
  for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
  goto error;
}

// Let's fourier transform
for(ifr=0;ifr<=nn/2;ifr++) freq[ifr] = ifr / tot_time_sec;
// energy dependent Fourier transform
if(abs(photar_sw) <= 10){
  for(ie=0;ie<ne;ie++){
    for(it=0;it<nn;it++){
     data_r[it] = far[ie+ne*it] / far_prim[ie] / flare_duration_sec;
     data_im[it] = 0.;
    }
    four(data_r, data_im, nn, ISGN);
    k1=0;
    k2=0;
    for(ifr=0;ifr<=nn/2;ifr++){
      data_r[ifr]=data_r[ifr]*tot_time_sec;
      data_im[ifr]=data_im[ifr]*tot_time_sec;
// Let's correct for the Box shape:
      Im_sinc = PI * flare_duration_sec * freq[ifr];
      if(Im_sinc==0.) Re_sinc = 1.;
      else Re_sinc = Im_sinc / tan(Im_sinc);
      r_part[ie+ne*ifr] = data_r[ifr] * Re_sinc -
                          data_im[ifr] * Im_sinc;
      im_part[ie+ne*ifr] = data_im[ifr] * Re_sinc + 
                           data_r[ifr] * Im_sinc;
      ampl[ie+ne*ifr] = sqrt(r_part[ie+ne*ifr]*r_part[ie+ne*ifr] +
                             im_part[ie+ne*ifr]*im_part[ie+ne*ifr]);
      phase[ie+ne*ifr] = atan2(im_part[ie+ne*ifr], r_part[ie+ne*ifr]);
// Let's correct for the shift in time
      phase[ie+ne*ifr] -= PI2 * freq[ifr] * t0;
      r_part[ie+ne*ifr] = ampl[ie+ne*ifr] * cos(phase[ie+ne*ifr]);
      im_part[ie+ne*ifr] = ampl[ie+ne*ifr] * sin(phase[ie+ne*ifr]);
      phase[ie+ne*ifr] = atan2(im_part[ie+ne*ifr], r_part[ie+ne*ifr]);
      phase_u1[ie+ne*ifr] = atan2(im_part[ie+ne*ifr], r_part[ie+ne*ifr]) + k1*PI2;
      if(ifr>0 && (phase_u1[ie+ne*ifr]-phase_u1[ie+ne*(ifr-1)])<-PI)k1++;
      if(ifr>0 && (phase_u1[ie+ne*ifr]-phase_u1[ie+ne*(ifr-1)])>PI)k1--;
      phase_u1[ie+ne*ifr] = atan2(im_part[ie+ne*ifr], r_part[ie+ne*ifr]) + k1*PI2;
      phase_u2[ie+ne*ifr] = atan2(im_part[ie+ne*ifr], r_part[ie+ne*ifr]) + k2*PI2;
      if(ie>0 && (phase_u2[ie+ne*ifr]-phase_u2[ie-1+ne*ifr])<-PI)k2++;
      if(ie>0 && (phase_u2[ie+ne*ifr]-phase_u2[ie-1+ne*ifr])>PI)k2--;
      phase_u2[ie+ne*ifr] = atan2(im_part[ie+ne*ifr], r_part[ie+ne*ifr]) + k2*PI2;
    }
  }
//compute frequency wrapping for XSPEC
  ifr = 1;
  nonwrap = 1;
  while( ifr<nn/2 && nonwrap ){
    ifr++;
    for(ie=0;ie<ne;ie++)
      nonwrap *= (im_part[ie+ne*ifr]*im_part[ie+ne*(ifr-1)] > 0.);
  }
  fwrap = freq[ifr-1];
// energy dependent Fourier transform when primary is included
#ifndef OUTSIDE_XSPEC
  if( photar_sw > 0){
#endif
    for(ie=0;ie<ne;ie++){
      k1=0;
      k2=0;
      for(ifr=0;ifr<=nn/2;ifr++){
        r_part_tot[ie+ne*ifr] = far_prim[ie] * (1.+r_part[ie+ne*ifr]);
        im_part_tot[ie+ne*ifr] = far_prim[ie] * im_part[ie+ne*ifr];
        ampl_tot[ie+ne*ifr] = far_prim[ie] *
                              sqrt((1.+r_part[ie+ne*ifr])*(1.+r_part[ie+ne*ifr]) +
                                   im_part[ie+ne*ifr]*im_part[ie+ne*ifr]);
        phase_tot[ie+ne*ifr] = atan2(im_part[ie+ne*ifr],
                                     1.+r_part[ie+ne*ifr]);
        phase_tot_u1[ie+ne*ifr] = atan2(im_part[ie+ne*ifr],
                                  1.+r_part[ie+ne*ifr])+k1*PI2;
        if(ifr>0 && (phase_tot_u1[ie+ne*ifr]-phase_tot_u1[ie+ne*(ifr-1)])<-PI)k1++;
        if(ifr>0 && (phase_tot_u1[ie+ne*ifr]-phase_tot_u1[ie+ne*(ifr-1)])>PI)k1--;
        phase_tot_u1[ie+ne*ifr] = atan2(im_part[ie+ne*ifr],
                                  1.+r_part[ie+ne*ifr])+k1*PI2;
        phase_tot_u2[ie+ne*ifr] = atan2(im_part[ie+ne*ifr],
                                  1.+r_part[ie+ne*ifr])+k2*PI2;
        if(ie>0 && (phase_tot_u2[ie+ne*ifr]-phase_tot_u2[ie-1+ne*ifr])<-PI)k2++;
        if(ie>0 && (phase_tot_u2[ie+ne*ifr]-phase_tot_u2[ie-1+ne*ifr])>PI)k2--;
        phase_tot_u2[ie+ne*ifr] = atan2(im_part[ie+ne*ifr],
                                  1.+r_part[ie+ne*ifr])+k2*PI2;
      }
    }
#ifndef OUTSIDE_XSPEC
  }
#endif
}
//Fourier transform for chosen energy bands
//we do not need the following condition since for photar_sw<=10 we also need
//to compute the FFT of the reference band if nbands>0!!!
//#ifndef OUTSIDE_XSPEC
//if(abs(photar_sw) > 10){
//#endif
for(iband=0;iband<nbands;iband++){
  for(it=0;it<nn;it++){
    if(it<nn) data_r[it] = flux_bands[iband+nbands*it] / 
                           flux_bands_prim[iband] / flare_duration_sec;
    else data_r[it] = 0.;
    data_im[it] = 0.;
  }
  four(data_r, data_im, nn, ISGN);
  k1=0;
//  k2=0;
  for(ifr=0;ifr<=nn/2;ifr++){
    data_r[ifr] = data_r[ifr]*tot_time_sec;
//    r_part_bands[iband+nbands*ifr] = data_r[ifr];
    data_im[ifr] = data_im[ifr]*tot_time_sec;
//    im_part_bands[iband+nbands*ifr] = data_im[ifr];
// Let's correct for the Box shape:
    Im_sinc = PI * flare_duration_sec * freq[ifr];
    if(Im_sinc==0.) Re_sinc = 1.;
    else Re_sinc = Im_sinc / tan(Im_sinc);
    r_part_bands[iband+nbands*ifr] = 
      data_r[ifr] * Re_sinc - data_im[ifr] * Im_sinc;
    im_part_bands[iband+nbands*ifr] = 
      data_im[ifr] * Re_sinc + data_r[ifr] * Im_sinc;
    ampl_bands[iband+nbands*ifr] =
      sqrt(r_part_bands[iband+nbands*ifr]*r_part_bands[iband+nbands*ifr]+
           im_part_bands[iband+nbands*ifr]*im_part_bands[iband+nbands*ifr]);
    phase_bands[iband+nbands*ifr] = 
      atan2(im_part_bands[iband+nbands*ifr],r_part_bands[iband+nbands*ifr]);
// Let's correct for the shift in time
    phase_bands[iband+nbands*ifr] -= PI2 * freq[ifr] * t0;
    r_part_bands[iband+nbands*ifr] = ampl_bands[iband+nbands*ifr] * 
                                     cos(phase_bands[iband+nbands*ifr]);
    im_part_bands[iband+nbands*ifr] = ampl_bands[iband+nbands*ifr] * 
                                      sin(phase_bands[iband+nbands*ifr]);
    phase_bands[iband+nbands*ifr] = 
      atan2(im_part_bands[iband+nbands*ifr],r_part_bands[iband+nbands*ifr]);
    phase_bands_u1[iband+nbands*ifr] = atan2(im_part_bands[iband+nbands*ifr], r_part_bands[iband+nbands*ifr])+k1*PI2;
    if(ifr>0 && (phase_bands_u1[iband+nbands*ifr]-phase_bands_u1[iband+nbands*(ifr-1)])<-PI)k1++;
    if(ifr>0 && (phase_bands_u1[iband+nbands*ifr]-phase_bands_u1[iband+nbands*(ifr-1)])>PI)k1--;
    phase_bands_u1[iband+nbands*ifr] = atan2(im_part_bands[iband+nbands*ifr], r_part_bands[iband+nbands*ifr])+k1*PI2;
//    phase_bands_u2[iband+nbands*ifr] = atan2(im_part_bands[iband+nbands*ifr], r_part_bands[iband+nbands*ifr])+k2*PI2;
//    if(iband>0 && (phase_bands_u2[iband+nbands*ifr]-phase_bands_u2[iband-1+nbands*ifr])<-PI)k2++;
//    if(iband>0 && (phase_bands_u2[iband+nbands*ifr]-phase_bands_u2[iband-1+nbands*ifr])>PI)k2--;
//    phase_bands_u2[iband+nbands*ifr] = atan2(im_part_bands[iband+nbands*ifr], r_part_bands[iband+nbands*ifr])+k2*PI2;
  }
}
// compute amplitudes and phases for the whole energy band excluding current energy bin
if( nbands>0 && exclude_energy && ( photar_sw == -7 || photar_sw == -8 ) )
  for(ie=0;ie<ne;ie++)
    for(ifr=0;ifr<=nn/2;ifr++){
      ttmp = flux_bands_prim[nbands-1] * r_part_bands[nbands-1+nbands*ifr] - 
             far_prim[ie] * r_part[ie+ne*ifr];
      ttmp1 = flux_bands_prim[nbands-1] * im_part_bands[nbands-1+nbands*ifr] - 
              far_prim[ie] * im_part[ie+ne*ifr];
      if(photar_sw == -7) 
        ampl_etot2[ie+ne*ifr] = sqrt( ttmp*ttmp + ttmp1*ttmp1 ) / 
                                ( flux_bands_prim[nbands-1] - far_prim[ie] );
      if(photar_sw == -8) phase_etot2[ie+ne*ifr] = atan2( ttmp1, ttmp );
    }
//compute frequency wrapping for XSPEC
if(abs(photar_sw) > 10){
  ifr = 1;
  nonwrap = 1;
  while( ifr<nn/2 && nonwrap ){
    ifr++;
    for(iband=0;iband<nbands;iband++)
      nonwrap *= (im_part_bands[iband+nbands*ifr]*im_part_bands[iband+nbands*(ifr-1)] > 0.);
  }
  fwrap = freq[ifr-1];
}
// Fourier transform for chosen energy bands when primary is included
#ifndef OUTSIDE_XSPEC
if( photar_sw > 0){
#endif
  for(iband=0;iband<nbands;iband++){
    k1=0;
//    k2=0;
    for(ifr=0;ifr<=nn/2;ifr++){
      r_part_bands_tot[iband+nbands*ifr] = flux_bands_prim[iband] *
                                           (1.+r_part_bands[iband+nbands*ifr]);
      im_part_bands_tot[iband+nbands*ifr] = flux_bands_prim[iband] *
                                            im_part_bands[iband+nbands*ifr];
      ampl_bands_tot[iband+nbands*ifr] = flux_bands_prim[iband]*
        sqrt((1.+r_part_bands[iband+nbands*ifr]) *
             (1.+r_part_bands[iband+nbands*ifr]) +
             im_part_bands[iband+nbands*ifr]*im_part_bands[iband+nbands*ifr]);
      phase_bands_tot[iband+nbands*ifr]=
        atan2(im_part_bands[iband+nbands*ifr], 1.+r_part_bands[iband+nbands*ifr]);
      phase_bands_tot_u1[iband+nbands*ifr]=
        atan2(im_part_bands[iband+nbands*ifr], 1.+r_part_bands[iband+nbands*ifr])+
              k1*PI2;
      if(ifr>0 && (phase_bands_tot_u1[iband+nbands*ifr]-phase_bands_tot_u1[iband+nbands*(ifr-1)])<-PI)k1++;
      if(ifr>0 && (phase_bands_tot_u1[iband+nbands*ifr]-phase_bands_tot_u1[iband+nbands*(ifr-1)])>PI)k1--;
      phase_bands_tot_u1[iband+nbands*ifr]=atan2(im_part_bands[iband+nbands*ifr],
        1.+r_part_bands[iband+nbands*ifr])+k1*PI2;
//      phase_bands_tot_u2[iband+nbands*ifr]=atan2(im_part_bands[iband+nbands*ifr],
//        1.+r_part_bands[iband+nbands*ifr])+k2*PI2;
//      if(iband>0 && (phase_bands_tot_u2[iband+nbands*ifr]-phase_bands_tot_u2[iband-1+nbands*ifr])<-PI)k2++;
//      if(iband>0 && (phase_bands_tot_u2[iband+nbands*ifr]-phase_bands_tot_u2[iband-1+nbands*ifr])>PI)k2--;
//      phase_bands_tot_u2[iband+nbands*ifr]=atan2(im_part_bands[iband+nbands*ifr],
//        1.+r_part_bands[iband+nbands*ifr])+k2*PI2;
    }
  }
// compute amplitudes and phases for the whole energy band excluding current energy bin
  if( nbands>0 && exclude_energy && ( photar_sw == 7 || photar_sw == 8 ) )
    for(ie=0;ie<ne;ie++)
      for(ifr=0;ifr<=nn/2;ifr++){
        ttmp = r_part_bands_tot[nbands-1+nbands*ifr] - r_part_tot[ie+ne*ifr];
        ttmp1 = im_part_bands_tot[nbands-1+nbands*ifr] - im_part_tot[ie+ne*ifr];
        if(photar_sw == 7) ampl_tot_etot2[ie+ne*ifr] = sqrt( ttmp*ttmp + ttmp1*ttmp1 );
        if(photar_sw == 8) phase_tot_etot2[ie+ne*ifr] = atan2( ttmp1, ttmp );
      }
#ifndef OUTSIDE_XSPEC
}
//}
#endif
// integrated Fourier transform up to the first wrapping zero
if( abs(photar_sw) <= 10 ){
#ifndef OUTSIDE_XSPEC
  if( frequency1 == 0. ){
#endif
//just reflection
    if( photar_sw < 0){
      for(ie=0;ie<ne;ie++){
        r_part_int[ie] = r_part[ie];
        im_part_int[ie] = im_part[ie];
        ampl2_int[ie] = 0;
        delay2_int[ie] = 0;
      }
      ifr = 2;
      nonwrap = 1;
      for(ie=0;ie<ne;ie++) 
//        nonwrap *= !(r_part[ie+ne*ifr] > 0. && 
//                     im_part[ie+ne*ifr]*im_part[ie+ne*(ifr-1)] < 0.);
        nonwrap *= (im_part[ie+ne*ifr]*im_part[ie+ne*(ifr-1)] > 0.);
      while( ifr<=nn/2 && nonwrap ){
        for(ie=0;ie<ne;ie++){
          r_part_int[ie] += r_part[ie+ne*ifr];
          im_part_int[ie] += im_part[ie+ne*ifr];
          ampl2_int[ie] += ampl[ie+ne*ifr];
          delay2_int[ie] += phase[ie+ne*ifr] / PI2 / freq[ifr];
        }
        ifr++;
        if(ifr<=nn/2){
          nonwrap = 1;
//          for(ie=0;ie<ne;ie++) nonwrap *= !(r_part[ie+ne*ifr] > 0. && 
//            im_part[ie+ne*ifr]*im_part[ie+ne*(ifr-1)] < 0.);
          for(ie=0;ie<ne;ie++) 
            nonwrap *= (im_part[ie+ne*ifr]*im_part[ie+ne*(ifr-1)] > 0.);
        }
      }
      ifr--;
      freq_wrap = freq[ifr];
      for(ie=0;ie<ne;ie++){
        r_part_int[ie] /= ( ifr + 1 );
        im_part_int[ie] /= ( ifr + 1 );
        ampl2_int[ie] /= ( ifr + 1 );
        delay2_int[ie] /= ( ifr + 1 );
      }
      for(ie=0;ie<ne;ie++){
        k2=0;
        ampl_int[ie] = sqrt( r_part_int[ie]  * r_part_int[ie] +
                             im_part_int[ie] * im_part_int[ie] );
        phase_int[ie] = atan2(im_part_int[ie], r_part_int[ie]);
        phase_int_u2[ie] = atan2(im_part_int[ie], r_part_int[ie])+k2*PI2;
        if(ie>0 && (phase_int_u2[ie]-phase_int_u2[ie-1]) < -PI) k2++;
        if(ie>0 && (phase_int_u2[ie]-phase_int_u2[ie-1]) > PI) k2--;
        phase_int_u2[ie] = atan2(im_part_int[ie], r_part_int[ie])+k2*PI2;
      }
//    reference energy band Fourier transform integrated up to the first wrapping 
//    zero
      if(nbands>0){
        r_part_int_etot = r_part_bands[nbands-1];
        im_part_int_etot = im_part_bands[nbands-1];
        ampl2_int_etot = 0;
        delay2_int_etot = 0;
        ifr0=ifr;
        for( ifr=1;ifr<=ifr0;ifr++ ){
          r_part_int_etot += r_part_bands[nbands-1+nbands*ifr];
          im_part_int_etot += im_part_bands[nbands-1+nbands*ifr];
          ampl2_int_etot += ampl_bands[nbands-1+nbands*ifr];
          delay2_int_etot += phase_bands[nbands-1+nbands*ifr] / PI2 / freq[ifr];
        }
        r_part_int_etot /= ( ifr0+1 );
        im_part_int_etot /= ( ifr0+1 );
        ampl2_int_etot /= ( ifr0+1 );
        delay2_int_etot /= ( ifr0+1 );
        ampl_int_etot = sqrt( r_part_int_etot  * r_part_int_etot +
                              im_part_int_etot * im_part_int_etot );
        phase_int_etot = atan2(im_part_int_etot, r_part_int_etot);
        if( exclude_energy && ( photar_sw == -7 || photar_sw == -8 ) ){
          for(ie=0;ie<ne;ie++){
            ampl2_int_etot2[ie] = 0;
            delay2_int_etot2[ie] = 0;
            ifr0=ifr;
            for( ifr=1;ifr<=ifr0;ifr++ ){
              ampl2_int_etot2[ie] += ampl_etot2[ie+ne*ifr];
              delay2_int_etot2[ie] += phase_etot2[ie+ne*ifr] / PI2 / freq[ifr];
            }
            ampl2_int_etot2[ie] /= ( ifr0+1 );
            delay2_int_etot2[ie] /= ( ifr0+1 );
          }
        }
      }else{
        ampl2_int_etot = 1.;
        delay2_int_etot = 0.;
        ampl_int_etot = 1.;
        phase_int_etot = 0.;
        if( photar_sw == -7 || photar_sw == -8 )
          for(ie=0;ie<ne;ie++){
            ampl2_int_etot2[ie] = 1.;
            delay2_int_etot2[ie] = 0.;
          }
      }
    }
//including primary  
#ifndef OUTSIDE_XSPEC
    if( photar_sw > 0){
#endif
      for(ie=0;ie<ne;ie++){
        r_part_tot_int[ie] = r_part_tot[ie];
        im_part_tot_int[ie] = im_part_tot[ie];
        ampl2_tot_int[ie] = 0;
        delay2_tot_int[ie] = 0;
      }
      ifr = 2;
      nonwrap = 1;
      for(ie=0;ie<ne;ie++) 
//        nonwrap *= !(r_part_tot[ie+ne*ifr] > 0. && 
//                     im_part_tot[ie+ne*ifr]*im_part_tot[ie+ne*(ifr-1)] < 0.);
        nonwrap *= (im_part_tot[ie+ne*ifr]*im_part_tot[ie+ne*(ifr-1)] > 0.);
      while( ifr<=nn/2 && nonwrap ){
        for(ie=0;ie<ne;ie++){
          r_part_tot_int[ie] += r_part_tot[ie+ne*ifr];
          im_part_tot_int[ie] += im_part_tot[ie+ne*ifr];
          ampl2_tot_int[ie] += ampl_tot[ie+ne*ifr];
          delay2_tot_int[ie] += phase_tot[ie+ne*ifr] / PI2 / freq[ifr];
        }
        ifr++;
        if(ifr<=nn/2){
          nonwrap = 1;
//          for(ie=0;ie<ne;ie++) nonwrap *= !(r_part_tot[ie+ne*ifr] > 0. && 
//            im_part_tot[ie+ne*ifr]*im_part_tot[ie+ne*(ifr-1)] < 0.);
          for(ie=0;ie<ne;ie++) 
            nonwrap *= (im_part_tot[ie+ne*ifr]*im_part_tot[ie+ne*(ifr-1)] > 0.);
        }
      }
      ifr--;
      freq_wrap = freq[ifr];
      freq_wrap0 = freq_wrap;
      for(ie=0;ie<ne;ie++){
        r_part_tot_int[ie] /= ( ifr+1 );
        im_part_tot_int[ie] /= ( ifr+1 );
        ampl2_tot_int[ie] /= ( ifr+1 );
        delay2_tot_int[ie] /= ( ifr+1 );
      }
      for(ie=0;ie<ne;ie++){
        k2=0;
        ampl_tot_int[ie] = sqrt( r_part_tot_int[ie]  * r_part_tot_int[ie] +
                                 im_part_tot_int[ie] * im_part_tot_int[ie] );
        phase_tot_int[ie] = atan2(im_part_tot_int[ie], r_part_tot_int[ie]);
        phase_tot_int_u2[ie] = atan2(im_part_tot_int[ie], r_part_tot_int[ie])+k2*PI2;
        if(ie>0 && (phase_tot_int_u2[ie]-phase_tot_int_u2[ie-1]) < -PI) k2++;
        if(ie>0 && (phase_tot_int_u2[ie]-phase_tot_int_u2[ie-1]) > PI) k2--;
        phase_tot_int_u2[ie] = atan2(im_part_tot_int[ie], r_part_tot_int[ie])+k2*PI2;
      }
//    reference energy band Fourier transform integrated up to the first wrapping 
//    zero
      if(nbands>0){
        r_part_tot_int_etot = r_part_bands_tot[nbands-1];
        im_part_tot_int_etot = im_part_bands_tot[nbands-1];
        ampl2_tot_int_etot = 0;
        delay2_tot_int_etot = 0;
        ifr0=ifr;
        for( ifr=1;ifr<=ifr0;ifr++){
          r_part_tot_int_etot += r_part_bands_tot[nbands-1+nbands*ifr];
          im_part_tot_int_etot += im_part_bands_tot[nbands-1+nbands*ifr];
          ampl2_tot_int_etot += ampl_bands_tot[nbands-1+nbands*ifr];
          delay2_tot_int_etot += phase_bands_tot[nbands-1+nbands*ifr] / PI2 / freq[ifr];
        }
        r_part_tot_int_etot /= ( ifr0+1 );
        im_part_tot_int_etot /= ( ifr0+1 );
        ampl2_tot_int_etot /= ( ifr0+1 );
        delay2_tot_int_etot /= ( ifr0+1 );
        ampl_tot_int_etot = sqrt( r_part_tot_int_etot  * r_part_tot_int_etot +
                                  im_part_tot_int_etot * im_part_tot_int_etot );
        phase_tot_int_etot = atan2(im_part_tot_int_etot, r_part_tot_int_etot);
        if( exclude_energy && ( photar_sw == 7 || photar_sw == 8 ) )
          for(ie=0;ie<ne;ie++){
            ampl2_tot_int_etot2[ie] = 0;
            delay2_tot_int_etot2[ie] = 0;
            ifr0=ifr;
            for( ifr=1;ifr<=ifr0;ifr++ ){
              ampl2_tot_int_etot2[ie] += ampl_tot_etot2[ie+ne*ifr];
              delay2_tot_int_etot2[ie] += phase_tot_etot2[ie+ne*ifr] / PI2 / freq[ifr];
            }
            ampl2_tot_int_etot2[ie] /= ( ifr0+1 );
            delay2_tot_int_etot2[ie] /= ( ifr0+1 );
          }
      }else{
        ampl2_tot_int_etot = 1.;
        delay2_tot_int_etot = 0.;
        ampl_tot_int_etot = 1.;
        phase_tot_int_etot = 0.;
        if( photar_sw == 7 || photar_sw == 8 ){
          for(ie=0;ie<ne;ie++){
            ampl2_tot_int_etot2[ie] = 1.;
            delay2_tot_int_etot2[ie] = 0.;
          }
        }
      }
#ifndef OUTSIDE_XSPEC
    }
  }
#endif
}
// frequency band integrated energy dependend Fourier transform
if( abs(photar_sw) <= 10 && frequency1 > 0 && frequency1 < frequency2 ){
//given frequency1 and frequency2, find the corresponding index in freq[]:
  ifr0 = (int) ceil( frequency1 / freq[1] );
  if(ifr0 < 2) ifr0 = 2;
  else if(ifr0 > nn/2) ifr0 = nn/2;
  ifrn = (int) ceil(frequency2 / freq[1]);
  if(ifrn < 2) ifrn = 2;
  else if(ifrn > nn/2) ifrn = nn/2;
//the first and last partial bin
  ttmp = (frequency1 - freq[ifr0 - 1]) / (freq[ifr0] - freq[ifr0 - 1]);
  ttmp1 = 1. - ttmp;
  utmp = (frequency2 - freq[ifrn - 1]) / (freq[ifrn] - freq[ifrn - 1]);
  utmp1 = 1. - utmp;
//just reflection
  if( photar_sw < 0){
    for(ie=0;ie<ne;ie++){
//  real part
      y1 = r_part[ie+ne*(ifr0-1)];
      y2 = r_part[ie+ne*ifr0];
      y3 = r_part[ie+ne*(ifrn-1)];
      y4 = r_part[ie+ne*ifrn];
      if( ifr0 < ifrn ){
        r_part_fband[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                             r_part[ie+ne*ifr0] ) / 2. *
                           ( freq[ifr0] - frequency1 );
        r_part_fband[ie] += ( r_part[ie+ne*(ifrn-1)] + 
                              (utmp1 * y3 + utmp * y4) ) / 2. *
                            ( frequency2 - freq[ifrn-1] );
      }else if( ifr0 == ifrn ){
        r_part_fband[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                             (utmp1 * y3 + utmp * y4) ) / 2. *
                           ( frequency2 - frequency1 );
      }
//    ...all the whole bins
      for(ifr=ifr0+1;ifr<ifrn;ifr++) 
        r_part_fband[ie] += ( r_part[ie+ne*ifr] + 
                              r_part[ie+ne*(ifr-1)] ) / 2. * freq[1];
//  imaginary part
      y1 = im_part[ie+ne*(ifr0-1)];
      y2 = im_part[ie+ne*ifr0];
      y3 = im_part[ie+ne*(ifrn-1)];
      y4 = im_part[ie+ne*ifrn];
      if( ifr0 < ifrn ){
        im_part_fband[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                              im_part[ie+ne*ifr0] ) / 2. *
                            ( freq[ifr0] - frequency1 );
        im_part_fband[ie] += ( im_part[ie+ne*(ifrn-1)] + 
                               (utmp1 * y3 + utmp * y4) ) / 2. *
                             ( frequency2 - freq[ifrn-1] );
      }else if( ifr0 == ifrn ){
        im_part_fband[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                              (utmp1 * y3 + utmp * y4) ) / 2. *
                            ( frequency2 - frequency1 );
      }
//    ...all the whole bins
      for(ifr=ifr0+1;ifr<ifrn;ifr++) 
        im_part_fband[ie] += ( im_part[ie+ne*ifr] + 
                               im_part[ie+ne*(ifr-1)]) / 2. * freq[1];
//  amplitude (averaged over frequency, i.e. not from averaged Re and Im parts!)
      y1 = ampl[ie+ne*(ifr0-1)];
      y2 = ampl[ie+ne*ifr0];
      y3 = ampl[ie+ne*(ifrn-1)];
      y4 = ampl[ie+ne*ifrn];
      if( ifr0 < ifrn ){
        ampl2_fband[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                              ampl[ie+ne*ifr0] ) / 2. *
                            ( freq[ifr0] - frequency1 );
        ampl2_fband[ie] += ( ampl[ie+ne*(ifrn-1)] + 
                               (utmp1 * y3 + utmp * y4) ) / 2. *
                             ( frequency2 - freq[ifrn-1] );
      }else if( ifr0 == ifrn ){
        ampl2_fband[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                              (utmp1 * y3 + utmp * y4) ) / 2. *
                            ( frequency2 - frequency1 );
      }
//    ...all the whole bins
      for(ifr=ifr0+1;ifr<ifrn;ifr++) 
        ampl2_fband[ie] += ( ampl[ie+ne*ifr] + 
                               ampl[ie+ne*(ifr-1)]) / 2. * freq[1];
//  delay (averaged over frequency, i.e. not from averaged Re and Im parts!)
      y1 = phase[ie+ne*(ifr0-1)] / freq[ifr0-1];
      y2 = phase[ie+ne*ifr0] / freq[ifr0];
      y3 = phase[ie+ne*(ifrn-1)] / freq[ifrn-1];
      y4 = phase[ie+ne*ifrn] / freq[ifrn];
      if( ifr0 < ifrn ){
        delay2_fband[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                            phase[ie+ne*ifr0] / freq[ifr0] ) / 2. *
                          ( freq[ifr0] - frequency1 );
        delay2_fband[ie] += ( phase[ie+ne*(ifrn-1)] / freq[ifrn-1] + 
                             (utmp1 * y3 + utmp * y4) ) / 2. *
                           ( frequency2 - freq[ifrn-1] );
      }else if( ifr0 == ifrn ){
        delay2_fband[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                            (utmp1 * y3 + utmp * y4) ) / 2. *
                          ( frequency2 - frequency1 );
      }
//    ...all the whole bins
      for(ifr=ifr0+1;ifr<ifrn;ifr++) 
        delay2_fband[ie] += ( phase[ie+ne*ifr] / freq[ifr] + 
                             phase[ie+ne*(ifr-1)] / freq[ifr-1] ) / 2. 
                             * freq[1];
    }
    for(ie=0;ie<ne;ie++){
      r_part_fband[ie] /= (frequency2 - frequency1);
      im_part_fband[ie] /= (frequency2 - frequency1);
      ampl2_fband[ie] /= (frequency2 - frequency1);
      delay2_fband[ie] /= ( (frequency2 - frequency1) * PI2 );
    }
    for(ie=0;ie<ne;ie++){
      k2=0;
      ampl_fband[ie] = sqrt( r_part_fband[ie]  * r_part_fband[ie] +
                             im_part_fband[ie] * im_part_fband[ie] );
      phase_fband[ie] = atan2(im_part_fband[ie], r_part_fband[ie]);
      phase_fband_u2[ie] = atan2(im_part_fband[ie], r_part_fband[ie])+k2*PI2;
      if(ie>0 && (phase_fband_u2[ie]-phase_fband_u2[ie-1]) < -PI) k2++;
      if(ie>0 && (phase_fband_u2[ie]-phase_fband_u2[ie-1]) > PI) k2--;
      phase_fband_u2[ie] = atan2(im_part_fband[ie], r_part_fband[ie])+k2*PI2;
    }
//reference energy band Fourier transform integrated in the given frequency band
    if(nbands>0){
//    real part
      y1 = r_part_bands[nbands-1+nbands*(ifr0-1)];
      y2 = r_part_bands[nbands-1+nbands*ifr0];
      y3 = r_part_bands[nbands-1+nbands*(ifrn-1)];
      y4 = r_part_bands[nbands-1+nbands*ifrn];
      if( ifr0 < ifrn ){
        r_part_fband_etot = ( (ttmp1 * y1 + ttmp * y2) + 
                              r_part_bands[nbands-1+nbands*ifr0] ) / 2. *
                            ( freq[ifr0] - frequency1 );
        r_part_fband_etot += ( r_part_bands[nbands-1+nbands*(ifrn-1)] + 
                               (utmp1 * y3 + utmp * y4) ) / 2. *
                             ( frequency2 - freq[ifrn-1] );
      }else if( ifr0 == ifrn ){
        r_part_fband_etot = ( (ttmp1 * y1 + ttmp * y2) + 
                              (utmp1 * y3 + utmp * y4) ) / 2. *
                            ( frequency2 - frequency1 );
      }
//    ...all the whole bins
      for(ifr=ifr0+1;ifr<ifrn;ifr++) 
        r_part_fband_etot += ( r_part_bands[nbands-1+nbands*ifr] + 
                               r_part_bands[nbands-1+nbands*(ifr-1)]) / 2. * freq[1];
//    imaginary part
      y1 = im_part_bands[nbands-1+nbands*(ifr0-1)];
      y2 = im_part_bands[nbands-1+nbands*ifr0];
      y3 = im_part_bands[nbands-1+nbands*(ifrn-1)];
      y4 = im_part_bands[nbands-1+nbands*ifrn];
      if( ifr0 < ifrn ){
        im_part_fband_etot = ( (ttmp1 * y1 + ttmp * y2) + 
                               im_part_bands[nbands-1+nbands*ifr0] ) / 2. *
                             ( freq[ifr0] - frequency1 );
        im_part_fband_etot += ( im_part_bands[nbands-1+nbands*(ifrn-1)] + 
                                (utmp1 * y3 + utmp * y4) ) / 2. *
                              ( frequency2 - freq[ifrn-1] );
      }else if( ifr0 == ifrn ){
        im_part_fband_etot = ( (ttmp1 * y1 + ttmp * y2) + 
                               (utmp1 * y3 + utmp * y4) ) / 2. *
                             ( frequency2 - frequency1 );
      }
//    ...all the whole bins
      for(ifr=ifr0+1;ifr<ifrn;ifr++) 
        im_part_fband_etot += ( im_part_bands[nbands-1+nbands*ifr] + 
                                im_part_bands[nbands-1+nbands*(ifr-1)]) / 2. * freq[1];
//    amplitude (averaged over frequency, i.e. not from averaged Re and Im parts!)
      y1 = ampl_bands[nbands-1+nbands*(ifr0-1)];
      y2 = ampl_bands[nbands-1+nbands*ifr0];
      y3 = ampl_bands[nbands-1+nbands*(ifrn-1)];
      y4 = ampl_bands[nbands-1+nbands*ifrn];
      if( ifr0 < ifrn ){
        ampl2_fband_etot = ( (ttmp1 * y1 + ttmp * y2) + 
                               ampl_bands[nbands-1+nbands*ifr0] ) / 2. *
                             ( freq[ifr0] - frequency1 );
        ampl2_fband_etot += ( ampl_bands[nbands-1+nbands*(ifrn-1)] + 
                                (utmp1 * y3 + utmp * y4) ) / 2. *
                              ( frequency2 - freq[ifrn-1] );
      }else if( ifr0 == ifrn ){
        ampl2_fband_etot = ( (ttmp1 * y1 + ttmp * y2) + 
                               (utmp1 * y3 + utmp * y4) ) / 2. *
                             ( frequency2 - frequency1 );
      }
//    ...all the whole bins
      for(ifr=ifr0+1;ifr<ifrn;ifr++) 
        ampl2_fband_etot += ( ampl_bands[nbands-1+nbands*ifr] + 
                              ampl_bands[nbands-1+nbands*(ifr-1)]) / 2. * freq[1];
//    delay (averaged over frequency, i.e. not from averaged Re and Im parts!)
      y1 = phase_bands[nbands-1+nbands*(ifr0-1)] / freq[ifr0-1];
      y2 = phase_bands[nbands-1+nbands*ifr0] / freq[ifr0];
      y3 = phase_bands[nbands-1+nbands*(ifrn-1)] / freq[ifrn-1];
      y4 = phase_bands[nbands-1+nbands*ifrn] / freq[ifrn];
      if( ifr0 < ifrn ){
        delay2_fband_etot = ( (ttmp1 * y1 + ttmp * y2) + 
                             phase_bands[nbands-1+nbands*ifr0] / freq[ifr0] ) / 2. *
                           ( freq[ifr0] - frequency1 );
        delay2_fband_etot += ( phase_bands[nbands-1+nbands*(ifrn-1)] / freq[ifrn-1] + 
                              (utmp1 * y3 + utmp * y4) ) / 2. *
                            ( frequency2 - freq[ifrn-1] );
      }else if( ifr0 == ifrn ){
        delay2_fband_etot = ( (ttmp1 * y1 + ttmp * y2) + 
                             (utmp1 * y3 + utmp * y4) ) / 2. *
                           ( frequency2 - frequency1 );
      }
//    ...all the whole bins
      for(ifr=ifr0+1;ifr<ifrn;ifr++)
        delay2_fband_etot += ( phase_bands[nbands-1+nbands*ifr] / freq[ifr] + 
                              phase_bands[nbands-1+nbands*(ifr-1)] / freq[ifr-1] ) / 2. 
                              * freq[1];
      r_part_fband_etot /= (frequency2 - frequency1);
      im_part_fband_etot /= (frequency2 - frequency1);
      ampl2_fband_etot /= (frequency2 - frequency1);
      delay2_fband_etot /= ( (frequency2 - frequency1) * PI2 );
      ampl_fband_etot  = sqrt( r_part_fband_etot  * r_part_fband_etot +
                               im_part_fband_etot * im_part_fband_etot );
      phase_fband_etot = atan2(im_part_fband_etot, r_part_fband_etot);
//    amplitude (averaged over frequency, i.e. not from averaged Re and Im parts!)
//    computed for the whole energy range except current bin
      if( exclude_energy && photar_sw == -7 ){
        for(ie=0;ie<ne;ie++){
          y1 = ampl_etot2[ie+ne*(ifr0-1)];
          y2 = ampl_etot2[ie+ne*ifr0];
          y3 = ampl_etot2[ie+ne*(ifrn-1)];
          y4 = ampl_etot2[ie+ne*ifrn];
          if( ifr0 < ifrn ){
            ampl2_fband_etot2[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                                   ampl_etot2[ie+ne*ifr0] ) / 2. *
                                 ( freq[ifr0] - frequency1 );
            ampl2_fband_etot2[ie] += ( ampl_etot2[ie+ne*(ifrn-1)] + 
                                    (utmp1 * y3 + utmp * y4) ) / 2. *
                                  ( frequency2 - freq[ifrn-1] );
          }else if( ifr0 == ifrn ){
            ampl2_fband_etot2[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                                   (utmp1 * y3 + utmp * y4) ) / 2. *
                                 ( frequency2 - frequency1 );
          }
//        ...all the whole bins
          for(ifr=ifr0+1;ifr<ifrn;ifr++) 
            ampl2_fband_etot2[ie] += ( ampl_etot2[ie+ne*ifr] + 
                                  ampl_etot2[ie+ne*(ifr-1)]) / 2. * freq[1];
          ampl2_fband_etot2[ie] /= (frequency2 - frequency1);
        }
      }
//    delay (averaged over frequency, i.e. not from averaged Re and Im parts!)
//    computed for the whole energy range except current bin
      if( exclude_energy && photar_sw == -8 ){
        for(ie=0;ie<ne;ie++){
          y1 = phase_etot2[ie+ne*(ifr0-1)] / freq[ifr0-1];
          y2 = phase_etot2[ie+ne*ifr0] / freq[ifr0];
          y3 = phase_etot2[ie+ne*(ifrn-1)] / freq[ifrn-1];
          y4 = phase_etot2[ie+ne*ifrn] / freq[ifrn];
          if( ifr0 < ifrn ){
            delay2_fband_etot2[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                                 phase_etot2[ie+ne*ifr0] / freq[ifr0] ) / 2. *
                               ( freq[ifr0] - frequency1 );
            delay2_fband_etot2[ie] += ( phase_etot2[ie+ne*(ifrn-1)] / freq[ifrn-1] + 
                                  (utmp1 * y3 + utmp * y4) ) / 2. *
                                ( frequency2 - freq[ifrn-1] );
          }else if( ifr0 == ifrn ){
            delay2_fband_etot2[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                                 (utmp1 * y3 + utmp * y4) ) / 2. *
                               ( frequency2 - frequency1 );
          }
//        ...all the whole bins
          for(ifr=ifr0+1;ifr<ifrn;ifr++)
            delay2_fband_etot2[ie] += ( phase_etot2[ie+ne*ifr] / freq[ifr] + 
                                  phase_etot2[ie+ne*(ifr-1)] / freq[ifr-1] ) / 2. 
                                  * freq[1];
          delay2_fband_etot2[ie] /= ( (frequency2 - frequency1) * PI2 );
        }
      }
    }else{
      ampl2_fband_etot = 1.;
      delay2_fband_etot = 0.;
      ampl_fband_etot = 1.;
      phase_fband_etot = 0.;
      if( photar_sw == -7 || photar_sw == -8 )
        for(ie=0;ie<ne;ie++){
          ampl2_fband_etot2[ie] = 1.;
          delay2_fband_etot2[ie] = 0.;
        }
    }
  }
//including primary
#ifndef OUTSIDE_XSPEC
  if( photar_sw > 0){
#endif
    for(ie=0;ie<ne;ie++){
//  real part
      y1 = r_part_tot[ie+ne*(ifr0-1)];
      y2 = r_part_tot[ie+ne*ifr0];
      y3 = r_part_tot[ie+ne*(ifrn-1)];
      y4 = r_part_tot[ie+ne*ifrn];
      if( ifr0 < ifrn ){
        r_part_tot_fband[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                                 r_part_tot[ie+ne*ifr0] ) / 2. *
                               ( freq[ifr0] - frequency1 );
        r_part_tot_fband[ie] += ( r_part_tot[ie+ne*(ifrn-1)] + 
                                  (utmp1 * y3 + utmp * y4) ) / 2. *
                                ( frequency2 - freq[ifrn-1] );
      }else if( ifr0 == ifrn ){
        r_part_tot_fband[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                                 (utmp1 * y3 + utmp * y4) ) / 2. *
                               ( frequency2 - frequency1 );
      }
//    ...all the whole bins
      for(ifr=ifr0+1;ifr<ifrn;ifr++) 
        r_part_tot_fband[ie] += ( r_part_tot[ie+ne*ifr] + 
                                  r_part_tot[ie+ne*(ifr-1)] ) / 2. * freq[1];
//    imaginary part
      y1 = im_part_tot[ie+ne*(ifr0-1)];
      y2 = im_part_tot[ie+ne*ifr0];
      y3 = im_part_tot[ie+ne*(ifrn-1)];
      y4 = im_part_tot[ie+ne*ifrn];
      if( ifr0 < ifrn ){
        im_part_tot_fband[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                                  im_part_tot[ie+ne*ifr0] ) / 2. *
                                ( freq[ifr0] - frequency1 );
        im_part_tot_fband[ie] += ( im_part_tot[ie+ne*(ifrn-1)] + 
                                   (utmp1 * y3 + utmp * y4) ) / 2. *
                                 ( frequency2 - freq[ifrn-1] );
      }else if( ifr0 == ifrn ){
        im_part_tot_fband[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                                  (utmp1 * y3 + utmp * y4) ) / 2. *
                                ( frequency2 - frequency1 );
      }
//    ...all the whole bins
      for(ifr=ifr0+1;ifr<ifrn;ifr++) 
        im_part_tot_fband[ie] += ( im_part_tot[ie+ne*ifr] + 
                                   im_part_tot[ie+ne*(ifr-1)]) / 2. * freq[1];
//    amplitude (averaged over frequency, i.e. not from averaged Re and Im parts!)
      y1 = ampl_tot[ie+ne*(ifr0-1)];
      y2 = ampl_tot[ie+ne*ifr0];
      y3 = ampl_tot[ie+ne*(ifrn-1)];
      y4 = ampl_tot[ie+ne*ifrn];
      if( ifr0 < ifrn ){
        ampl2_tot_fband[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                                  ampl_tot[ie+ne*ifr0] ) / 2. *
                                ( freq[ifr0] - frequency1 );
        ampl2_tot_fband[ie] += ( ampl_tot[ie+ne*(ifrn-1)] + 
                                   (utmp1 * y3 + utmp * y4) ) / 2. *
                                 ( frequency2 - freq[ifrn-1] );
      }else if( ifr0 == ifrn ){
        ampl2_tot_fband[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                                  (utmp1 * y3 + utmp * y4) ) / 2. *
                                ( frequency2 - frequency1 );
      }
//    ...all the whole bins
      for(ifr=ifr0+1;ifr<ifrn;ifr++) 
        ampl2_tot_fband[ie] += ( ampl_tot[ie+ne*ifr] + 
                                 ampl_tot[ie+ne*(ifr-1)]) / 2. * freq[1];
//    delay (averaged over frequency, i.e. not from averaged Re and Im parts!)
      y1 = phase_tot[ie+ne*(ifr0-1)] / freq[ifr0-1];
      y2 = phase_tot[ie+ne*ifr0] / freq[ifr0];
      y3 = phase_tot[ie+ne*(ifrn-1)] / freq[ifrn-1];
      y4 = phase_tot[ie+ne*ifrn] / freq[ifrn];
      if( ifr0 < ifrn ){
        delay2_tot_fband[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                                phase_tot[ie+ne*ifr0] / freq[ifr0] ) / 2. *
                              ( freq[ifr0] - frequency1 );
        delay2_tot_fband[ie] += ( phase_tot[ie+ne*(ifrn-1)] / freq[ifrn-1] + 
                                 (utmp1 * y3 + utmp * y4) ) / 2. *
                               ( frequency2 - freq[ifrn-1] );
      }else if( ifr0 == ifrn ){
        delay2_tot_fband[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                                (utmp1 * y3 + utmp * y4) ) / 2. *
                               ( frequency2 - frequency1 );
      }
//    ...all the whole bins
      for(ifr=ifr0+1;ifr<ifrn;ifr++) 
        delay2_tot_fband[ie] += ( phase_tot[ie+ne*ifr] / freq[ifr] + 
                                  phase_tot[ie+ne*(ifr-1)] / freq[ifr-1] ) / 2. 
                                  * freq[1];
    }
    for(ie=0;ie<ne;ie++){
      r_part_tot_fband[ie] /= (frequency2 - frequency1);
      im_part_tot_fband[ie] /= (frequency2 - frequency1);
      ampl2_tot_fband[ie] /= (frequency2 - frequency1);
      delay2_tot_fband[ie] /= ( (frequency2 - frequency1) * PI2 );
    }
    for(ie=0;ie<ne;ie++){
      k2=0;
      ampl_tot_fband[ie] = sqrt( r_part_tot_fband[ie]  * r_part_tot_fband[ie] +
                           im_part_tot_fband[ie] * im_part_tot_fband[ie] );
      phase_tot_fband[ie] = atan2(im_part_tot_fband[ie], r_part_tot_fband[ie]);
      phase_tot_fband_u2[ie] = atan2(im_part_tot_fband[ie], r_part_tot_fband[ie])+k2*PI2;
      if(ie>0 && (phase_tot_fband_u2[ie]-phase_tot_fband_u2[ie-1]) < -PI) k2++;
      if(ie>0 && (phase_tot_fband_u2[ie]-phase_tot_fband_u2[ie-1]) > PI) k2--;
      phase_tot_fband_u2[ie] = atan2(im_part_tot_fband[ie], r_part_tot_fband[ie])+k2*PI2;
    }
//reference energy band Fourier transform integrated in the given frequency band
    if(nbands>0){
//    real part
      y1 = r_part_bands_tot[nbands-1+nbands*(ifr0-1)];
      y2 = r_part_bands_tot[nbands-1+nbands*ifr0];
      y3 = r_part_bands_tot[nbands-1+nbands*(ifrn-1)];
      y4 = r_part_bands_tot[nbands-1+nbands*ifrn];
      if( ifr0 < ifrn ){
        r_part_tot_fband_etot = ( (ttmp1 * y1 + ttmp * y2) + 
                                  r_part_bands_tot[nbands-1+nbands*ifr0] ) / 2. *
                                ( freq[ifr0] - frequency1 );
        r_part_tot_fband_etot += ( r_part_bands_tot[nbands-1+nbands*(ifrn-1)] + 
                                   (utmp1 * y3 + utmp * y4) ) / 2. *
                                 ( frequency2 - freq[ifrn-1] );
      }else if( ifr0 == ifrn ){
        r_part_tot_fband_etot = ( (ttmp1 * y1 + ttmp * y2) + 
                                  (utmp1 * y3 + utmp * y4) ) / 2. *
                                ( frequency2 - frequency1 );
      }
//    ...all the whole bins
      for(ifr=ifr0+1;ifr<ifrn;ifr++) 
        r_part_tot_fband_etot += ( r_part_bands_tot[nbands-1+nbands*ifr] + 
                                   r_part_bands_tot[nbands-1+nbands*(ifr-1)]) / 2. * freq[1];
//    imaginary part
      y1 = im_part_bands_tot[nbands-1+nbands*(ifr0-1)];
      y2 = im_part_bands_tot[nbands-1+nbands*ifr0];
      y3 = im_part_bands_tot[nbands-1+nbands*(ifrn-1)];
      y4 = im_part_bands_tot[nbands-1+nbands*ifrn];
      if( ifr0 < ifrn ){
        im_part_tot_fband_etot = ( (ttmp1 * y1 + ttmp * y2) + 
                                   im_part_bands_tot[nbands-1+nbands*ifr0] ) / 2. *
                                 ( freq[ifr0] - frequency1 );
        im_part_tot_fband_etot += ( im_part_bands_tot[nbands-1+nbands*(ifrn-1)] + 
                                    (utmp1 * y3 + utmp * y4) ) / 2. *
                                  ( frequency2 - freq[ifrn-1] );
      }else if( ifr0 == ifrn ){
        im_part_tot_fband_etot = ( (ttmp1 * y1 + ttmp * y2) + 
                                   (utmp1 * y3 + utmp * y4) ) / 2. *
                                 ( frequency2 - frequency1 );
      }
//    ...all the whole bins
      for(ifr=ifr0+1;ifr<ifrn;ifr++) 
        im_part_tot_fband_etot += ( im_part_bands_tot[nbands-1+nbands*ifr] + 
                                    im_part_bands_tot[nbands-1+nbands*(ifr-1)]) / 2. * freq[1];
//    amplitude (averaged over frequency, i.e. not from averaged Re and Im parts!)
      y1 = ampl_bands_tot[nbands-1+nbands*(ifr0-1)];
      y2 = ampl_bands_tot[nbands-1+nbands*ifr0];
      y3 = ampl_bands_tot[nbands-1+nbands*(ifrn-1)];
      y4 = ampl_bands_tot[nbands-1+nbands*ifrn];
      if( ifr0 < ifrn ){
        ampl2_tot_fband_etot = ( (ttmp1 * y1 + ttmp * y2) + 
                                   ampl_bands_tot[nbands-1+nbands*ifr0] ) / 2. *
                                 ( freq[ifr0] - frequency1 );
        ampl2_tot_fband_etot += ( ampl_bands_tot[nbands-1+nbands*(ifrn-1)] + 
                                    (utmp1 * y3 + utmp * y4) ) / 2. *
                                  ( frequency2 - freq[ifrn-1] );
      }else if( ifr0 == ifrn ){
        ampl2_tot_fband_etot = ( (ttmp1 * y1 + ttmp * y2) + 
                                   (utmp1 * y3 + utmp * y4) ) / 2. *
                                 ( frequency2 - frequency1 );
      }
//    ...all the whole bins
      for(ifr=ifr0+1;ifr<ifrn;ifr++) 
        ampl2_tot_fband_etot += ( ampl_bands_tot[nbands-1+nbands*ifr] + 
                                  ampl_bands_tot[nbands-1+nbands*(ifr-1)]) / 2. * freq[1];
//    delay (averaged over frequency, i.e. not from averaged Re and Im parts!)
      y1 = phase_bands_tot[nbands-1+nbands*(ifr0-1)] / freq[ifr0-1];
      y2 = phase_bands_tot[nbands-1+nbands*ifr0] / freq[ifr0];
      y3 = phase_bands_tot[nbands-1+nbands*(ifrn-1)] / freq[ifrn-1];
      y4 = phase_bands_tot[nbands-1+nbands*ifrn] / freq[ifrn];
      if( ifr0 < ifrn ){
        delay2_tot_fband_etot = ( (ttmp1 * y1 + ttmp * y2) + 
                                 phase_bands_tot[nbands-1+nbands*ifr0] / freq[ifr0] ) / 2. *
                               ( freq[ifr0] - frequency1 );
        delay2_tot_fband_etot += ( phase_bands_tot[nbands-1+nbands*(ifrn-1)] / freq[ifrn-1] + 
                                  (utmp1 * y3 + utmp * y4) ) / 2. *
                                ( frequency2 - freq[ifrn-1] );
      }else if( ifr0 == ifrn ){
        delay2_tot_fband_etot = ( (ttmp1 * y1 + ttmp * y2) + 
                                 (utmp1 * y3 + utmp * y4) ) / 2. *
                               ( frequency2 - frequency1 );
      }
//    ...all the whole bins
      for(ifr=ifr0+1;ifr<ifrn;ifr++)
        delay2_tot_fband_etot += ( phase_bands_tot[nbands-1+nbands*ifr] / freq[ifr] + 
                                  phase_bands_tot[nbands-1+nbands*(ifr-1)] / freq[ifr-1] ) / 2. 
                                  * freq[1];
      r_part_tot_fband_etot /= (frequency2 - frequency1);
      im_part_tot_fband_etot /= (frequency2 - frequency1);
      ampl2_tot_fband_etot /= (frequency2 - frequency1);
      delay2_tot_fband_etot /= ( (frequency2 - frequency1) * PI2 );
      ampl_tot_fband_etot  = sqrt( r_part_tot_fband_etot  * r_part_tot_fband_etot +
                                   im_part_tot_fband_etot * im_part_tot_fband_etot );
      phase_tot_fband_etot = atan2(im_part_tot_fband_etot, r_part_tot_fband_etot);
//    amplitude (averaged over frequency, i.e. not from averaged Re and Im parts!)
//    computed for the whole energy range except current bin
      if( exclude_energy && photar_sw == 7){
        for(ie=0;ie<ne;ie++){
          y1 = ampl_tot_etot2[ie+ne*(ifr0-1)];
          y2 = ampl_tot_etot2[ie+ne*ifr0];
          y3 = ampl_tot_etot2[ie+ne*(ifrn-1)];
          y4 = ampl_tot_etot2[ie+ne*ifrn];
          if( ifr0 < ifrn ){
            ampl2_tot_fband_etot2[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                                       ampl_tot_etot2[ie+ne*ifr0] ) / 2. *
                                     ( freq[ifr0] - frequency1 );
            ampl2_tot_fband_etot2[ie] += ( ampl_tot_etot2[ie+ne*(ifrn-1)] + 
                                        (utmp1 * y3 + utmp * y4) ) / 2. *
                                      ( frequency2 - freq[ifrn-1] );
          }else if( ifr0 == ifrn ){
            ampl2_tot_fband_etot2[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                                       (utmp1 * y3 + utmp * y4) ) / 2. *
                                     ( frequency2 - frequency1 );
          }
//        ...all the whole bins
          for(ifr=ifr0+1;ifr<ifrn;ifr++) 
            ampl2_tot_fband_etot2[ie] += ( ampl_tot_etot2[ie+ne*ifr] + 
                                           ampl_tot_etot2[ie+ne*(ifr-1)]) / 2. * freq[1];
          ampl2_tot_fband_etot2[ie] /= (frequency2 - frequency1);
        }
      }
//    delay (averaged over frequency, i.e. not from averaged Re and Im parts!)
//    computed for the whole energy range except current bin
      if( exclude_energy && photar_sw == 8 ){
        for(ie=0;ie<ne;ie++){
          y1 = phase_tot_etot2[ie+ne*(ifr0-1)] / freq[ifr0-1];
          y2 = phase_tot_etot2[ie+ne*ifr0] / freq[ifr0];
          y3 = phase_tot_etot2[ie+ne*(ifrn-1)] / freq[ifrn-1];
          y4 = phase_tot_etot2[ie+ne*ifrn] / freq[ifrn];
          if( ifr0 < ifrn ){
            delay2_tot_fband_etot2[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                                     phase_tot_etot2[ie+ne*ifr0] / freq[ifr0] ) / 2. *
                                   ( freq[ifr0] - frequency1 );
            delay2_tot_fband_etot2[ie] += ( phase_tot_etot2[ie+ne*(ifrn-1)] / freq[ifrn-1] + 
                                      (utmp1 * y3 + utmp * y4) ) / 2. *
                                    ( frequency2 - freq[ifrn-1] );
          }else if( ifr0 == ifrn ){
            delay2_tot_fband_etot2[ie] = ( (ttmp1 * y1 + ttmp * y2) + 
                                     (utmp1 * y3 + utmp * y4) ) / 2. *
                                   ( frequency2 - frequency1 );
          }
//        ...all the whole bins
          for(ifr=ifr0+1;ifr<ifrn;ifr++)
            delay2_tot_fband_etot2[ie] += ( phase_tot_etot2[ie+ne*ifr] / freq[ifr] + 
                                      phase_tot_etot2[ie+ne*(ifr-1)] / freq[ifr-1] ) / 2. 
                                      * freq[1];
          delay2_tot_fband_etot2[ie] /= ( (frequency2 - frequency1) * PI2 );
        }
      }
    }else{
      ampl2_tot_fband_etot = 1.;
      delay2_tot_fband_etot = 0.;
      ampl_tot_fband_etot = 1.;
      phase_tot_fband_etot = 0.;
      if( photar_sw == 7 || photar_sw == 8 )
        for(ie=0;ie<ne;ie++){
          ampl2_tot_fband_etot2[ie] = 1.;
          delay2_tot_fband_etot2[ie] = 0.;
        }
    }
#ifndef OUTSIDE_XSPEC
  }
#endif
}
//Let's rebin the frequency dependent FT to the xspec frequencies
//Note that ear[] contains frequencies and ne denotes number of frequency bins!
if( abs(photar_sw) > 10 && nbands > 0){
  iiband=0;
  for( iband = 0; iband < nbands; iband += ( nbands > 1 ? (nbands-1) : 1 ) ){
    for(ifreq=0;ifreq<ne;ifreq++){
      if(photar_sw < 0){
        if(photar_sw == -17)ampl2_bands_rebin[iiband+2*ifreq]=0.;
        else if(photar_sw == -18)delay2_bands_rebin[iiband+2*ifreq]=0.;
        else{
          r_part_bands_rebin[iiband+2*ifreq]=0.;
          im_part_bands_rebin[iiband+2*ifreq]=0.;
        }
      }else if(photar_sw > 0){
        if(photar_sw == 17)ampl2_bands_tot_rebin[iiband+2*ifreq]=0.;
        else if(photar_sw == 18)delay2_bands_tot_rebin[iiband+2*ifreq]=0.;
        else{
          r_part_bands_tot_rebin[iiband+2*ifreq]=0.;
          im_part_bands_tot_rebin[iiband+2*ifreq]=0.;
        }
      }
//    given frequency1 given by ear[ifreq] and frequency2 at ear[ifreq+1], 
//    find the corresponding index in freq[]:
      frequency1=ear[ifreq];
      frequency2=ear[ifreq+1];
      if(frequency2 <= freq[1] || frequency1 >= freq[nn/2])continue;
      if(frequency1 < freq[1])frequency1=freq[1];
      if(frequency2 > freq[nn/2])frequency2=freq[nn/2];
      ifr0 = (int) ceil( frequency1 / freq[1] );
      if(ifr0 < 2) ifr0 = 2;
//      else if(ifr0 > nn/2) ifr0 = nn/2;
      ifrn = (int) ceil(frequency2 / freq[1]);
//      if(ifrn < 2) ifrn = 2;
//      else 
      if(ifrn > nn/2) ifrn = nn/2;
//    the first and last partial bin
      ttmp = (frequency1 - freq[ifr0 - 1]) / (freq[ifr0] - freq[ifr0 - 1]);
      ttmp1 = 1. - ttmp;
      utmp = (frequency2 - freq[ifrn - 1]) / (freq[ifrn] - freq[ifrn - 1]);
      utmp1 = 1. - utmp;
//    just reflection
      if(photar_sw < 0){
        if(photar_sw == -17){
//      direct amplitude rebinnig
          y1 = ampl_bands[iband+nbands*(ifr0-1)];
          y2 = ampl_bands[iband+nbands*ifr0];
          y3 = ampl_bands[iband+nbands*(ifrn-1)];
          y4 = ampl_bands[iband+nbands*ifrn];
          if( ifr0 < ifrn ){
            ampl2_bands_rebin[iiband+2*ifreq] += 
              ( (ttmp1 * y1 + ttmp * y2) + ampl_bands[iband+nbands*ifr0] )
              / 2. * ( freq[ifr0] - frequency1 );
            ampl2_bands_rebin[iiband+2*ifreq] += 
              ( ampl_bands[iband+nbands*(ifrn-1)] + (utmp1 * y3 + utmp * y4) )
              / 2. * ( frequency2 - freq[ifrn-1] );
          }else if( ifr0 == ifrn ){
            ampl2_bands_rebin[iiband+2*ifreq] += 
              ( (ttmp1 * y1 + ttmp * y2) + (utmp1 * y3 + utmp * y4) )
              / 2. * ( frequency2 - frequency1 );
          }
//        ...all the whole bins
          for(ifr=ifr0+1;ifr<ifrn;ifr++) 
            ampl2_bands_rebin[iiband+2*ifreq] += 
              ( ampl_bands[iband+nbands*ifr] + ampl_bands[iband+nbands*(ifr-1)]) 
              / 2. * freq[1];
        }else if(photar_sw == -18){
//      direct delay rebinnig
          y1 = phase_bands[iband+nbands*(ifr0-1)] / PI2 / freq[ifr0-1];
          y2 = phase_bands[iband+nbands*ifr0] / PI2 / freq[ifr0];
          y3 = phase_bands[iband+nbands*(ifrn-1)] / PI2 / freq[ifrn-1];
          y4 = phase_bands[iband+nbands*ifrn] / PI2 / freq[ifrn];
          if( ifr0 < ifrn ){
            delay2_bands_rebin[iiband+2*ifreq] += 
              ( (ttmp1 * y1 + ttmp * y2) + phase_bands[iband+nbands*ifr0] / PI2 / freq[ifr0])
              / 2. * ( freq[ifr0] - frequency1 );
            delay2_bands_rebin[iiband+2*ifreq] += 
              ( phase_bands[iband+nbands*(ifrn-1)] / PI2 / freq[ifrn-1] + 
              (utmp1 * y3 + utmp * y4) ) / 2. * ( frequency2 - freq[ifrn-1] );
          }else if( ifr0 == ifrn ){
            delay2_bands_rebin[iiband+2*ifreq] += 
              ( (ttmp1 * y1 + ttmp * y2) + (utmp1 * y3 + utmp * y4) )
              / 2. * ( frequency2 - frequency1 );
          }
//        ...all the whole bins
          for(ifr=ifr0+1;ifr<ifrn;ifr++) 
            delay2_bands_rebin[iiband+2*ifreq] += 
              ( phase_bands[iband+nbands*ifr] / freq[ifr] + phase_bands[iband+nbands*(ifr-1)] / freq[ifr-1]) / PI2 
              / 2. * freq[1];
        }else{
//        real part
          y1 = r_part_bands[iband+nbands*(ifr0-1)];
          y2 = r_part_bands[iband+nbands*ifr0];
          y3 = r_part_bands[iband+nbands*(ifrn-1)];
          y4 = r_part_bands[iband+nbands*ifrn];
          if( ifr0 < ifrn ){
            r_part_bands_rebin[iiband+2*ifreq] += 
              ( (ttmp1 * y1 + ttmp * y2) + r_part_bands[iband+nbands*ifr0] )
              / 2. * ( freq[ifr0] - frequency1 );
            r_part_bands_rebin[iiband+2*ifreq] += 
              ( r_part_bands[iband+nbands*(ifrn-1)] + (utmp1 * y3 + utmp * y4) )
              / 2. * ( frequency2 - freq[ifrn-1] );
          }else if( ifr0 == ifrn ){
            r_part_bands_rebin[iiband+2*ifreq] += 
              ( (ttmp1 * y1 + ttmp * y2) + (utmp1 * y3 + utmp * y4) )
              / 2. * ( frequency2 - frequency1 );
          }
//        ...all the whole bins
          for(ifr=ifr0+1;ifr<ifrn;ifr++) 
            r_part_bands_rebin[iiband+2*ifreq] += 
              ( r_part_bands[iband+nbands*ifr] + r_part_bands[iband+nbands*(ifr-1)]) 
              / 2. * freq[1];
//        imaginary part
          y1 = im_part_bands[iband+nbands*(ifr0-1)];
          y2 = im_part_bands[iband+nbands*ifr0];
          y3 = im_part_bands[iband+nbands*(ifrn-1)];
          y4 = im_part_bands[iband+nbands*ifrn];
          if( ifr0 < ifrn ){
            im_part_bands_rebin[iiband+2*ifreq] += 
              ( (ttmp1 * y1 + ttmp * y2) + im_part_bands[iband+nbands*ifr0] )
              / 2. * ( freq[ifr0] - frequency1 );
            im_part_bands_rebin[iiband+2*ifreq] += 
              ( im_part_bands[iband+nbands*(ifrn-1)] + (utmp1 * y3 + utmp * y4) )
              / 2. * ( frequency2 - freq[ifrn-1] );
          }else if( ifr0 == ifrn ){
            im_part_bands_rebin[iiband+2*ifreq] += 
            ( (ttmp1 * y1 + ttmp * y2) + (utmp1 * y3 + utmp * y4) )
            / 2. * ( frequency2 - frequency1 );
          }
//        ...all the whole bins
          for(ifr=ifr0+1;ifr<ifrn;ifr++) 
            im_part_bands_rebin[iiband+2*ifreq] += 
              ( im_part_bands[iband+nbands*ifr] + im_part_bands[iband+nbands*(ifr-1)])
              / 2. * freq[1];
          ampl_bands_rebin[iiband+2*ifreq]  = 
            sqrt( r_part_bands_rebin[iiband+2*ifreq]  * r_part_bands_rebin[iiband+2*ifreq] +
                  im_part_bands_rebin[iiband+2*ifreq] * im_part_bands_rebin[iiband+2*ifreq] );
          phase_bands_rebin[iiband+2*ifreq] = atan2(im_part_bands_rebin[iiband+2*ifreq], r_part_bands_rebin[iiband+2*ifreq]);
        }
//    including primary
      }else if(photar_sw > 0){
        if(photar_sw == 17){
//      direct amplitude rebinnig
          y1 = ampl_bands_tot[iband+nbands*(ifr0-1)];
          y2 = ampl_bands_tot[iband+nbands*ifr0];
          y3 = ampl_bands_tot[iband+nbands*(ifrn-1)];
          y4 = ampl_bands_tot[iband+nbands*ifrn];
          if( ifr0 < ifrn ){
            ampl2_bands_tot_rebin[iiband+2*ifreq] += 
              ( (ttmp1 * y1 + ttmp * y2) + ampl_bands_tot[iband+nbands*ifr0] )
              / 2. * ( freq[ifr0] - frequency1 );
            ampl2_bands_tot_rebin[iiband+2*ifreq] += 
              ( ampl_bands_tot[iband+nbands*(ifrn-1)] + (utmp1 * y3 + utmp * y4) )
              / 2. * ( frequency2 - freq[ifrn-1] );
          }else if( ifr0 == ifrn ){
            ampl2_bands_tot_rebin[iiband+2*ifreq] += 
              ( (ttmp1 * y1 + ttmp * y2) + (utmp1 * y3 + utmp * y4) )
              / 2. * ( frequency2 - frequency1 );
          }
//        ...all the whole bins
          for(ifr=ifr0+1;ifr<ifrn;ifr++) 
            ampl2_bands_tot_rebin[iiband+2*ifreq] += 
              ( ampl_bands_tot[iband+nbands*ifr] + ampl_bands_tot[iband+nbands*(ifr-1)]) 
              / 2. * freq[1];
        }else if(photar_sw == 18){
//      direct delay rebinnig
          y1 = phase_bands_tot[iband+nbands*(ifr0-1)] / PI2 / freq[ifr0-1];
          y2 = phase_bands_tot[iband+nbands*ifr0] / PI2 / freq[ifr0];
          y3 = phase_bands_tot[iband+nbands*(ifrn-1)] / PI2 / freq[ifrn-1];
          y4 = phase_bands_tot[iband+nbands*ifrn] / PI2 / freq[ifrn];
          if( ifr0 < ifrn ){
            delay2_bands_tot_rebin[iiband+2*ifreq] += 
              ( (ttmp1 * y1 + ttmp * y2) + phase_bands_tot[iband+nbands*ifr0] / PI2 / freq[ifr0] )
              / 2. * ( freq[ifr0] - frequency1 );
            delay2_bands_tot_rebin[iiband+2*ifreq] += 
              ( phase_bands_tot[iband+nbands*(ifrn-1)] / PI2 / freq[ifrn-1] + (utmp1 * y3 + utmp * y4) )
              / 2. * ( frequency2 - freq[ifrn-1] );
          }else if( ifr0 == ifrn ){
            delay2_bands_tot_rebin[iiband+2*ifreq] += 
              ( (ttmp1 * y1 + ttmp * y2) + (utmp1 * y3 + utmp * y4) )
              / 2. * ( frequency2 - frequency1 );
          }
//        ...all the whole bins
          for(ifr=ifr0+1;ifr<ifrn;ifr++) 
            delay2_bands_tot_rebin[iiband+2*ifreq] += 
              ( phase_bands_tot[iband+nbands*ifr] / freq[ifr] + phase_bands_tot[iband+nbands*(ifr-1)] / freq[ifr-1] ) / PI2 
              / 2. * freq[1];
        }else{
//        real part
          y1 = r_part_bands_tot[iband+nbands*(ifr0-1)];
          y2 = r_part_bands_tot[iband+nbands*ifr0];
          y3 = r_part_bands_tot[iband+nbands*(ifrn-1)];
          y4 = r_part_bands_tot[iband+nbands*ifrn];
          if( ifr0 < ifrn ){
            r_part_bands_tot_rebin[iiband+2*ifreq] += 
              ( (ttmp1 * y1 + ttmp * y2) + r_part_bands_tot[iband+nbands*ifr0] )
              / 2. * ( freq[ifr0] - frequency1 );
            r_part_bands_tot_rebin[iiband+2*ifreq] += 
              ( r_part_bands_tot[iband+nbands*(ifrn-1)] + (utmp1 * y3 + utmp * y4) )
              / 2. * ( frequency2 - freq[ifrn-1] );
          }else if( ifr0 == ifrn ){
            r_part_bands_tot_rebin[iiband+2*ifreq] += 
              ( (ttmp1 * y1 + ttmp * y2) + (utmp1 * y3 + utmp * y4) )
              / 2. * ( frequency2 - frequency1 );
          }
//        ...all the whole bins
          for(ifr=ifr0+1;ifr<ifrn;ifr++) 
            r_part_bands_tot_rebin[iiband+2*ifreq] += 
              ( r_part_bands_tot[iband+nbands*ifr] + r_part_bands_tot[iband+nbands*(ifr-1)]) 
              / 2. * freq[1];
//        imaginary part
          y1 = im_part_bands_tot[iband+nbands*(ifr0-1)];
          y2 = im_part_bands_tot[iband+nbands*ifr0];
          y3 = im_part_bands_tot[iband+nbands*(ifrn-1)];
          y4 = im_part_bands_tot[iband+nbands*ifrn];
          if( ifr0 < ifrn ){
            im_part_bands_tot_rebin[iiband+2*ifreq] += 
              ( (ttmp1 * y1 + ttmp * y2) + im_part_bands_tot[iband+nbands*ifr0] )
              / 2. * ( freq[ifr0] - frequency1 );
            im_part_bands_tot_rebin[iiband+2*ifreq] += 
              ( im_part_bands_tot[iband+nbands*(ifrn-1)] + (utmp1 * y3 + utmp * y4) )
              / 2. * ( frequency2 - freq[ifrn-1] );
          }else if( ifr0 == ifrn ){
            im_part_bands_tot_rebin[iiband+2*ifreq] += 
            ( (ttmp1 * y1 + ttmp * y2) + (utmp1 * y3 + utmp * y4) )
            / 2. * ( frequency2 - frequency1 );
          }
//        ...all the whole bins
          for(ifr=ifr0+1;ifr<ifrn;ifr++) 
            im_part_bands_tot_rebin[iiband+2*ifreq] += 
              ( im_part_bands_tot[iband+nbands*ifr] + im_part_bands_tot[iband+nbands*(ifr-1)])
              / 2. * freq[1];
          ampl_bands_tot_rebin[iiband+2*ifreq]  = 
            sqrt( r_part_bands_tot_rebin[iiband+2*ifreq]  * r_part_bands_tot_rebin[iiband+2*ifreq] +
                  im_part_bands_tot_rebin[iiband+2*ifreq] * im_part_bands_tot_rebin[iiband+2*ifreq] );
          phase_bands_tot_rebin[iiband+2*ifreq] = atan2(im_part_bands_tot_rebin[iiband+2*ifreq], r_part_bands_tot_rebin[iiband+2*ifreq]);
        }
      }
    }
    iiband++;
  }
}
/******************************************************************************/
#ifdef OUTSIDE_XSPEC
// Write the tables
// energy dependent Fourier transform
if(abs(photar_sw)<=10){
  sprintf(filename, "kynrefrev_%s_real.dat",cparam);
  fw1 = fopen(filename, "w");
  sprintf(filename, "kynrefrev_%s_imag.dat",cparam);
  fw2 = fopen(filename, "w");
  sprintf(filename, "kynrefrev_%s_ampl.dat",cparam);
  fw3 = fopen(filename, "w");
  sprintf(filename, "kynrefrev_%s_phase.dat",cparam);
  fw4 = fopen(filename, "w");
  sprintf(filename, "kynrefrev_%s_phase_u1.dat",cparam);
  fw5 = fopen(filename, "w");
  sprintf(filename, "kynrefrev_%s_phase_u2.dat",cparam);
  fw6 = fopen(filename, "w");
  for(ifr=0;ifr<=nn/2;ifr++){
    for(ie=0;ie<ne;ie++){
      fprintf(fw1, "%E\t", r_part[ie+ne*ifr]);
      fprintf(fw2, "%E\t", im_part[ie+ne*ifr]);
      fprintf(fw3, "%E\t", ampl[ie+ne*ifr]);
      fprintf(fw4, "%E\t", phase[ie+ne*ifr]);
      fprintf(fw5, "%E\t", phase_u1[ie+ne*ifr]);
      fprintf(fw6, "%E\t", phase_u2[ie+ne*ifr]);
    }
    fprintf(fw1, "\n");
    fprintf(fw2, "\n");
    fprintf(fw3, "\n");
    fprintf(fw4, "\n");
    fprintf(fw5, "\n");
    fprintf(fw6, "\n");
  }
  fclose(fw1);
  fclose(fw2);
  fclose(fw3);
  fclose(fw4);
  fclose(fw5);
  fclose(fw6);
}
// Fourier transform for chosen energy bands
sprintf(filename, "kynrefrev_%s_bands_real.dat",cparam);
fw1 = fopen(filename, "w");
sprintf(filename, "kynrefrev_%s_bands_imag.dat",cparam);
fw2 = fopen(filename, "w");
sprintf(filename, "kynrefrev_%s_bands_ampl.dat",cparam);
fw3 = fopen(filename, "w");
sprintf(filename, "kynrefrev_%s_bands_phase.dat",cparam);
fw4 = fopen(filename, "w");
sprintf(filename, "kynrefrev_%s_bands_phase_u1.dat",cparam);
fw5 = fopen(filename, "w");
//sprintf(filename, "kynrefrev_%s_bands_phase_u2.dat",cparam);
//fw4 = fopen(filename, "w");
for(ifr=0;ifr<=nn/2;ifr++){
  fprintf(fw1, "%E\t", freq[ifr]);
  fprintf(fw2, "%E\t", freq[ifr]);
  fprintf(fw3, "%E\t", freq[ifr]);
  fprintf(fw4, "%E\t", freq[ifr]);
  fprintf(fw5, "%E\t", freq[ifr]);
  for(iband=0;iband<nbands;iband++){
    fprintf(fw1, "%E\t", r_part_bands[iband+nbands*ifr]);
    fprintf(fw2, "%E\t", im_part_bands[iband+nbands*ifr]);
    fprintf(fw3, "%E\t", ampl_bands[iband+nbands*ifr]);
    fprintf(fw4, "%E\t", phase_bands[iband+nbands*ifr]);
    fprintf(fw5, "%E\t", phase_bands_u1[iband+nbands*ifr]);
  }
  fprintf(fw1, "\n");
  fprintf(fw2, "\n");
  fprintf(fw3, "\n");
  fprintf(fw4, "\n");
  fprintf(fw5, "\n");
}
fclose(fw1);
fclose(fw2);
fclose(fw3);
fclose(fw4);
fclose(fw5);
// energy dependent Fourier transform when primary is included
if(abs(photar_sw)<=10){
  sprintf(filename, "kynrefrev_%s_real_tot.dat",cparam);
  fw1 = fopen(filename, "w");
  sprintf(filename, "kynrefrev_%s_imag_tot.dat",cparam);
  fw2 = fopen(filename, "w");
  sprintf(filename, "kynrefrev_%s_ampl_tot.dat",cparam);
  fw3 = fopen(filename, "w");
  sprintf(filename, "kynrefrev_%s_phase_tot.dat",cparam);
  fw4 = fopen(filename, "w");
  sprintf(filename, "kynrefrev_%s_phase_tot_u1.dat",cparam);
  fw5 = fopen(filename, "w");
  sprintf(filename, "kynrefrev_%s_phase_tot_u2.dat",cparam);
  fw6 = fopen(filename, "w");
  for(ifr=0;ifr<=nn/2;ifr++){
    for(ie=0;ie<ne;ie++){
      fprintf(fw1, "%E\t", r_part_tot[ie+ne*ifr]);
      fprintf(fw2, "%E\t", im_part_tot[ie+ne*ifr]);
      fprintf(fw3, "%E\t", ampl_tot[ie+ne*ifr]);
      fprintf(fw4, "%E\t", phase_tot[ie+ne*ifr]);
      fprintf(fw5, "%E\t", phase_tot_u1[ie+ne*ifr]);
      fprintf(fw6, "%E\t", phase_tot_u2[ie+ne*ifr]);
    }
    fprintf(fw1, "\n");
    fprintf(fw2, "\n");
    fprintf(fw3, "\n");
    fprintf(fw4, "\n");
    fprintf(fw5, "\n");
    fprintf(fw6, "\n");
  }
  fclose(fw1);
  fclose(fw2);
  fclose(fw3);
  fclose(fw4);
  fclose(fw5);
  fclose(fw6);
}
// Fourier transform for chosen energy bands when primary is included
sprintf(filename, "kynrefrev_%s_bands_real_tot.dat",cparam);
fw1 = fopen(filename, "w");
sprintf(filename, "kynrefrev_%s_bands_imag_tot.dat",cparam);
fw2 = fopen(filename, "w");
sprintf(filename, "kynrefrev_%s_bands_ampl_tot.dat",cparam);
fw3 = fopen(filename, "w");
sprintf(filename, "kynrefrev_%s_bands_phase_tot.dat",cparam);
fw4 = fopen(filename, "w");
sprintf(filename, "kynrefrev_%s_bands_phase_tot_u1.dat",cparam);
fw5 = fopen(filename, "w");
for(ifr=0;ifr<=nn/2;ifr++){
  fprintf(fw1, "%E\t", freq[ifr]);
  fprintf(fw2, "%E\t", freq[ifr]);
  fprintf(fw3, "%E\t", freq[ifr]);
  fprintf(fw4, "%E\t", freq[ifr]);
  fprintf(fw5, "%E\t", freq[ifr]);
  for(iband=0;iband<nbands;iband++){
    fprintf(fw1, "%E\t", r_part_bands_tot[iband+nbands*ifr]);
    fprintf(fw2, "%E\t", im_part_bands_tot[iband+nbands*ifr]);
    fprintf(fw3, "%E\t", ampl_bands_tot[iband+nbands*ifr]);
    fprintf(fw4, "%E\t", phase_bands_tot[iband+nbands*ifr]);
    fprintf(fw5, "%E\t", phase_bands_tot_u1[iband+nbands*ifr]);
  }
  fprintf(fw1, "\n");
  fprintf(fw2, "\n");
  fprintf(fw3, "\n");
  fprintf(fw4, "\n");
  fprintf(fw5, "\n");
}
fclose(fw1);
fclose(fw2);
fclose(fw3);
fclose(fw4);
fclose(fw5);
// frequency integrated energy dependent functions
if(abs(photar_sw)<=10){
  sprintf(filename, "kynrefrev_%s_fft_tot_int.dat",cparam);
  fw1 = fopen(filename, "w");
  for(ie=0;ie<ne;ie++)
    fprintf(fw1, "%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n", 
            (ear[ie]+ear[ie+1])/2.,
            r_part_tot_int[ie], im_part_tot_int[ie], ampl_tot_int[ie], 
            phase_tot_int[ie], phase_tot_int_u2[ie], delay2_tot_int[ie],
            ampl_tot_int[ie]/ampl_tot_int_etot,
// we will use relative delay with respect to the whole energy band and redefine to be positive
            -(phase_tot_int[ie]-phase_tot_int_etot)/PI/freq_wrap0,
            -(phase_tot_int_u2[ie]-phase_tot_int_etot)/PI/freq_wrap0,
            -(delay2_tot_int[ie]-delay2_tot_int_etot),
// and we will also exclude the current energy bin from total
            ampl_tot_int[ie]/sqrt(
              ( r_part_tot_int_etot- r_part_tot_int[ie])*( r_part_tot_int_etot- r_part_tot_int[ie])+
              (im_part_tot_int_etot-im_part_tot_int[ie])*(im_part_tot_int_etot-im_part_tot_int[ie])),
            -(phase_tot_int[ie]-atan2(
               im_part_tot_int_etot-im_part_tot_int[ie],
                r_part_tot_int_etot- r_part_tot_int[ie]))/PI/freq_wrap0);
  fclose(fw1);
  sprintf(filename, "kynrefrev_%s_freq_wrap.dat",cparam);
  fw1 = fopen(filename, "w");
  fprintf(fw1, "%E", freq_wrap0);
  fclose(fw1);
// frequency band integrated Fourier transform
  if( frequency1 > 0. && frequency1 < frequency2){
  sprintf(filename, "kynrefrev_%s_fft_tot_fband.dat",cparam);
  fw1 = fopen(filename, "w");
  for(ie=0;ie<ne;ie++)
    fprintf(fw1, "%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n", 
            (ear[ie]+ear[ie+1])/2.,
            r_part_tot_fband[ie], im_part_tot_fband[ie], ampl_tot_fband[ie], 
            phase_tot_fband[ie], phase_tot_fband_u2[ie], delay2_tot_fband[ie],
            ampl_tot_fband[ie]/ampl_tot_fband_etot,
// we will use relative delay with respect to the whole energy band and redefine to be positive
            -(phase_tot_fband[ie]-phase_tot_fband_etot)/PI/(frequency1+frequency2),
            -(phase_tot_fband_u2[ie]-phase_tot_fband_etot)/PI/(frequency1+frequency2), 
            -(delay2_tot_fband[ie]-delay2_tot_fband_etot),
// and we will also exclude the current energy bin from total
            ampl_tot_fband[ie]/sqrt(
              ( r_part_tot_fband_etot- r_part_tot_fband[ie])*( r_part_tot_fband_etot- r_part_tot_fband[ie])+
              (im_part_tot_fband_etot-im_part_tot_fband[ie])*(im_part_tot_fband_etot-im_part_tot_fband[ie])),
            -(phase_tot_fband[ie]-atan2(
               im_part_tot_fband_etot-im_part_tot_fband[ie],
                r_part_tot_fband_etot- r_part_tot_fband[ie]))/PI/(frequency1+frequency2));
  fclose(fw1);
  }
}
#endif
/******************************************************************************/
// store chosen function into XSPEC photar array
if(abs(photar_sw)<=10){
  if(frequency1 == 0.){
    for(ie=0;ie<ne;ie++){
      if(photar_sw == -1) photar[ie] = r_part_int[ie];
      else if(photar_sw == -2) photar[ie] = im_part_int[ie];
      else if(photar_sw == -3) photar[ie] = ampl_int[ie];
      else if(photar_sw == -4) photar[ie] = phase_int[ie];
      else if(photar_sw == -5){
        if(exclude_energy){
          utmp = flux_bands_prim[nbands-1] * r_part_int_etot -
                 far_prim[ie] * r_part_int[ie];
          utmp1 = flux_bands_prim[nbands-1] * im_part_int_etot -
                  far_prim[ie] * im_part_int[ie];
          photar[ie] = ampl_int[ie] / 
                       ( sqrt( utmp * utmp + utmp1 * utmp1 ) / 
                         ( flux_bands_prim[nbands-1] - far_prim[ie] ) );
        }else photar[ie] = ampl_int[ie] / ampl_int_etot;
      }else if(photar_sw == -6){
        if(exclude_energy){
          utmp = flux_bands_prim[nbands-1] * r_part_int_etot -
                 far_prim[ie] * r_part_int[ie];
          utmp1 = flux_bands_prim[nbands-1] * im_part_int_etot -
                  far_prim[ie] * im_part_int[ie];
          photar[ie] = - (phase_int[ie] - atan2(utmp1, utmp)) / PI / freq_wrap;
        }else photar[ie] = - (phase_int[ie] - phase_int_etot) / PI / freq_wrap;
      }else if(photar_sw == -7){
        if(exclude_energy) photar[ie] = ampl2_int[ie] / ampl2_int_etot2[ie];
        else photar[ie] = ampl2_int[ie] / ampl2_int_etot;
      }else if(photar_sw == -8){
        if(exclude_energy)photar[ie] = - (delay2_int[ie] - delay2_int_etot2[ie]);
        else photar[ie] = - (delay2_int[ie] - delay2_int_etot);
      }else if(photar_sw == 1) photar[ie] = r_part_tot_int[ie];
      else if(photar_sw == 2) photar[ie] = im_part_tot_int[ie];
      else if(photar_sw == 3) photar[ie] = ampl_tot_int[ie];
      else if(photar_sw == 4) photar[ie] = phase_tot_int[ie];
      else if(photar_sw == 5){
        if(exclude_energy){
          utmp=r_part_tot_int_etot -r_part_tot_int[ie];
          utmp1=im_part_tot_int_etot-im_part_tot_int[ie];
          photar[ie] = ampl_tot_int[ie] / sqrt( utmp * utmp + utmp1 * utmp1 );
        }else photar[ie] = ampl_tot_int[ie] / ampl_tot_int_etot;
      }else if(photar_sw == 6){
        if(exclude_energy){
          utmp=r_part_tot_int_etot -r_part_tot_int[ie];
          utmp1=im_part_tot_int_etot-im_part_tot_int[ie];
          photar[ie] = - (phase_tot_int[ie] - atan2(utmp1, utmp)) / PI / freq_wrap;
        }else photar[ie] = - (phase_tot_int[ie] - phase_tot_int_etot) / PI / freq_wrap;
      }else if(photar_sw == 7){
        if(exclude_energy) photar[ie] = ampl2_tot_int[ie] / ampl2_tot_int_etot2[ie];
        else photar[ie] = ampl2_tot_int[ie] / ampl2_tot_int_etot;
      }else if(photar_sw == 8){
        if(exclude_energy) photar[ie] = - (delay2_tot_int[ie] - delay2_tot_int_etot2[ie]);
        else photar[ie] = - (delay2_tot_int[ie] - delay2_tot_int_etot);
      }
      photar[ie] *= (ear[ie+1] - ear[ie]);
    }
  }else if(frequency1 >= frequency2){
    if(frequency1 > freq[nn/2]){
      sprintf(errortxt, "fft_reverberation: the frequency must be lower than or equal to %lg Hz.",
                        freq[nn/2]);
      xs_write(errortxt, 5);
      for (ie = 0; ie < ne; ie++) photar[ie] = 0.;
      goto error;
    }
//  given frequency, find the corresponding index in freq[]:
    ifr0 = (int) ceil(frequency1 / freq[1]);
    ttmp = (frequency1 - freq[ifr0 - 1]) / (freq[ifr0] - freq[ifr0 - 1]);
    ttmp1 = 1. - ttmp;
    for(ie=0;ie<ne;ie++){
      if(photar_sw == -1){
        y1 = r_part[ie+ne*(ifr0-1)];
        y2 = r_part[ie+ne*ifr0];
      }else if(photar_sw == -2){
        y1 = im_part[ie+ne*(ifr0-1)];
        y2 = im_part[ie+ne*ifr0];
      }else if(photar_sw == -3 || photar_sw == -5 || photar_sw == -7){
        y1 = ampl[ie+ne*(ifr0-1)];
        y2 = ampl[ie+ne*ifr0];
      }else if(photar_sw == -4 || photar_sw == -6 || photar_sw == -8){
        y1 = phase[ie+ne*(ifr0-1)];
        y2 = phase[ie+ne*ifr0];
      }else if(photar_sw == 1){
        y1 = r_part_tot[ie+ne*(ifr0-1)];
        y2 = r_part_tot[ie+ne*ifr0];
      }else if(photar_sw == 2){
        y1 = im_part_tot[ie+ne*(ifr0-1)];
        y2 = im_part_tot[ie+ne*ifr0];
      }else if(photar_sw == 3 || photar_sw == 5 || photar_sw == 7){
        y1 = ampl_tot[ie+ne*(ifr0-1)];
        y2 = ampl_tot[ie+ne*ifr0];
      }else if(photar_sw == 4 || photar_sw == 6 || photar_sw == 8){
        y1 = phase_tot[ie+ne*(ifr0-1)];
        y2 = phase_tot[ie+ne*ifr0];
      }
      photar[ie] = (ttmp1 * y1 + ttmp * y2);
      if( (photar_sw ==-5 || photar_sw ==-6 || photar_sw ==-7 || photar_sw ==-8) && nbands > 0){
        y1 = r_part_bands[nbands-1+nbands*(ifr0-1)];
        y2 = r_part_bands[nbands-1+nbands*ifr0];
        if(exclude_energy){
          y1 = ( flux_bands_prim[nbands-1] * y1 - 
                 far_prim[ie] * r_part[ie+ne*(ifr0-1)] ) / flux_bands_prim[nbands-1];
          y2 = ( flux_bands_prim[nbands-1] * y2 - 
                 far_prim[ie] * r_part[ie+ne*ifr0] ) / flux_bands_prim[nbands-1];
        }
        r_part_etot = (ttmp1 * y1 + ttmp * y2);
        y1 = im_part_bands[nbands-1+nbands*(ifr0-1)];
        y2 = im_part_bands[nbands-1+nbands*ifr0];
        if(exclude_energy){
          y1 = ( flux_bands_prim[nbands-1] * y1 - 
                 far_prim[ie] * im_part[ie+ne*(ifr0-1)] ) / flux_bands_prim[nbands-1];
          y2 = ( flux_bands_prim[nbands-1] * y2 - 
                 far_prim[ie] * im_part[ie+ne*ifr0] ) / flux_bands_prim[nbands-1];
        }
        im_part_etot = (ttmp1 * y1 + ttmp * y2);
        ampl_etot  = sqrt( r_part_etot  * r_part_etot +
                           im_part_etot * im_part_etot );
        phase_etot = atan2(im_part_etot, r_part_etot);
      }
      if( (photar_sw ==5 || photar_sw ==6 || photar_sw ==7 || photar_sw ==8) && nbands > 0){
        y1 = r_part_bands_tot[nbands-1+nbands*(ifr0-1)];
        y2 = r_part_bands_tot[nbands-1+nbands*ifr0];
        if(exclude_energy){
          y1 -= r_part_tot[ie+ne*(ifr0-1)];
          y2 -= r_part_tot[ie+ne*ifr0];
        }
        r_part_tot_etot=(ttmp1 * y1 + ttmp * y2);
        y1 = im_part_bands_tot[nbands-1+nbands*(ifr0-1)];
        y2 = im_part_bands_tot[nbands-1+nbands*ifr0];
        if(exclude_energy){
          y1 -= im_part_tot[ie+ne*(ifr0-1)];
          y2 -= im_part_tot[ie+ne*ifr0];
        }
        im_part_tot_etot=(ttmp1 * y1 + ttmp * y2);
        ampl_tot_etot  = sqrt( r_part_tot_etot  * r_part_tot_etot +
                               im_part_tot_etot * im_part_tot_etot );
        phase_tot_etot = atan2(im_part_tot_etot, r_part_tot_etot);
      }
      if( ( photar_sw == -5 || photar_sw == -7 ) && nbands > 0)photar[ie] /= ampl_etot;
      if( ( photar_sw == 5 || photar_sw == 7 ) && nbands > 0)photar[ie] /= ampl_tot_etot;
      if( photar_sw == -6 || photar_sw == -8 )
        photar[ie] = - (photar[ie] - ( nbands > 0 ? phase_etot : 0.) ) / (PI2 * frequency1);
      if( photar_sw == 6 || photar_sw == 8 )
        photar[ie] = - (photar[ie] - ( nbands > 0 ? phase_tot_etot : 0.) ) / (PI2 * frequency1);
      photar[ie] *= (ear[ie+1] - ear[ie]);
    }
  }else{
    for(ie=0;ie<ne;ie++){
      if(photar_sw == -1) photar[ie] = r_part_fband[ie];
      else if(photar_sw == -2) photar[ie] = im_part_fband[ie];
      else if(photar_sw == -3) photar[ie] = ampl_fband[ie];
      else if(photar_sw == -4) photar[ie] = phase_fband[ie];
      else if(photar_sw == -5){
        if(exclude_energy){
          utmp = flux_bands_prim[nbands-1] * r_part_fband_etot -
                 far_prim[ie] * r_part_fband[ie];
          utmp1 = flux_bands_prim[nbands-1] * im_part_fband_etot -
                  far_prim[ie] * im_part_fband[ie];
          photar[ie] = ampl_fband[ie] / 
                       ( sqrt( utmp * utmp + utmp1 * utmp1 ) / 
                         ( flux_bands_prim[nbands-1] - far_prim[ie] ) );
        }else photar[ie] = ampl_fband[ie] / ampl_fband_etot;
// we will use relative delay with respect to the whole energy band and redefine to be positive
      }else if(photar_sw == -6){
        if(exclude_energy){
          utmp = flux_bands_prim[nbands-1] * r_part_fband_etot -
                 far_prim[ie] * r_part_fband[ie];
          utmp1 = flux_bands_prim[nbands-1] * im_part_fband_etot -
                  far_prim[ie] * im_part_fband[ie];
          photar[ie] =  -(phase_fband[ie]-atan2(utmp1,utmp))/PI/(frequency2+frequency1);
        }else
          photar[ie] = -(phase_fband[ie]-phase_fband_etot)/PI/(frequency2+frequency1);
      }else if(photar_sw == -7){
        if(exclude_energy) photar[ie] = ampl2_fband[ie] / ampl2_fband_etot2[ie];
        else photar[ie] = ampl2_fband[ie] / ampl2_fband_etot;
// we will use relative delay with respect to the whole energy band and redefine to be positive
      }else if(photar_sw == -8){
        if(exclude_energy) photar[ie] = - ( delay2_fband[ie] - delay2_fband_etot2[ie] );
        else photar[ie] = - ( delay2_fband[ie] - delay2_fband_etot );
      }else if(photar_sw == 1) photar[ie] = r_part_tot_fband[ie];
      else if(photar_sw == 2) photar[ie] = im_part_tot_fband[ie];
      else if(photar_sw == 3) photar[ie] = ampl_tot_fband[ie];
      else if(photar_sw == 4)
        photar[ie] = phase_tot_fband[ie];
      else if(photar_sw == 5){
        if(exclude_energy){
          utmp=r_part_tot_fband_etot -r_part_tot_fband[ie];
          utmp1=im_part_tot_fband_etot-im_part_tot_fband[ie];
          photar[ie] = ampl_tot_fband[ie] / sqrt( utmp * utmp + utmp1 * utmp1 );
        }else photar[ie] = ampl_tot_fband[ie]/ampl_tot_fband_etot;
// we will use relative delay with respect to the whole energy band and redefine to be positive
      }else if(photar_sw == 6){
        if(exclude_energy){
          utmp=r_part_tot_fband_etot -r_part_tot_fband[ie];
          utmp1=im_part_tot_fband_etot-im_part_tot_fband[ie];
          photar[ie] =  -(phase_tot_fband[ie]-atan2(utmp1,utmp))/PI/(frequency2+frequency1);
        }else
          photar[ie] = -(phase_tot_fband[ie]-phase_tot_fband_etot)/PI/(frequency2+frequency1);
      }else if(photar_sw == 7){
        if(exclude_energy) photar[ie] = ampl2_tot_fband[ie]/ampl2_tot_fband_etot2[ie];
        else photar[ie] = ampl2_tot_fband[ie]/ampl2_tot_fband_etot;
// we will use relative delay with respect to the whole energy band and redefine to be positive
      }else if(photar_sw == 8){
        if(exclude_energy) photar[ie] =  -(delay2_tot_fband[ie]-delay2_tot_fband_etot2[ie]);
        else photar[ie] = -(delay2_tot_fband[ie]-delay2_tot_fband_etot);
      }
      photar[ie] *= (ear[ie+1] - ear[ie]);
    }
  }
}else{
//we do not need to multiply by the frequency bin width since we already did it!
  for(ie=0;ie<ne;ie++){
    if(photar_sw == -11) photar[ie] = r_part_bands_rebin[2*ie];
    else if(photar_sw == -12) photar[ie] = im_part_bands_rebin[2*ie];
    else if(photar_sw == -13) photar[ie] = ampl_bands_rebin[2*ie];
    else if(photar_sw == -14) 
      photar[ie] = phase_bands_rebin[2*ie] * (ear[ie+1] - ear[ie]);
    else if(photar_sw == -15){
      photar[ie] = ampl_bands_rebin[2*ie];
      if(nbands > 1)photar[ie] *= (ear[ie+1] - ear[ie]) / ampl_bands_rebin[1+2*ie];
    }else if(photar_sw == -16){
      photar[ie] = phase_bands_rebin[2*ie];
      if(nbands > 1){
        photar[ie] -= phase_bands_rebin[1+2*ie];
        if(photar[ie] < -PI)photar[ie]+=PI2;
        if(photar[ie] >  PI)photar[ie]-=PI2;
      }
      photar[ie] *=  (ear[ie+1] - ear[ie]) / (PI * (ear[ie]+ear[ie+1]));
    }else if(photar_sw == -17){
      photar[ie] = ampl2_bands_rebin[2*ie];
      if(nbands > 1)photar[ie] *= (ear[ie+1] - ear[ie]) / ampl2_bands_rebin[1+2*ie];
    }else if(photar_sw == -18){
      photar[ie] = delay2_bands_rebin[2*ie];
      if(nbands > 1)photar[ie] -= delay2_bands_rebin[1+2*ie];
    }else if(photar_sw == 11) photar[ie] = r_part_bands_tot_rebin[2*ie];
    else if(photar_sw == 12) photar[ie] = im_part_bands_tot_rebin[2*ie];
    else if(photar_sw == 13) photar[ie] = ampl_bands_tot_rebin[2*ie];
    else if(photar_sw == 14) 
      photar[ie] = phase_bands_tot_rebin[2*ie] * (ear[ie+1] - ear[ie]);
    else if(photar_sw == 15){
      photar[ie] = ampl_bands_tot_rebin[2*ie];
      if(nbands > 1)photar[ie] *= (ear[ie+1] - ear[ie]) / ampl_bands_tot_rebin[1+2*ie];
    }else if(photar_sw == 16){
      photar[ie] = phase_bands_tot_rebin[2*ie];
      if(nbands > 1){
        photar[ie] -= phase_bands_tot_rebin[1+2*ie];
        if(photar[ie] < -PI)photar[ie]+=PI2;
        if(photar[ie] >  PI)photar[ie]-=PI2;
      }
      photar[ie] *= (ear[ie+1] - ear[ie]) / (PI * (ear[ie]+ear[ie+1]));
    }else if(photar_sw == 17){
      photar[ie] = ampl2_bands_tot_rebin[2*ie];
      if(nbands > 1)photar[ie] *= (ear[ie+1] - ear[ie]) / ampl2_bands_tot_rebin[1+2*ie];
    }else if(photar_sw == 18){
      photar[ie] = delay2_bands_tot_rebin[2*ie];
      if(nbands > 1)photar[ie] -= delay2_bands_tot_rebin[1+2*ie];
    }
  }
}

error:
//Let's free all allocated arrays
if (data_r  != NULL){free((void *) data_r);  data_r  = NULL;}
if (data_im != NULL){free((void *) data_im); data_im = NULL;}
if (freq != NULL){free((void *) freq); freq = NULL;}
if (r_part  != NULL){free((void *) r_part);  r_part  = NULL;}
if (im_part != NULL){free((void *) im_part); im_part = NULL;}
if (ampl    != NULL){free((void *) ampl);    ampl    = NULL;}
if (phase   != NULL){free((void *) phase);   phase   = NULL;}
if (phase_u1!= NULL){free((void *) phase_u1);phase_u1 = NULL;}
if (phase_u2!= NULL){free((void *) phase_u2);phase_u2 = NULL;}
if (ampl_etot2  != NULL){free((void *) ampl_etot2);  ampl_etot2  = NULL;}
if (phase_etot2 != NULL){free((void *) phase_etot2); phase_etot2 = NULL;}
if (r_part_tot  != NULL){free((void *) r_part_tot);  r_part_tot  = NULL;}
if (im_part_tot != NULL){free((void *) im_part_tot); im_part_tot = NULL;}
if (ampl_tot    != NULL){free((void *) ampl_tot);    ampl_tot    = NULL;}
if (phase_tot   != NULL){free((void *) phase_tot);   phase_tot   = NULL;}
if (phase_tot_u1!= NULL){free((void *) phase_tot_u1);phase_tot_u1= NULL;}
if (phase_tot_u2!= NULL){free((void *) phase_tot_u2);phase_tot_u2= NULL;}
if (ampl_tot_etot2    != NULL){free((void *) ampl_tot_etot2);    ampl_tot_etot2    = NULL;}
if (phase_tot_etot2   != NULL){free((void *) phase_tot_etot2);   phase_tot_etot2   = NULL;}
if (r_part_bands  != NULL){free((void *) r_part_bands);  r_part_bands  = NULL;}
if (im_part_bands != NULL){free((void *) im_part_bands); im_part_bands = NULL;}
if (ampl_bands    != NULL){free((void *) ampl_bands);    ampl_bands    = NULL;}
if (phase_bands   != NULL){free((void *) phase_bands);   phase_bands   = NULL;}
if (phase_bands_u1!= NULL){free((void *) phase_bands_u1);phase_bands_u1= NULL;}
if (r_part_bands_tot  != NULL){free((void *) r_part_bands_tot);  r_part_bands_tot  = NULL;}
if (im_part_bands_tot != NULL){free((void *) im_part_bands_tot); im_part_bands_tot = NULL;}
if (ampl_bands_tot    != NULL){free((void *) ampl_bands_tot);    ampl_bands_tot    = NULL;}
if (phase_bands_tot   != NULL){free((void *) phase_bands_tot);   phase_bands_tot   = NULL;}
if (phase_bands_tot_u1!= NULL){free((void *) phase_bands_tot_u1);phase_bands_tot_u1= NULL;}

return;
}
/*******************************************************************************
 ******************************************************************************/
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four(double *data_r, double *data_im, long n1, int isign)
{
unsigned long nn,n,mmax,m,j,istep,i;
double wtemp,wr,wpr,wpi,wi,theta;
double tempr,tempi;
double *data;

nn=(unsigned long) n1;
data=(double *)malloc((size_t) ((2 * n1 + 1)*sizeof(double)));
for (i=0;i<n1;i++) {
 data[2*i+1]=data_r[i];
 data[2*i+2]=data_im[i];
}
n=nn << 1;
j=1;
for (i=1;i<n;i+=2) {
  if (j > i) {
    SWAP(data[j],data[i]);
    SWAP(data[j+1],data[i+1]);
  }
  m=nn;
  while (m >= 2 && j > m) {
    j -= m;
    m >>= 1;
  }
  j += m;
}
mmax=2;
while (n > mmax) {
  istep=mmax << 1;
  theta=isign*(6.28318530717959/mmax);
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0;
  wi=0.0;
  for (m=1;m<mmax;m+=2) {
    for (i=m;i<=n;i+=istep) {
      j=i+mmax;
      tempr=wr*data[j]-wi*data[j+1];
      tempi=wr*data[j+1]+wi*data[j];
      data[j]=data[i]-tempr;
      data[j+1]=data[i+1]-tempi;
      data[i] += tempr;
      data[i+1] += tempi;
    }
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  mmax=istep;
}
for (i=0;i<n1;i++) {
  data_r[i]=data[2*i+1]/(isign == -1 ? n1 : 1.);
  data_im[i]=data[2*i+2]/(isign == -1 ? n1 : 1.);
}
free((char*) (data)); 
}
#undef SWAP
