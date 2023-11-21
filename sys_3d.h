#ifndef SYS_H
#define SYS_H

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <string>

#define N_THREADS 272

/*************** contants *****************/

const double OMEGAM = 0.305147;
const double OMEGAB = 0.0482266;
const double HUBBLE_h = 0.68;
const double FBARYON = OMEGAB/OMEGAM;
const double thompson = 6.28e-18; //cm2

/*************** Simulation Parameters *******************/

//Keep it as const as data type in header files
//0.256
const double BOXSIZE = 1.024; //L in Mpc/h

//direction vectors, 4 pixel connectivity
const int dx[] = {1,-1, 0, 0, 0, 0};
const int dy[] = {0, 0, 1,-1, 0, 0};
const int dz[] = {0, 0, 0, 0, 1,-1};

//PDF red/blue pops declaration
const double ratio_max_thres = 1.5;

/*********************************************************/

struct units_block {
  double X_UNIT;
  double TIME_UNIT;
  double VEL_UNIT;
  double RHO_UNIT;
  double H_UNIT;
  double He_UNIT;
  double EDEN_UNIT;
  double TEMPMIN;
};

int  find_systems(float*, int *, double, const long, const long, const long);
void find_total(float *data, int *labels, double *total, int comp, double dV, const long);
void find_size(int *labels, double *, int, double, const long);

void find_max(float *, int *, double *, int, const long);
void find_argmax(float *data, int *labels, double *max, int comp, int *x, int *y, int *z, const long NGRID, const long NGRIDR, const long NSLICE);
void find_min(float *, int *, double *, int, const long);

//get weightes average quantity over the skewers for each systems, that includes NHI and Delta_gas
void skewers (float *delta, float *temp, float *nHI, float *xHI, float *gamma, double *delta_w, double *temp_w, double *NHI, double *xHI_delta_w, double *gamma_w, int *labels, int comp, int sk, int *x, int *y, int *z, double dx, const long NGRID, const long NSLICE);

void systems(std::string, std::string, std::string, double, const long, const long);
void dfs(int, int, int, int, bool *, int *, const long, const long);
double alphaB(double);

/****************** Physical Constants *******************/
const double SIGMA_ALPHA[2] = {4.45e-18, 1.12e-18}; //cm^-2;  Lyman alpha cross sections of HI and HeII
const double CRITDENSITY = 1.8791e-29;   // Current critical density, in g/cm^3, not including h^2!
const double CRITDENMSOLMPC = 2.7755e11; // Current critical density, in Msol/Mpc^3, not including h^2!
const double MPROTON = 1.6726e-24;           // in grams
const double KBOLTZ = 1.38066e-16;            // in erg/K
const double MYRTOSEC = 3.15569e13;
const double KPCTOCM = 3.08560e21;           // cm/kpc
const double MPCTOCM = 3.08560e24;           // cm/Mpc
const double SOLARTOGRAMS = 1.989e33;
const double TCMB = 2.726;                   // CMB temperature at z=0
const double MUB = 1.22;           // Mean molecular weight, for primordial
const double MELECTRON = 9.11e-28; // in [grams]
const double EELECTRON = 1.60217662e-19; // Couloumbs
const double LIGHTSPEED = 3.0e5;  //km s^-1
const double LIGHTSPEED_CGS = 3.0e10;  //cm s^-1
const double mH_CGS = 1.67223e-24; //mass of Hydrogen in grams
const double mHe_CGS = 4*mH_CGS;
const double Tcmb0 = 2.7255; //Present day CMB temp in K
const double GAMMA_IDEAL = 5.0/3.0;
const double EION_HEII = 54.4; //ionization potential of HeII in eV
const double SIGMA_HEII = 6.3e-18/4.*1.21; //cm^2 -- photoionization cross section for HeII
const double SIGMA_HI = 6.3463e-18; //cm^2 -- photoionization cross section for HI
const double SIGMA_T = 6.6524e-25; //cm^2 -- thomson scattering constant
const double EVTOERG = 1.60217646e-12;
const double H0_CGS   = 3.24086e-18; //in 1/s.  Does not include h.
const double K2eV = KBOLTZ/EVTOERG;
const double eV2K = 1.0/K2eV;
const double HI_eV = 13.6;
const double HI_K = HI_eV*eV2K;
const double Planck_h = 6.626075540e-27; //Planck's constant in erg/s
#endif 
