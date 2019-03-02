#ifndef WAVELET_H
#define WAVELET_H
#include "stdafx.h"

// 1D Functions
void* dwt1(std::string, std::vector<double> &, std::vector<double> &, std::vector<double> &);
void* dyadic_zpad_1d(std::vector<double> &);
double convol(std::vector<double> &, std::vector<double> &, std::vector<double> &);
int filtcoef(std::string, std::vector<double> &, std::vector<double> &, std::vector<double> &,std::vector<double> &);
void downsamp(std::vector<double> &, int, std::vector<double> &);
void upsamp(std::vector<double> &, int, std::vector<double> &);
void circshift(std::vector<double> &, int);
int sign(int);
void* idwt1(std::string wname, std::vector<double> &, std::vector<double> &, std::vector<double> &);
int vecsum(std::vector<double> &, std::vector<double> &, std::vector<double> &);

// 1D Symmetric Extension DWT Functions
void* dwt_sym(std::vector<double> &, int, std::string, std::vector<double> &, std::vector<double> &,	std::vector<int> &);
void* dwt1_sym(std::string, std::vector<double> &, std::vector<double> &, std::vector<double> &);
void* idwt_sym(std::vector<double> &, std::vector<double> &, std::string, std::vector<double> &, std::vector<int> &);
void* symm_ext(std::vector<double> &, int);
void* idwt1_sym(std::string, std::vector<double> &, std::vector<double> &, std::vector<double> &); // Not Tested

// 1D Stationary Wavelet Transform
std::vector<double> swt(std::vector<double> signal1, int J, std::string filterCoefficientName, std::vector<double> swt_output, int length);
void* iswt(std::vector<double> &, int, std::string, std::vector<double> &);
void* per_ext(std::vector<double> &, int);

// 2D Functions
void* branch_lp_dn(std::string, std::vector<double> &, std::vector<double> &);
void* branch_hp_dn(std::string, std::vector<double> &, std::vector<double> &);
void* branch_lp_hp_up(std::string, std::vector<double> &, std::vector<double> &, std::vector<double> &);
//void* dwt_2d(std::vector<std::vector<double> > &, int , std::string , std::vector<std::vector<double> > &, std::vector<double> &) ;
//void* idwt_2d(std::vector<std::vector<double> > &,std::vector<double> &, std::string ,std::vector<std::vector<double> > &);
void* dyadic_zpad_2d(std::vector<std::vector<double> > &, std::vector<std::vector<double> > &);
void* dwt_output_dim(std::vector<std::vector<double> >&, int &, int &);
void* zero_remove(std::vector<std::vector<double> > &, std::vector<std::vector<double> > &);
void* getcoeff2d(std::vector<std::vector<double> > &, std::vector<std::vector<double> > &,std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, std::vector<double> &, int &);
void* idwt2(std::string, std::vector<std::vector<double> > &, std::vector<std::vector<double> >  &,std::vector<std::vector<double> >  &, std::vector<std::vector<double> >  &, std::vector<std::vector<double> > &);
void* dwt2(std::string, std::vector<std::vector<double> > &, std::vector<std::vector<double> >  &,std::vector<std::vector<double> >  &, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &);
void* downsamp2(std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, int, int);
void* upsamp2(std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, int, int);

// 2D DWT (Symmetric Extension) Functions
void* dwt_2d_sym(std::vector<std::vector<double> > &, int, std::string, std::vector<double> &, std::vector<double> &,	std::vector<int> &);
void* dwt2_sym(std::string, std::vector<std::vector<double> > &, std::vector<std::vector<double> >  &,std::vector<std::vector<double> >  &, std::vector<std::vector<double> >  &, std::vector<std::vector<double> > &);
void* idwt_2d_sym(std::vector<double>  &, std::vector<double> &, std::string, std::vector<std::vector<double> > &,std::vector<int> &);
void* circshift2d(std::vector<std::vector<double> > &, int, int);
void symm_ext2d(std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, int);
void* dispDWT(std::vector<double> &, std::vector<std::vector<double> > &, std::vector<int> &, std::vector<int> &, int);
void* dwt_output_dim_sym(std::vector<int> &, std::vector<int> &, int);

//2D Stationary Wavelet Transform
void* swt_2d(std::vector<std::vector<double> > &, int, std::string, std::vector<double> &);
void* per_ext2d(std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, int);

// FFT functions
double convfft(std::vector<double> &, std::vector<double> &, std::vector<double> &);
double convfftm(std::vector<double> &, std::vector<double> &, std::vector<double> &);
void* fft(std::vector<std::complex<double> > &, int, unsigned int);
void* bitreverse(std::vector<std::complex<double> > &);
void* freq(std::vector<double> &, std::vector<double> &);

//New
void* dwt1_sym_m(std::string wname, std::vector<double> &signal, std::vector<double> &cA, std::vector<double> &cD);//FFTW3 for 2D
void* idwt1_sym_m(std::string wname, std::vector<double> &X, std::vector<double> &app, std::vector<double> &detail);
void* dwt(std::vector<double> &sig, int J, std::string nm, std::vector<double> &dwt_output, std::vector<double> &flag, std::vector<int> &length);
void* idwt(std::vector<double> &, std::vector<double> &, std::string, std::vector<double> &, std::vector<int> &);
void* dwt_2d(std::vector<std::vector<double> > &, int, std::string, std::vector<double> &, std::vector<double> &,	std::vector<int> &);
void* dwt1_m(std::string wname, std::vector<double> &signal, std::vector<double> &cA, std::vector<double> &cD);
void* idwt_2d(std::vector<double>  &dwtop, std::vector<double> &flag, std::string nm,std::vector<std::vector<double> > &idwt_output, std::vector<int> &length);
void* idwt1_m(std::string wname, std::vector<double> &X, std::vector<double> &cA, std::vector<double> &cD);
void* dwt_output_dim2(std::vector<int> &length, std::vector<int> &length2, int J);
#endif