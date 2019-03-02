#ifndef FILTER_H
#define FILTER_H
#pragma once
#include "stdafx.h"

#define MAX_NUM_FILTER_TAPS 1000

enum filterType { LPF, HPF, BPF };//generic filter class can be set to be a high pass, low pass, or bandpass filter

class FIRFilter {//an FIR filter is defined by the sum y[n] = b_0*x[n] + ... + b_N*x[n-N], where y[n] = output, x[n] = input value
private:
	filterType m_filt_t;
	int m_num_taps;	//# of taps is the # of coefficients in the filter's output sum
	int m_error_flag;
	double m_Fs;
	double m_Fx;
	double m_lambda;
	double *m_taps;
	double *m_sr;
	void designLPF();
	void designHPF();

	// Only needed for the bandpass filter case
	double m_Fu, m_phi;
	void designBPF();

public:
	FIRFilter(filterType filt_t, int num_taps, double Fs, double Fx);
	FIRFilter(filterType filt_t, int num_taps, double Fs, double Fl, double Fu);
	~FIRFilter();
	void init();
	double do_sample(double data_sample);
	int get_error_flag() { return m_error_flag; };
	void get_taps(double *taps);
	int write_taps_to_file(char* filename);
	int write_freqres_to_file(char* filename);
};



class OnePole {//a one-pole filter
public:
	OnePole();
	OnePole(double Fc);
	~OnePole();
	void setFc(double Fc);
	float process(float in);
protected:
	double a0, b1, z1;
};
#endif
