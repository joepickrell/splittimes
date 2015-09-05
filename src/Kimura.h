/*
 * Kimura.h
 *
 *  Created on: Jul 11, 2012
 *      Author: pickrell
 */

#ifndef KIMURA_H_
#define KIMURA_H_

extern "C"{
	#include "kimsubs.h"
}
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <gsl/gsl_sf.h>
#include <cmath>
#include <sstream>
//#include "nicksrc/nicklib.h"
using namespace std;

class Kimura{
public:
	Kimura(double);
	int N;
	vector<double> freqs;
	double tau;
	vector<vector<double> > transition_probs;
	void set_transition_probs();
	void print_transition_probs(string);
	double a, b, c, lambda;
	double epsilon;
	double phi, resphi;
	double get_poly_dens(double, double, double);
	vector<double> get_all_spectrum(double, double);
	vector<double> get_spectrum(vector<double>);
	vector<double> get_anc_spectrum();
	double llk(int, int, vector<double>);
	double all_llk(int, vector<int>, vector<double>);
	double llk_anc(int, vector<int>);
	double optim_anc(int, vector<int>, bool);
	double optim_anc_h0(int, vector<int>);
	int golden_section_a(double, double, double, double, int, vector<int>, double *);
	int golden_section_b(double, double, double, double, int, vector<int>, double *);
	int golden_section_c(double, double, double, double, int, vector<int>, double *);
	int golden_section_lambda(double, double, double, double, int, vector<int>, double *);
	void print_spec_compare(int, vector<int>);
	void print_spec_compare(string, int, vector<int>);
	void print_params(string);
	void reset();
	//pair<vector<double>, vector<double> > get_complete_spectrum(vector<double>, double);
};


#endif /* KIMURA_H_ */
