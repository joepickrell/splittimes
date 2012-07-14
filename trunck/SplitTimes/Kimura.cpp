/*
 * Kimura.cpp
 *
 *  Created on: Jul 11, 2012
 *      Author: pickrell
 */
#include "Kimura.h"
using namespace std;

Kimura::Kimura(double drift){
	N = 100;
	freqs.push_back(0);
	for (int i = 0; i < N; i++){
		double tmpf1 = (double) i / (double) N;
		double tmpf2 = (double) (i+1) / (double) N;
		double tmpf3 = (tmpf1+tmpf2)/2;
		freqs.push_back(tmpf3);
	}
	freqs.push_back(1);
	tau = drift;
	set_transition_probs();
	a = 0.0;
	b = -1.0;
	c = 1;
	epsilon = 1e-7;
	lambda = 0.001;
    phi = (1+sqrt(5))/2;
    resphi = 2-phi;

}

double Kimura::get_poly_dens(double y, double x, double t){
	// density * 1000;
	double toreturn = evalkim(y, x, t, 1.0) ;
	return toreturn;
}


vector<double> Kimura::get_all_spectrum(double x, double t){
	vector<double> toreturn(N+2, 0.0);
	double sum = 0;
	double m = 0;
	for (int i = 1; i <= N ; i++){

		double f = freqs[i];
		//cout << f << "\n";
		//f = (double)i / (double) N;
		double d = get_poly_dens(f, x, t);
		sum += d;
		m += d*f;
		toreturn[i] += d;
	}

	double fixup = x*(double)N - m;
	double fixdown =  (double) N - sum - fixup;
	//cout << x << " "<< m << " "<< sum << " "<< fixup << " "<< fixdown << "\n";
	toreturn[0] += fixdown;
	toreturn[N+1] += fixup;
	return toreturn;
}

vector<double> Kimura::get_spectrum(vector<double> st){
	vector<double> toreturn(N+2, 0.0);

	// integrate
	for (int i = 0; i <= N; i++){
		//cout << i  << "\n"; cout.flush();
		double f = freqs[i];
		double std = st[i]/(double)N;
		vector<double> tmp = transition_probs[i];
		for (int j = 0; j <= N+1; j++){
			toreturn[j] += tmp[j]/(double)N* std;
		}
	}
	//toreturn[0] = st[0]/ (double)N;
	return toreturn;
}

vector<double> Kimura::get_anc_spectrum(){
	vector<double> toreturn(N+2, 0.0);
	double sum = 0;
	for (int i = 1 ; i <= N; i++){
		double f = freqs[i];
		double add  = a*f*f + b*f +c;
		sum += add;
		toreturn[i] = add;
	}
	for (int i = 0 ; i <= N+1; i++){
		toreturn[i] = (1-lambda)*(toreturn[i] / sum)* (double)N;
	}
	toreturn[0] = lambda*(double)N;
	return toreturn;
}

double Kimura::llk(int n, int k, vector<double> spec){
	// integrate binomial
	double toreturn = 0;
	for (int i = 0; i<=N+1; i++){
		double add = 0;
		double f = freqs[i];
		double d = spec[i];
		add+= d* gsl_sf_fact(n)/(gsl_sf_fact(k)*gsl_sf_fact(n-k)) * pow(f, k)*pow(1-f, n-k);
		toreturn += add;
	}
	//if (k < 1e-8) toreturn = lambda + (1-lambda)*toreturn/ (double)N;
	toreturn = toreturn/(double)N;
	return log(toreturn);
}

double Kimura::all_llk(int m, vector<int> obs, vector<double> spec){
	double toreturn = 0;
	for (int i = 0; i <= m; i++){
		int count = obs[i];
		int k = i;
		double single = llk(m, k, spec);
		//lambda = 0.3;
		//double testsingle = llk(m, k, spec);
		//lambda = 0.0001;
		//cout << i << " "<< single << " "<< testsingle << "\n";
		double toadd = count*single;
		toreturn += toadd;
	}
	return toreturn;
}
void Kimura::reset(){
	a = 0.0;
	b = -1.0;
	c = 1.0;
	lambda = 0.001;
}
double Kimura::llk_anc(int m, vector<int> obs){
	double toreturn;

	vector<double> anc_spec = get_anc_spectrum();

	vector<double> now_spec = get_spectrum(anc_spec);
	toreturn = all_llk(m, obs, now_spec);
	return toreturn;
}

void Kimura::set_transition_probs(){
	transition_probs.clear();
	for (int i = 0; i <= N+1; i++){
		//cout << i  << "\n";
		double f = freqs[i];
		vector<double> tmpprobs = get_all_spectrum(f, tau);
		transition_probs.push_back(tmpprobs);
	}
	cout << "\n";

}

void Kimura::print_transition_probs(string outfile){
	ofstream out(outfile.c_str());
	for (int i = 0; i < transition_probs.size(); i++){
		for(int j = 0; j < transition_probs.size(); j++) out<< transition_probs[i][j] << " ";
		out <<  "\n";
	}
}

double Kimura::optim_anc(int m, vector<int> obs){
	double start_llik = llk_anc(m, obs);
	//cout << "start "<< start_llik << "\n";
	double current_llk = start_llik;
	bool done = false;
	int nit = 0;
	while(!done){

		double guessa = a;
		double guessb = b;
		double guessc = c;
		double guesslambda = log(lambda);
		double min = -10.0;
		double max = 10.0;
		//cout << "guessl "<< guesslambda << "\n";
		golden_section_a(min, guessa, max, 0.0001, m, obs, &current_llk);
		//cout << "c_a "<< current_llk << "\n";
		golden_section_b(min, guessb, max, 0.0001, m, obs, &current_llk);
		//cout << "c_b "<< current_llk << "\n";
		golden_section_c(min, guessc, max, 0.0001, m, obs, &current_llk);
		//cout << "c_c "<< current_llk << "\n";
		golden_section_lambda(-20.0, guesslambda, 0, 0.0001, m, obs, &current_llk);
		double new_llik = llk_anc(m, obs);
		if (new_llik < start_llik+ epsilon) done = true;
		else start_llik = new_llik;

	}
	return current_llk;
}


double Kimura::optim_anc_h0(int m, vector<int> obs){
	double start_llik = llk_anc(m, obs);
	//cout << "start "<< start_llik << "\n";
	double current_llk = start_llik;
	bool done = false;
	int nit = 0;
	lambda = 0;
	while(!done){

		double guessa = a;
		double guessb = b;
		double guessc = c;
		//double guesslambda = log(lambda);
		double min = -10.0;
		double max = 10.0;
		//cout << "guessl "<< guesslambda << "\n";
		golden_section_a(min, guessa, max, 0.0001, m, obs, &current_llk);
		//cout << "c_a "<< current_llk << "\n";
		golden_section_b(min, guessb, max, 0.0001, m, obs, &current_llk);
		//cout << "c_b "<< current_llk << "\n";
		golden_section_c(min, guessc, max, 0.0001, m, obs, &current_llk);
		//cout << "c_c "<< current_llk << "\n";
		//golden_section_lambda(-20.0, guesslambda, -1, 0.0001, m, obs, &current_llk);
		double new_llik = llk_anc(m, obs);
		if (new_llik < start_llik+ epsilon) done = true;
		else start_llik = new_llik;

	}
	return current_llk;
}

int Kimura::golden_section_a(double min, double guess, double max, double tau, int m, vector<int> obs, double * current_llk){
	double x;

	if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
	else x = guess - resphi *(guess-min);
	if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
		double new_a = (min+max)/2;
		a = new_a;
		*current_llk = llk_anc(m, obs);
		return 0;
	}

	a = x;
	double f_x = -llk_anc(m, obs);

	a = guess;
	double f_guess = -llk_anc(m, obs);

	if (f_x < f_guess){
		if ( (max-guess) > (guess-min) )	return golden_section_a(guess, x, max, tau, m, obs, current_llk);
		else return golden_section_a(min, x, guess, tau, m, obs, current_llk);
	}
	else{
		if ( (max - guess) > (guess - min)  ) return golden_section_a(min, guess, x, tau, m, obs, current_llk);
		else return golden_section_a(x, guess, max, tau, m, obs, current_llk);
	}
}

int Kimura::golden_section_b(double min, double guess, double max, double tau, int m, vector<int> obs, double * current_llk){
	double x;

	if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
	else x = guess - resphi *(guess-min);
	if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
		double new_b = (min+max)/2;
		b = new_b;
		*current_llk = llk_anc(m, obs);
		return 0;
	}

	b = x;
	double f_x = -llk_anc(m, obs);

	b = guess;
	double f_guess = -llk_anc(m, obs);

	if (f_x < f_guess){
		if ( (max-guess) > (guess-min) )	return golden_section_b(guess, x, max, tau, m, obs, current_llk);
		else return golden_section_b(min, x, guess, tau, m, obs, current_llk);
	}
	else{
		if ( (max - guess) > (guess - min)  ) return golden_section_b(min, guess, x, tau, m, obs, current_llk);
		else return golden_section_b(x, guess, max, tau, m, obs, current_llk);
	}
}


int Kimura::golden_section_c(double min, double guess, double max, double tau, int m, vector<int> obs, double * current_llk){
	double x;

	if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
	else x = guess - resphi *(guess-min);
	if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
		double new_c = (min+max)/2;
		c = new_c;
		*current_llk = llk_anc(m, obs);
		return 0;
	}

	c = x;
	double f_x = -llk_anc(m, obs);

	c = guess;
	double f_guess = -llk_anc(m, obs);

	if (f_x < f_guess){
		if ( (max-guess) > (guess-min) )	return golden_section_c(guess, x, max, tau, m, obs, current_llk);
		else return golden_section_c(min, x, guess, tau, m, obs, current_llk);
	}
	else{
		if ( (max - guess) > (guess - min)  ) return golden_section_c(min, guess, x, tau, m, obs, current_llk);
		else return golden_section_c(x, guess, max, tau, m, obs, current_llk);
	}
}


int Kimura::golden_section_lambda(double min, double guess, double max, double tau, int m, vector<int> obs, double * current_llk){
	double x;

	//cout << "setting x "<< guess << " "<< min << " "<< max << "\n"; cout.flush();
	if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
	else x = guess - resphi *(guess-min);
	//cout << "x set to "<< x << "\n"; cout.flush();
	if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
		double new_lambda = (min+max)/2;
		lambda = exp(new_lambda);
		*current_llk = llk_anc(m, obs);
		//cout << "returning\n";
		return 0;
	}
	//cout << "x before "<< x << "\n"; cout.flush();
	lambda = exp(x);
	double f_x = -llk_anc(m, obs);
	//cout << "x "<< lambda << " "<< f_x << "\n";  cout.flush();

	//cout << "guess before "<< guess << "\n"; cout.flush();
	lambda = exp(guess);

	double f_guess = -llk_anc(m, obs);
	//cout << "guess "<< lambda << " "<< f_guess << "\n"; cout.flush();

	if (f_x < f_guess){
		if ( (max-guess) > (guess-min) )	return golden_section_lambda(guess, x, max, tau, m, obs, current_llk);
		else return golden_section_lambda(min, x, guess, tau, m, obs, current_llk);
	}
	else{
		if ( (max - guess) > (guess - min)  ) return golden_section_lambda(min, guess, x, tau, m, obs, current_llk);
		else return golden_section_lambda(x, guess, max, tau, m, obs, current_llk);
	}
}

void Kimura::print_spec_compare(int m, vector<int> emp_spec){
	int total = 0;
	vector<double> theory_spec(m+1, 0.0);
	for (vector<int>::iterator it = emp_spec.begin(); it != emp_spec.end(); it++)total += *it;
	vector<double> anc_spec = get_anc_spectrum();
	for (int i = 0; i <= m ;i++){
		for (int j = 0; j < freqs.size(); j++){
			double initf = freqs[j];
			double initdens = anc_spec[j]/ (double) N;
			double tmpsum = 0;
			for (int k = 0; k < freqs.size(); k++){
				double f = freqs[k];
				double trans = transition_probs[j][k]/(double) N;
				double binom_prob = gsl_sf_fact(m)/(gsl_sf_fact(i)*gsl_sf_fact(m-i)) * pow(f, i)*pow(1-f, m-i);
				//cout << i << " "<< initf<< " "<< f << " "<< initdens << " "<< trans << " "<< binom_prob << "\n";
				tmpsum+= trans*binom_prob;

			}
			theory_spec[i] += initdens*tmpsum;
		}
	}
	//for (int i = 0; i <=m; i++){
	//	double tmp = theory_spec[i];
	//	if (i <1e-8) theory_spec[i] = lambda + (1-lambda) * tmp;
	//	else theory_spec[i] = (1-lambda) *tmp;
	//}
	for (int i = 0; i <=m; i++){
		cout << i << " "<< theory_spec[i]<< " "<<  (double) emp_spec[i]/(double) total << "\n";
	}
}


void Kimura::print_spec_compare(string f, int m, vector<int> emp_spec){
	ofstream outf(f.c_str());
	int total = 0;
	vector<double> theory_spec(m+1, 0.0);
	for (vector<int>::iterator it = emp_spec.begin(); it != emp_spec.end(); it++)total += *it;
	vector<double> anc_spec = get_anc_spectrum();
	for (int i = 0; i <= m ;i++){
		for (int j = 0; j < freqs.size(); j++){
			double initf = freqs[j];
			double initdens = anc_spec[j]/ (double) N;
			double tmpsum = 0;
			for (int k = 0; k < freqs.size(); k++){
				double f = freqs[k];
				double trans = transition_probs[j][k]/(double) N;
				double binom_prob = gsl_sf_fact(m)/(gsl_sf_fact(i)*gsl_sf_fact(m-i)) * pow(f, i)*pow(1-f, m-i);
				tmpsum+= trans*binom_prob;

			}
			theory_spec[i] += initdens*tmpsum;
		}
	}
	for (int i = 0; i <=m; i++){
		outf << i << " "<< theory_spec[i]<< " "<<  (double) emp_spec[i]/(double) total << "\n";
	}
}

void Kimura::print_params(string f){
	ofstream outf(f.c_str());
	outf << "a "<< a << "\n";
	outf << "b "<< b << "\n";
	outf << "c "<< c << "\n";
	outf << "lambda "<< lambda << "\n";
}
