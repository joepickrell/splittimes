/*
 * test.cpp
 *
 *  Created on: Jul 11, 2012
 *      Author: pickrell
 */

#include "Kimura.h"
#include <cstdlib>
#include <iostream>

int main(){
	Kimura t(0.03);
	//t.print_transition_probs("transprobs");
	vector<double> s = t.get_anc_spectrum();

	vector<double> n = t.get_spectrum(s);
	//for (int i = 0; i < n.size(); i++){
	//	cout << i <<  " "<< s[i]<< " "<< n[i] << "\n";
	//}
	//for (int i = 1 ; i <= 1000; i++){
	//	cout << t.freqs[i] << " "<< s[i] << "\n";
	//}
	//double lk = t.llk(10, 3, s);
	//cout << lk << "\n";
	//vector<int> tmp;
	//for (int i = 0; i <=10; i++){
	//	int toadd = 0;
	//	if (i == 3) toadd = 40;
	//	tmp.push_back(toadd);
	//}
	//double lk2 = t.llk_anc(10, tmp, 0.0001);

	//double lk2 = t.all_llk(10, tmp, s);
	//cout << lk2 << "\n";

	//vector<double> tmpspec = t.get_all_spectrum(0.1, 0.01);
	//for (int i = 0; i <= t.N; i++){
	//		cout << t.freqs[i] << " "<< s[i] << " "<< tmpspec[i]<< "\n";
	//}
	int m = 10;
	vector<int> spec;
	spec.push_back(36092); spec.push_back(24056); spec.push_back(19682); spec.push_back(15889); spec.push_back(12696);
	spec.push_back(10614);spec.push_back(8588); spec.push_back(6976); spec.push_back(5782);spec.push_back(4305); spec.push_back(3168 );
	//t.print_spec_compare(m, spec);
	t.optim_anc(m, spec);
	t.print_spec_compare(m, spec);
	t.print_params("params");
	//cout << t.llk_anc(m, spec)<< "\n";
	//for (int i = 18; i < 19; i++){
	//	double f = (double) i /1000.0;
	//	t.lambda = f;
	//	t.optim_anc(m, spec);
	//	cout << f << " "<< t.llk_anc(m, spec)<< "\n";
	//}
	//t.print_spec_compare(m, spec);
	//cout << t.a << " "<< t.b << " "<< t.c << "\n";
	//t.optim_anc(m, spec);
	//t.print_spec_compare(m, spec);
	//cout << t.a << " "<< t.b << " "<< t.c << "\n";
	//cout << t.llk_anc(m, spec) << "\n";
	//vector<double> now_spec = t.get_spectrum(s);
	//for (int i = 0; i <= t.N+1; i++){
	//		cout << t.freqs[i] << " "<< s[i] << " "<< now_spec[i]<< "\n";
	//}

	return 0;
}


