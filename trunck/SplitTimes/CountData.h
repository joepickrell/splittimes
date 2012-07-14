/*
 * CountData.h
 *
 *  Created on: Apr 1, 2011
 *      Author: pickrell
 */

#ifndef COUNTDATA_H_
#define COUNTDATA_H_
#include "Settings.hpp"
using namespace std;

class CountData{
public:

	CountData(string, int);
	void read_counts(string);
	vector<int> get_der_counts(int);
	vector<int> get_der_counts_jackknife(int, int);
	int K, nblock;
	map<string, int> pop2id;
	map<int, string> id2pop;
	map<int, double> mean_hzy;
	map<int, double> mean_ninds;

	vector<vector<pair<int, int> > > allele_counts;

	int npop, nsnp;
	string get_pops(); //in Newick format

	//vector<int> get_der_counts(int);

};


#endif /* COUNTDATA_H_ */
