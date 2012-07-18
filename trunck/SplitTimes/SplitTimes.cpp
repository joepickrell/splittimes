/*
 * SplitTimes.cpp
 *
 *  Created on: Jul 12, 2012
 *      Author: pickrell
 */





#include "Kimura.h"
#include "CountData.h"
#include "CmdLine.h"
using namespace std;


void printopts(){
        cout << "Options:\n";
        cout << "-i [file name] input file\n";
        cout << "-p population\n";
        cout << "-d drift length\n";
        cout << "-o [file name] outfile file\n";
        cout << "\n";
}

string infile;
string outfile = "split_out";
int blocksize = 500;
string population;
double d;

int main(int argc, char *argv[]){

    CCmdLine cmdline;
    if (cmdline.SplitLine(argc, argv) < 1){
        printopts();
        exit(1);
    }
    if (cmdline.HasSwitch("-i")) infile = cmdline.GetArgument("-i", 0).c_str();
    else{
    	printopts();
    	exit(1);
    }
    if (cmdline.HasSwitch("-p")) population = cmdline.GetArgument("-p", 0).c_str();
    else{
    	printopts();
    	exit(1);
    }
    if (cmdline.HasSwitch("-d")) d= atof(cmdline.GetArgument("-d", 0).c_str());
    else{
    	printopts();
    	exit(1);
    }
    if (cmdline.HasSwitch("-o")) outfile = cmdline.GetArgument("-o", 0).c_str();

    CountData cdata(infile, blocksize);
    int which = cdata.pop2id[population];
    vector<int> counts = cdata.get_der_counts(which);
    for (vector<int>::iterator it = counts.begin(); it != counts.end(); it++){
    	cout << *it << "\n";
    }
    /*
    vector<int> counts;
    ifstream in(infile.c_str());
    vector<string> line;
    string st;
    while(getline(in, st)){
            string buf;
            stringstream ss(st);
            line.clear();
            while (ss>> buf){
                    line.push_back(buf);
            }
            counts.push_back(atoi(line[0].c_str()));
    }
*/
    Kimura k(d);
    int m = counts.size()-1;
    k.optim_anc(m, counts);
    k.print_spec_compare(outfile+".spec", m, counts);
    k.print_params(outfile+".params");
    k.reset();
    k.optim_anc_h0(m, counts);
    k.print_spec_compare(outfile+".spec0", m, counts);
    double lambda = k.lambda;
    vector<double> lambdas;
    for (int i = 0; i < cdata.nblock; i++){
    	vector<int> tmpc = cdata.get_der_counts_jackknife(which, i);
    	//for (vector<int>::iterator it = tmpc.begin(); it != tmpc.end(); it++){
    	//	cout << *it << "\n";
    	//}
    	//k.reset();
    	k.optim_anc(m, tmpc);
    	cout << k.lambda << "\n";
    	lambdas.push_back(k.lambda);
    }
    double m_lambda = 0;
    for (vector<double>::iterator it = lambdas.begin(); it != lambdas.end(); it++) m_lambda+= *it;
    m_lambda = m_lambda/ (double) lambdas.size();
    double sum = 0;
    for (vector<double>::iterator it = lambdas.begin(); it != lambdas.end(); it++) sum += (*it - m_lambda)*(*it- m_lambda);
    double se = sqrt(  ((double) cdata.nblock - 1)/ (double) cdata.nblock * sum);
    string tmpoutf = outfile +".lambdase";
    ofstream lout (tmpoutf.c_str());
    lout << m_lambda << " "<< se << "\n";



    return 0;
}
