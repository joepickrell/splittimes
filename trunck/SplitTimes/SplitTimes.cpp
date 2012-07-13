/*
 * SplitTimes.cpp
 *
 *  Created on: Jul 12, 2012
 *      Author: pickrell
 */





#include "Kimura.h"
#include "CmdLine.h"
using namespace std;


void printopts(){
        cout << "Options:\n";
        cout << "-i [file name] input file\n";
        cout << "-d drift length\n";
        cout << "-o [file name] outfile file\n";
        cout << "\n";
}

string infile;
string outfile = "split_out";
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
    if (cmdline.HasSwitch("-d")) d= atof(cmdline.GetArgument("-d", 0).c_str());
    else{
    	printopts();
    	exit(1);
    }
    if (cmdline.HasSwitch("-o")) outfile = cmdline.GetArgument("-o", 0).c_str();

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

    Kimura k(d);
    int m = counts.size()-1;
    k.optim_anc(m, counts);
    //cout << "done?\n"; cout.flush();
    k.print_spec_compare(outfile+".spec", m, counts);
    k.print_params(outfile+".params");

    return 0;
}
