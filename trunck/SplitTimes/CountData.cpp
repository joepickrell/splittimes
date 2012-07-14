/*
 * CountData.cpp
 *
 *  Created on: Apr 1, 2011
 *      Author: pickrell
 */
#include "CountData.h"


CountData::CountData(string infile, int blocksize){
	K = blocksize;
	read_counts(infile);
	cout << nsnp << " SNPs in "<< nblock << " blocks\n";
}

vector<int> CountData::get_der_counts(int which){
	int nchrom = allele_counts[0][which].first +allele_counts[0][which].second+1;
	vector<int> toreturn(nchrom, 0);
	for (vector<vector<pair<int, int> > >::iterator it = allele_counts.begin(); it != allele_counts.end(); it++){
		int count = it->at(which).first;
		toreturn[count]+=1;
	}
	return toreturn;
}

vector<int> CountData::get_der_counts_jackknife(int which, int which_jack){
	int nchrom = allele_counts[0][which].first +allele_counts[0][which].second+1;
	vector<int> toreturn(nchrom, 0);
	for (int i = 0; i < nsnp; i++){
		if ( i >= which_jack*K && i < (which_jack+1)*K) continue;
		int count = allele_counts[i][which].first;
		toreturn[count]+=1;
	}
	return toreturn;
}

void CountData::read_counts(string infile){
    allele_counts.clear();
    pop2id.clear();
    id2pop.clear();
    npop = 0;
    nsnp = 0;

    string ext = infile.substr(infile.size()-3, 3);
    if (ext != ".gz"){
    	std::cerr << infile << " is not gzipped (only .gz files accepted)\n";
    	exit(1);
    }
	igzstream in(infile.c_str()); //only gzipped files
    vector<string> line;
    struct stat stFileInfo;
    int intStat;
    string st, buf;

    intStat = stat(infile.c_str(), &stFileInfo);
    if (intStat !=0){
            std::cerr<< "ERROR: cannot open file " << infile << "\n";
            exit(1);
    }

    /*
     * header contains population names
     */
    getline(in, st);
    stringstream ss(st);
    line.clear();
    while (ss>> buf){
    	line.push_back(buf);
     }
    /*
     * make map from header, number populations according to order
     */
    int start = 0;

    for(int i = start; i < line.size(); i++) {
    	pop2id.insert(make_pair(line[i], i-start));
    	id2pop.insert(make_pair(i-start, line[i]));
    	npop ++;
    }
    int headsize = line.size();
    /*
     * read counts, store in allele_counts
     */
    while(getline(in, st)){
            buf.clear();
            stringstream ss(st);
            line.clear();
            while (ss>> buf){
                    line.push_back(buf);
            }
            vector<pair<int, int> > topush;

            if (line.size() != headsize){
            	cerr << "ERROR: Line "<< nsnp <<" has "<< line.size() << " entries. Header has "<< headsize <<"\n";
            	exit(1);
            }
            for ( int i = start; i < line.size(); i++){
            	//cout <<  line[i] << "\n";
                typedef boost::tokenizer<boost::char_separator<char> >
                tokenizer;
                boost::char_separator<char> sep(",");
                tokenizer tokens(line[i], sep);
                vector<int> tmpcounts;
                for (tokenizer::iterator tok_iter = tokens.begin();  tok_iter != tokens.end(); ++tok_iter){
                        int tmp = atoi(tok_iter->c_str());
                        tmpcounts.push_back(tmp);
                }
                if (tmpcounts.size() != 2){
                	std::cerr << "ERROR: "<< line[i] << " does not have two alleles (expecting SNP data)\n";
                	exit(1);
                }
                topush.push_back(make_pair(tmpcounts[0], tmpcounts[1]));
            }
            allele_counts.push_back(topush);
            nsnp++;
    }
    nblock = nsnp/K;
}

/*

set<pair<string, pair<double, double> > > CountData::calculate_f4(int i0, int i1, int i2, int i3){
	set<pair<string, pair<double, double> > > toreturn;
	double mean1, se1, mean2, se2, mean3, se3;
	vector<double> f4_1;
	vector<double> f4_2;
	vector<double> f4_3;

	vector<double> f4_block_1;
	vector<double> f4_block_2;
	vector<double> f4_block_3;
	//calculate the covariance matrix in each block
	cout << "Estimating f_4 in "<< nblock << " blocks of size "<< params->window_size <<"\n"; cout.flush();
	mean1 = 0;
	mean2 = 0;
	mean3 =0;
	int total_nsnp = 0;
	for (int i = 0; i < nsnp; i++){
		if (isnan(gsl_matrix_get(alfreqs, i, i0))) continue;
		if (isnan(gsl_matrix_get(alfreqs, i, i1))) continue;
		if (isnan(gsl_matrix_get(alfreqs, i, i2))) continue;
		if (isnan(gsl_matrix_get(alfreqs, i, i3))) continue;
		double toadd = (gsl_matrix_get(alfreqs, i, i1) - gsl_matrix_get(alfreqs, i, i0))*(gsl_matrix_get(alfreqs, i, i3) - gsl_matrix_get(alfreqs, i, i2) );
		f4_1.push_back(toadd);
		mean1 += toadd;

		toadd = (gsl_matrix_get(alfreqs, i, i2) - gsl_matrix_get(alfreqs, i, i0))*(gsl_matrix_get(alfreqs, i, i3) - gsl_matrix_get(alfreqs, i, i1) );
		f4_2.push_back(toadd);
		mean2 += toadd;

		toadd = (gsl_matrix_get(alfreqs, i, i3) - gsl_matrix_get(alfreqs, i, i0))*(gsl_matrix_get(alfreqs, i, i2) - gsl_matrix_get(alfreqs, i, i1) );
		f4_3.push_back(toadd);
		mean3 += toadd;
		total_nsnp ++;
	}

	cout << "total_nsnp "<< total_nsnp << " nsnp "<< nsnp << "\n";
	mean1 = mean1 / (double) total_nsnp;
	mean2 = mean2 / (double) total_nsnp;
	mean3 = mean3 / (double) total_nsnp;
	for (int i = 0; i < nblock ; i++){
		double c1 = 0;
		double c2 = 0;
		double c3 = 0;
		int tmp_nsnp = 0;
		for (int n = 0; n < nsnp; n++){

			if (isnan(gsl_matrix_get(alfreqs, n, i0))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, i1))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, i2))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, i3))) continue;
			if ( n >= i*params->window_size && n < (i+1)*params->window_size) continue;
			double toadd = (gsl_matrix_get(alfreqs, n, i1) - gsl_matrix_get(alfreqs, n, i0))*(gsl_matrix_get(alfreqs, n, i3) - gsl_matrix_get(alfreqs, n, i2) );
			c1+= toadd;

			toadd = (gsl_matrix_get(alfreqs, n, i2) - gsl_matrix_get(alfreqs, n, i0))*(gsl_matrix_get(alfreqs, n, i3) - gsl_matrix_get(alfreqs, n, i1) );
			c2+= toadd;

			toadd = (gsl_matrix_get(alfreqs, n, i3) - gsl_matrix_get(alfreqs, n, i0))*(gsl_matrix_get(alfreqs, n, i2) - gsl_matrix_get(alfreqs, n, i1) );
			c3+= toadd;
			tmp_nsnp++;
		}
		double cov1 = c1/ (double) tmp_nsnp;
		double cov2 = c2/ (double) tmp_nsnp;
		double cov3 = c3/ (double) tmp_nsnp;
		//cout << i <<  " "<< cov1 << " "<< cov2 << " "<< cov3 << "\n";
		f4_block_1.push_back(cov1);
		f4_block_2.push_back(cov2);
		f4_block_3.push_back(cov3);

	}

	// and standard error
	double sum1 = 0;
	double sum2 = 0;
	double sum3 = 0;

	for (vector<double>::iterator it = f4_block_1.begin(); it != f4_block_1.end(); it++) sum1  += *it;
	mean1 = sum1/ (double) nblock;

	for (vector<double>::iterator it = f4_block_2.begin(); it != f4_block_2.end(); it++) sum2  += *it;
	mean2 = sum2/ (double) nblock;

	for (vector<double>::iterator it = f4_block_3.begin(); it != f4_block_3.end(); it++) sum3  += *it;
	mean3 = sum3/ (double) nblock;

	sum1 = 0;
	sum2 = 0;
	sum3 = 0;
	for (vector<double>::iterator it = f4_block_1.begin(); it != f4_block_1.end(); it++) {
		//cout << mean1 << " "<< *it << "\n";
		sum1+= (*it-mean1)*(*it-mean1);
	}
	se1 = ( (double) nblock- 1.0) / (double) nblock  * sum1;
	se1= sqrt(se1);

	for (vector<double>::iterator it = f4_block_2.begin(); it != f4_block_2.end(); it++) {
		//cout << mean2 << " "<< *it << "\n";
		sum2+= (*it-mean2)*(*it-mean2);
	}
	se2 = ( (double) nblock- 1.0) / (double) nblock  * sum2;
	se2 = sqrt(se2);

	for (vector<double>::iterator it = f4_block_3.begin(); it != f4_block_3.end(); it++) sum3+= (*it-mean3)*(*it-mean3);
	se3 = ( (double) nblock- 1.0) / (double) nblock  * sum3;
	se3 = sqrt(se3);

	pair<double, double> tmp = make_pair(mean1, se1);
	pair<string, pair<double, double> > tmp2 = make_pair( id2pop[i0] +","+id2pop[i1]+";"+id2pop[i2]+","+id2pop[i3], tmp );
	toreturn.insert(tmp2);

	tmp = make_pair(mean2, se2);
	tmp2 = make_pair( id2pop[i0] +","+id2pop[i2]+";"+id2pop[i1]+","+id2pop[i3], tmp );
	toreturn.insert(tmp2);

	tmp = make_pair(mean3, se3);
	tmp2 = make_pair( id2pop[i0] +","+id2pop[i3]+";"+id2pop[i1]+","+id2pop[i2], tmp );
	toreturn.insert(tmp2);


	return toreturn;
}

void CountData::set_cov_jackknife(int which){
	gsl_matrix_set_zero(cov);
	for(int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			string p1 = id2pop[i];
			string p2 = id2pop[j];
			double m = 0;
			for (int k = 0; k < nblock; k++){
				if (k == which) continue;
				m+= cov_samp[p1][p2].at(k);
			}
			m = m/ (double) (nblock-1);
			gsl_matrix_set(cov, i, j, m);
			gsl_matrix_set(cov, j, i, m);
		}
	}

}
*/

