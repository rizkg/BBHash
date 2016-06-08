#include "BooPHF.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <random>
#include <algorithm>

using namespace std;

typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
typedef boomphf::mphf<  u_int64_t, hasher_t  > boophf_t;

int main (int argc, char* argv[]){
	
	//PARAMETERS
	u_int64_t nelem = 1000000;
	uint nthreads = 1;

	if(argc !=3 ){
		printf("Usage :\n");
		printf("%s <nelem> <nthreads> \n",argv[0]);
		return EXIT_FAILURE;
	}
	
	if(argc ==3 ){
		nelem = strtoul(argv[1], NULL,0);
		nthreads = atoi(argv[2]);
	}
	
	uint64_t ii, jj;
	u_int64_t *data;

	/////  generation of random keys
	uint64_t rab = 100;
	static std::mt19937_64 rng;
	rng.seed(std::mt19937_64::default_seed); //default seed
	
	//rng.seed(seed2); //random seed from timer
	data = (u_int64_t * ) calloc(nelem+rab,sizeof(u_int64_t));
	
	for (u_int64_t i = 1; i < nelem+rab; i++){
		data[i] = rng();
	}
	printf("de-duplicating items \n");
	
	std::sort(data,data+nelem+rab);
	
	for (ii = 1, jj = 0; ii < nelem+rab; ii++) {
		if (data[ii] != data[jj])
			data[++jj] = data[ii];
	}
	printf("found %lli duplicated items  \n",nelem+rab-(jj + 1) );
	
	//////////////////
	// at this point, array data contains a set of nelem random unique keys
	
	
	boophf_t * bphf = NULL;
	double t_begin,t_end; struct timeval timet;
	
	
	printf("Construct a BooPHF with  %lli elements  \n",nelem);
	
	gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
	
	// mphf takes as input a c++ range. A simple array of keys can be wrapped with boomphf::range
	// but could be from a user defined iterator (enabling keys to be read from a file or from some complex non-contiguous structure)
	auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(data), static_cast<const u_int64_t*>(data+nelem));
	
	double gammaFactor = 1.0; // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query

	//build the mphf
	bphf = new boomphf::mphf<u_int64_t,hasher_t>(nelem,data_iterator,nthreads,gammaFactor);
	
	gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
	double elapsed = t_end - t_begin;
	
	
	printf("BooPHF constructed perfect hash for %llu keys in %.2fs\n", nelem,elapsed);
	printf("boophf  bits/elem : %f\n",(float) (bphf->totalBitSize())/nelem);
	
	//query mphf like this
	uint64_t  idx = bphf->lookup(data[0]);
	printf(" example query  %lli ----->  %llu \n",data[0],idx);
	
	free(data);
	delete bphf;
	return EXIT_SUCCESS;
}
