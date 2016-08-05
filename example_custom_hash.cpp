#include "BooPHF.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <random>
#include <algorithm>

using namespace std;



//example with user provided custom hasher for uint64_t type :

class Custom_uint64_Hasher
{
public:
	// the class should have operator () with this signature :
	// BBhash will use the 'seed' paramater to generate two different hash values form this key.
     //then it will generate internally a sequence of hash values using xorshifts, using these two first hash values as starting point.
	uint64_t operator ()   (uint64_t key, uint64_t seed=0) const
	{
		
		key ^= key >> 33;
		key *= 0xff51afd7ed558ccd;
		key ^= key >> 33;
		key *= 0xc4ceb9fe1a85ec53;
		key ^= key >> 33;
		
		key ^= seed;
		
		return key;
	}
};


//then tell BBhash to use this custom hash : (also appears below, line 104)
typedef boomphf::mphf<  u_int64_t, Custom_uint64_Hasher  > boophf_t;


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
	
	double gammaFactor = 2.0; // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
	// gamma = 2 is a good tradeoff (leads to approx 3.7 bits/key )

	//build the mphf
	bphf = new boomphf::mphf<u_int64_t,Custom_uint64_Hasher>(nelem,data_iterator,nthreads,gammaFactor);
	
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
