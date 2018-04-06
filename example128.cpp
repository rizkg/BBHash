#include "BooPHF.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <random>
#include <algorithm>

using namespace std;

  class SingleHashFunctor128
{
	typedef std::pair<uint64_t,uint64_t> hash_pair_t;
	
public:
	hash_pair_t operator ()  (const __uint128_t& key) const  {
		hash_pair_t result;
		result.first =  singleHasher ((uint64_t) (key >> 64), 0xAAAAAAAA55555555ULL)  ^  singleHasher ((uint64_t)key, 0xAAAAAAAA55555555ULL) ;
		result.second =  singleHasher ((uint64_t) (key >> 64), 0x33333333CCCCCCCCULL)  ^  singleHasher ((uint64_t)key, 0x33333333CCCCCCCCULL) ;
		
		return result;
	}
	
	uint64_t singleHasher(const uint64_t & key, uint64_t seed) const
	{
		
		uint64_t hash = seed;
		hash ^= (hash <<  7) ^  key * (hash >> 3) ^ (~((hash << 11) + (key ^ (hash >> 5))));
		hash = (~hash) + (hash << 21);
		hash = hash ^ (hash >> 24);
		hash = (hash + (hash << 3)) + (hash << 8);
		hash = hash ^ (hash >> 14);
		hash = (hash + (hash << 2)) + (hash << 4);
		hash = hash ^ (hash >> 28);
		hash = hash + (hash << 31);
		return hash;
	}
};


typedef SingleHashFunctor128 hasher_t;

typedef boomphf::mphf<  hasher_t  > boophf_t;

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
	u_int64_t *data64;
	__uint128_t *data;

	/////  generation of random keys
	static std::mt19937_64 rng;
	rng.seed(std::mt19937_64::default_seed); //default seed
	
	//rng.seed(seed2); //random seed from timer
	data64 = (u_int64_t * ) calloc(2*nelem,sizeof(u_int64_t));
	
	for (u_int64_t i = 0; i < 2*nelem; i++){ //generating 64 bits at a time
		data64[i] = rng();
	}
	data = (__uint128_t *)data64;
	//printf("%llx %llx \n",(uint64_t)(data[0]>>64),(uint64_t)data[0] );


	
	boophf_t * bphf = NULL;
	double t_begin,t_end; struct timeval timet;
	
	
	printf("Construct a BooPHF with  %lli elements  \n",nelem);
	
	gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
	
	// mphf takes as input a c++ range. A simple array of keys can be wrapped with boomphf::range
	// but could be from a user defined iterator (enabling keys to be read from a file or from some complex non-contiguous structure)
	// if your input data is in a std::vector, it is already a c++ range, you can provide it std::vector directly

	auto data_iterator = boomphf::range(static_cast<const __uint128_t*>(data), static_cast<const __uint128_t*>(data+nelem));
	
	double gammaFactor = 2.0; // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query

	//build the mphf
	bphf = new boomphf::mphf<hasher_t>(nelem,data_iterator,nthreads,gammaFactor);// gammaFactor param is optional, will be 2 by default if not provided
	
	gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
	double elapsed = t_end - t_begin;
	
	
	printf("BooPHF constructed perfect hash for %llu keys in %.2fs\n", nelem,elapsed);
	printf("boophf  bits/elem : %f\n",(float) (bphf->mem_size()*8)/nelem);
	
	//query mphf like this
	uint64_t  idx = bphf->lookup(data[0]);
	printf("example query %llx%llx   ----->  %llu\n",(uint64_t)(data[0]>>64),(uint64_t)data[0],idx );


	free(data);
	delete bphf;
	return EXIT_SUCCESS;
}
