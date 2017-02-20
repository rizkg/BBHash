#include "BooPHF.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <random>
#include <algorithm>
#include <fstream>

using namespace std;



//example with user provided custom hasher for uint64_t type :

class Custom_string_Hasher
{
public:
	// the class should have operator () with this signature :
	uint64_t operator ()   (std::string key, uint64_t seed=0) const
	{


		uint64_t hash  =  hash_fn(key);
		hash ^= seed;

		return hash;
	}

     std::hash<std::string> hash_fn;
};


//then tell BBhash to use this custom hash : (also appears below, line 104)
typedef boomphf::mphf<  std::string, Custom_string_Hasher  > boophf_t;





int main (int argc, char* argv[]){

	//PARAMETERS
	u_int64_t nelem = 0;
	uint nthreads = 1;
	uint gamma=2;

	if(argc !=4 ){
		printf("Usage :\n");
		printf("%s <file> <gamma> <nthreads> \n",argv[0]);
		return EXIT_FAILURE;
	}
	ifstream in;
	if(argc ==4 ){
		in.open(argv[1],std::ifstream::in);
		gamma = atoi(argv[2]);
		nthreads = atoi(argv[3]);
	}

	uint64_t ii, jj;

	boophf_t * bphf = NULL;
	double t_begin,t_end; struct timeval timet;


	printf("Construct a BooPHF with  %lli elements  \n",nelem);

	gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);

	std::vector<std::string> data;
	std::string line;
	while(not in.eof()){
		getline(in, line);
		data.push_back(line);
		nelem++;
	}


	//build the mphf
	bphf = new boomphf::mphf<std::string,Custom_string_Hasher>(nelem,data,nthreads,gamma);

	gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
	double elapsed = t_end - t_begin;


	printf("BooPHF constructed perfect hash for %llu keys in %.2fs\n", nelem,elapsed);
	printf("boophf  bits/elem : %f\n",(float) (bphf->totalBitSize())/nelem);
	gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);

	//query mphf like this
	for (u_int64_t i = 0; i < nelem; i++){
		uint64_t  idx = bphf->lookup(data[i]);
	}
	gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
	double elapsed2 = t_end - t_begin;

	printf("Query of %llu key  in %.2fs\n", nelem,elapsed2);

	delete bphf;
	return EXIT_SUCCESS;
}
