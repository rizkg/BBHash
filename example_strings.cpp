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
    typedef std::string Item;
    typedef std::pair<uint64_t,uint64_t> hash_pair_t;

    public:
    hash_pair_t operator ()  (const Item& key) const  {  
        hash_pair_t result;
        result.first =  singleHasher (key, 0xAAAAAAAA55555555ULL);;
        result.second =  singleHasher (key, 0x33333333CCCCCCCCULL);;

        return result;
    }

	uint64_t singleHasher (std::string key, uint64_t seed=0) const
	{
		

		uint64_t hash  =  hash_fn(key);
		hash ^= seed;
		
		return hash;
	}
	
     std::hash<std::string> hash_fn;
};


//then tell BBhash to use this custom hash : (also appears below, line 104)
typedef boomphf::mphf<  Custom_string_Hasher  > boophf_t;



//from http://stackoverflow.com/questions/440133/how-do-i-create-a-random-alpha-numeric-string-in-c
void gen_random(char *s, const int len) {
	static const char alphanum[] =
	"0123456789"
	"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	"abcdefghijklmnopqrstuvwxyz";
	
	for (int i = 0; i < len; ++i) {
		s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
	}
	
	s[len] = 0;
}



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
	std::vector<std::string> data;

	int string_size = 100;
	/////  generation of random strings

	char * tempchar = (char *) malloc(string_size*sizeof(char));

	//string lolstr;
	//ifstream inputfile("StringFile.txt",ios::in);
	for (u_int64_t i = 0; i < nelem; i++){
		//RANDOM STRINGS
		 gen_random(tempchar,string_size);
		 data.push_back((string)tempchar);
		
		//STRING READ FROM FILE
		//getline(inputfile,lolstr);
		//data.push_back(lolstr);
	}
	
	
	//////////////////
	// at this point, array data contains a set of nelem random unique keys
	
	boophf_t * bphf = NULL;
	double t_begin,t_end; struct timeval timet;
	
	
	printf("Construct a BooPHF with  %lli elements  \n",nelem);
	
	gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
	
	// mphf takes as input a c++ range. A std::vector is already a c++ range
	
	double gammaFactor = 2.0; // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
	// gamma = 2 is a good tradeoff (leads to approx 3.7 bits/key )

	//build the mphf
	bphf = new boomphf::mphf<Custom_string_Hasher>(nelem,data,nthreads,gammaFactor);
	
	gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
	double elapsed = t_end - t_begin;
	
	
	printf("BooPHF constructed perfect hash for %llu keys in %.2fs\n", nelem,elapsed);
	printf("boophf  bits/elem : %f\n",(float) (bphf->mem_size()*8)/nelem);
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
