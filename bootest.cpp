#include "BooPHF.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <random>
#include <algorithm> 
#include <climits>
u_int64_t *data;

//uncomment to check correctness of the func
//#define CHECK_MPHF

#define MAX_RANDOM 2147483648
#define srandomdev() srand((unsigned) time(NULL))

uint64_t random64 ()
{
	uint64_t low, high,res;
	low = random();
	high = random();
	
	res = (high << 32) + low;
	return res;
}



typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
typedef boomphf::mphf<  u_int64_t, hasher_t  > boophf_t;



u_int64_t nelem = 100000000;
int nthreads = 1;
int main (int argc, char* argv[])
{
	bool check_correctness = false;
	
	bool bench_lookup = false;
	
	
	if(argc >=2 )
	{
		nelem = strtoul(argv[1], NULL,0);
	}
	
	if(argc >=3 )
		nthreads = atoi(argv[2]);

	
	for (int ii=3; ii<argc; ii++)
	{
		if(!strcmp("-check",argv[ii])) check_correctness= true;
		if(!strcmp("-bench",argv[ii])) bench_lookup= true;
	}

	
	
	printf("Construct a BooPHF with  %lli elements (%lli MB for holding elems in ram) \n",nelem,nelem*sizeof(u_int64_t)/1024LL/1024LL);
	
	
	//creating some random input data
	//	srandomdev(); //init random generator
	//	data = (u_int64_t * ) calloc(nelem,sizeof(u_int64_t));
	//	for (u_int64_t i = 0; i < nelem; i++)
	//		data[i] = random64();
	
	uint64_t rab = 100;
//	
	static std::mt19937_64 rng;
	rng.seed(std::mt19937_64::default_seed);
	data = (u_int64_t * ) calloc(nelem+rab,sizeof(u_int64_t));

	for (u_int64_t i = 1; i < nelem+rab; i++)
		data[i] = rng();
//	
	
	
	
	//methode simple pas besoin de dup
//	data = (u_int64_t * ) calloc(nelem,sizeof(u_int64_t));
//		data[0] = 0;
//		u_int64_t step = ULLONG_MAX / nelem;
//			for (u_int64_t i = 1; i < nelem; i++)
//				data[i] = data[i-1] +step;
	
//	
	uint64_t ii, jj;
	printf("de-duplicating items \n");
	
	std::sort(data,data+nelem+rab);

	for (ii = 1, jj = 0; ii < nelem+rab; ii++) {
		if (data[ii] != data[jj])
			data[++jj] = data[ii];
	}

	printf("found %lli duplicated items  \n",nelem+rab-(jj + 1) );

	
	
	
	

	
	///create the boophf
	
	auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(data), static_cast<const u_int64_t*>(data+nelem));
	
	clock_t begin, end;
	
	double t_begin,t_end; struct timeval timet;
	
	gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
	
	bool fastmode = true; //build with fast mode, requires a little bit more ram ( 3 % of elems  are loaded in ram)
	boophf_t * bphf = new boomphf::mphf<u_int64_t,hasher_t>(nelem,data_iterator,nthreads,fastmode,2.0);
	
	gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);

	double elapsed = t_end - t_begin;

	
	
	printf("BooPHF constructed perfect hash for %llu keys in %.2fs\n", nelem,elapsed);
	
	u_int64_t mphf_value;
	
	printf("boophf  bits/elem : %f\n",(float) (bphf->totalBitSize())/nelem);
	
	
	if(check_correctness)	//test the mphf
	{
		u_int64_t nb_collision_detected = 0;
		u_int64_t range_problems = 0;
		char * check_table = (char * ) calloc(nelem,sizeof(char));
		
		
		for (u_int64_t i = 0; i < nelem; i++)
		{
			mphf_value = bphf->lookup(data[i]);
			if(mphf_value>=nelem)
			{
				range_problems++;
			}
			if(check_table[mphf_value]==0)
			{
				check_table[mphf_value]=1;
			}
			else
			{
				printf("collision for val %lli \n",mphf_value);
				nb_collision_detected++;
			}
		}
		
		if(nb_collision_detected ==  0 && range_problems ==0)
		{
			printf(" --- boophf working correctly --- \n");
		}
		else
		{
			printf("!!! problem, %llu collisions detected; %llu out of range !!!\n",nb_collision_detected,range_problems);
			return EXIT_FAILURE;

			
		}
		
		
		free(check_table);
	}
	
	
	
	if(bench_lookup)
	{
		u_int64_t dumb=0;
		begin = clock();
		for (u_int64_t i = 0; i < nelem; i++)
		{
			mphf_value = bphf->lookup(data[i]);
			
			//do some silly work
			dumb+= mphf_value;
			
		}
		
		end = clock();
		printf("BooPHF %llu lookups in  %.2fs,  approx  %.2f ns per lookup   (fingerprint %llu)  \n", nelem, (double)(end - begin) / CLOCKS_PER_SEC,  ((double)(end - begin) / CLOCKS_PER_SEC)*1000000000/nelem,dumb);
	}
	
	
	free(data);
	return EXIT_SUCCESS;
}


