#include "BooPHF.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <random>


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



typedef boomphf::HashFunctors<u_int64_t>  hasher_t;
typedef boomphf::mphf<  u_int64_t, hasher_t  > boophf_t;



u_int64_t nelem = 100000000;
int nthreads = 1;
int main (int argc, char* argv[])
{
	
	
	if(argc >=2 )
	{
		nelem = strtoul(argv[1], NULL,0);
	}
	if(argc ==3 )
	{
		nthreads = atoi(argv[2]);
	}
	
	printf("Construct a BooPHF with  %lli elements (%lli MB for holding elems in ram) \n",nelem,nelem*sizeof(u_int64_t)/1024LL/1024LL);

	
	
	//creating some random input data
//	srandomdev(); //init random generator
//	data = (u_int64_t * ) calloc(nelem,sizeof(u_int64_t));
//	for (u_int64_t i = 0; i < nelem; i++)
//		data[i] = random64();
	
	
	static std::mt19937_64 rng;
	rng.seed(std::mt19937_64::default_seed);
	data = (u_int64_t * ) calloc(nelem,sizeof(u_int64_t));
	for (u_int64_t i = 1; i < nelem; i++)
		data[i] = rng();
	
	///create the boophf
	
	auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(data), static_cast<const u_int64_t*>(data+nelem));

	clock_t begin, end;
	begin = clock();
	
	boophf_t * bphf = new boomphf::mphf<u_int64_t,hasher_t>(nelem,data_iterator,2.4);
	
	end = clock();

	printf("BooPHF constructed perfect hash for %llu keys in %.2fs\n", nelem, (double)(end - begin) / CLOCKS_PER_SEC);

	u_int64_t mphf_value;
	
	printf("boophf  bits/elem : %f\n",(float) (bphf->totalBitSize())/nelem);

	
#ifdef CHECK_MPHF
	//test the mphf
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
				//printf("collision for val %lli \n",mphf_value);
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

	}
	
	
	free(check_table);
#endif
	
	free(data);
	return EXIT_SUCCESS;
}

