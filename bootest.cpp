#include "BooPHF.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <random>
#include <algorithm>
#include <climits>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <thread>
#include <math.h>

//#include <chrono>


u_int64_t *data;


using namespace std;

//uncomment to check correctness of the func
//#define CHECK_MPHF


#define MAX_RANDOM 2147483648
#define srandomdev() srand((unsigned) time(NULL))


uint64_t korenXor(uint64_t x){
	x ^= (x << 21);
	x ^= (x >> 35);
	x ^= (x << 4);
	return x;
}


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


// iterator from disk file of u_int64_t  todo template
class bfile_iterator : public std::iterator<std::forward_iterator_tag, u_int64_t>{
	public:
	bfile_iterator()
	: _is(nullptr)
	, _pos(0)
	{}

	bfile_iterator(FILE* is)
	: _is(is)
	, _pos(0)
	{
		int reso = fseek(_is,0,SEEK_SET);
		advance();
	}
	u_int64_t const& operator*()  {  return _elem;  }

	bfile_iterator& operator++()
	{
		advance();
		return *this;
	}

	friend bool operator==(bfile_iterator const& lhs, bfile_iterator const& rhs)
	{
		if (!lhs._is || !rhs._is)  {  if (!lhs._is && !rhs._is) {  return true; } else {  return false;  } }
		assert(lhs._is == rhs._is);
		return rhs._pos == lhs._pos;
	}

	friend bool operator!=(bfile_iterator const& lhs, bfile_iterator const& rhs)  {  return !(lhs == rhs);  }
	private:
	void advance()
	{
		assert(_is);
		_pos++;
		int res = fread(&_elem,sizeof(u_int64_t),1,_is);
		if(res == 0)
		{
			_is = nullptr;
			_pos = 0;
			return;
		}
	}
	u_int64_t _elem;
	FILE * _is;
	unsigned long _pos;
};


class file_binary{
	public:

	file_binary(const char* filename)
	{
		_is = fopen(filename, "rb");
		if (!_is) {
			throw std::invalid_argument("Error opening " + std::string(filename));
		}
	}

	~file_binary()
	{
		fclose(_is);
	}

	bfile_iterator begin() const
	{
		return bfile_iterator(_is);
	}

	bfile_iterator end() const {return bfile_iterator(); }

	size_t        size () const  {  return 0;  }//todo ?

	private:
	FILE * _is;
};

//PARAMETERS
u_int64_t nelem = 1000*1000;
uint nthreads = 8; //warning must be a divisor of nBuckets
double gammaFactor = 1.0;
uint nBuckets = 96;
vector<FILE*> vFiles(nBuckets);
vector<uint> elinbuckets(nBuckets);
vector<boophf_t*> MPHFs(nBuckets);


void compactBucket(uint start, uint n ){
	for(uint i(start);i<start+n;++i){
		rewind(vFiles[i]);
		u_int64_t *data= (u_int64_t * ) malloc(elinbuckets[i]*sizeof(u_int64_t));
		// data = new u_int64_t(elinbuckets[i]);
		fread(data, sizeof(u_int64_t), elinbuckets[i], vFiles[i]);
		// for(uint ii(0);ii<elinbuckets[i];++ii){
		// 	printf("there is %llu \n", data[ii]);
		// }
		// printf("Go ? !!!\n");
		auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(data), static_cast<const u_int64_t*>(data+elinbuckets[i]));
		// printf("Go ?? !!!\n");
		MPHFs[i]= new boomphf::mphf<u_int64_t,hasher_t>(elinbuckets[i],data_iterator,1,gammaFactor,false);
		// printf("Go com???pactions !!!\n");
		// delete(data);
		free(data);
	}
}


int main (int argc, char* argv[])
{
	//if we want a random seed from timer
//	typedef std::chrono::high_resolution_clock myclock;
//	myclock::time_point beginning = myclock::now();

	bool check_correctness = false;
	bool bench_lookup = false;
	bool save_mphf = false;
	bool load_mphf = false;
	bool buckets = false;
	bool from_disk = false;

	if(argc >=2 )
	{
		nelem = strtoul(argv[1], NULL,0);
	}

	if(argc >=3 ){
		nthreads = atoi(argv[2]);
	}

	for (int ii=3; ii<argc; ii++)
	{
		if(!strcmp("-check",argv[ii])) check_correctness= true;
		if(!strcmp("-bench",argv[ii])) bench_lookup= true;
		if(!strcmp("-save",argv[ii])) save_mphf= true;
		if(!strcmp("-load",argv[ii])) load_mphf= true;
		if(!strcmp("-fromdisk",argv[ii])) from_disk= true;
		if(!strcmp("-buckets",argv[ii])) buckets= true;

	}

	if(buckets){
		clock_t begin, end;
		double t_begin,t_end; struct timeval timet;
		gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
		for(uint i(0);i<nBuckets;++i){
			vFiles[i]=fopen(("bucket"+to_string(i)).c_str(),"w+");
		}
		u_int64_t step = ULLONG_MAX / nelem;
		u_int64_t current = 0;
		for (u_int64_t i = 0; i < nelem; i++){
			current = current + step;
			// cout<<current<<endl;
			u_int64_t hash=korenXor(current)%nBuckets;
			//printf("hash %i  nBuckets %i\n",hash,nBuckets);
			fwrite(&current, sizeof(u_int64_t), 1, vFiles[hash]);
			++elinbuckets[hash];
		}

		vector<uint> nb_elem_in_previous_buckets (nBuckets);
		nb_elem_in_previous_buckets[0] = 0 ;
		for(int ii=1; ii<nBuckets; ii++ )
		{
			nb_elem_in_previous_buckets[ii] = nb_elem_in_previous_buckets[ii-1] + elinbuckets[ii-1];
		}
		
		printf("Go compactions !!!\n");
		
		double integ;
		assert(  modf((double)nBuckets/nthreads ,&integ) == 0  );
		vector<thread> threads;
		for(uint n(0);n<nthreads;++n){
			threads.push_back(thread(compactBucket,n*nBuckets/nthreads,nBuckets/nthreads));
		}
		for(auto &t : threads){t.join();}

		for(uint i(0);i<nBuckets;++i){
			fclose(vFiles[i]);
			remove(("bucket"+to_string(i)).c_str());
		}
		gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);

		double elapsed = t_end - t_begin;

		printf("BooPHF constructed perfect hash for %llu keys in %.2fs\n", nelem,elapsed);
		// cin.get();
		if(check_correctness){
			u_int64_t step2 = ULLONG_MAX / nelem;
			u_int64_t current2 = 0;
			u_int64_t range_problems(0);
			begin = clock();
			for (u_int64_t i = 0; i < nelem; i++){
				current2 += step2;
				// cout<<current2<<endl;
				uint64_t hash=korenXor(current2)%nBuckets;
				u_int64_t mphf_value = MPHFs[hash]->lookup(current2)+  nb_elem_in_previous_buckets [hash];
				// cout<<mphf_value<<endl;
				if(mphf_value>=nelem){
					range_problems++;
					printf("there is %llu problems \n", range_problems);
				}
			}
			printf("there is %llu problems\n", range_problems);

			end = clock();
			printf("BooPHF %llu lookups in  %.2fs,  approx  %.2f ns per lookup \n", nelem, (double)(end - begin) / CLOCKS_PER_SEC,  ((double)(end - begin) / CLOCKS_PER_SEC)*1000000000/nelem);
		}
		return EXIT_SUCCESS;
	}

	FILE * key_file = NULL;

	if(from_disk){
		key_file = fopen("keyfile","w+");
	}

	//creating some random input data
	//	srandomdev(); //init random generator
	//	data = (u_int64_t * ) calloc(nelem,sizeof(u_int64_t));
	//	for (u_int64_t i = 0; i < nelem; i++)
	//		data[i] = random64();

	uint64_t ii, jj;

	// obtain a seed from the timer

//	myclock::duration d = myclock::now() - beginning;
//	unsigned seed2 = d.count();

	if(!from_disk){
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

		//	data = (u_int64_t * ) calloc(nelem,sizeof(u_int64_t));
		//		data[0] = 0;
		//		u_int64_t step = ULLONG_MAX / nelem;
		//			for (u_int64_t i = 1; i < nelem; i++)
		//				data[i] = data[i-1] +step;
		//
	}
	else{
		//methode simple pas besoin de  de-dupliquer, mais pas "random"
		u_int64_t step = ULLONG_MAX / nelem;
		u_int64_t current = 0;
		fwrite(&current, sizeof(u_int64_t), 1, key_file);
		//printf("from disk  nelem %lli  step %llu \n",nelem,step);
		for (u_int64_t i = 1; i < nelem; i++)
		{
			current = current + step;
			fwrite(&current, sizeof(u_int64_t), 1, key_file);
		}
		fclose(key_file);
	}

	std::string output_filename;
	output_filename = "saved_mphf";
	boophf_t * bphf;

	clock_t begin, end;
	double t_begin,t_end; struct timeval timet;

	if(!load_mphf){
		printf("Construct a BooPHF with  %lli elements (%lli MB for holding elems in ram) \n",nelem,nelem*sizeof(u_int64_t)/1024LL/1024LL);

		///create the boophf

		gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);

		//MPHF CREATION
		if(from_disk)
		{

			auto data_iterator = file_binary("keyfile");
			bphf = new boomphf::mphf<u_int64_t,hasher_t>(nelem,data_iterator,nthreads,gammaFactor);
		}
		else
		{
			auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(data), static_cast<const u_int64_t*>(data+nelem));
			bphf = new boomphf::mphf<u_int64_t,hasher_t>(nelem,data_iterator,nthreads,gammaFactor);
		}

		gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);

		double elapsed = t_end - t_begin;

		printf("BooPHF constructed perfect hash for %llu keys in %.2fs\n", nelem,elapsed);

		printf("boophf  bits/elem : %f\n",(float) (bphf->totalBitSize())/nelem);
	}
	else{
		//assumes the mphf was saved before, reload it
		bphf = new boomphf::mphf<u_int64_t,hasher_t>();

		printf("Loading a BooPHF with  %lli elements (%lli MB for holding elems in ram) \n",nelem,nelem*sizeof(u_int64_t)/1024LL/1024LL);

		gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);

		std::ifstream is(output_filename, std::ios::binary);
		bphf->load(is);

		gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
		double elapsed = t_end - t_begin;

		printf("BooPHF re-loaded perfect hash for %llu keys in %.2fs\n", nelem,elapsed);
		printf("boophf  bits/elem : %f\n",(float) (bphf->totalBitSize())/nelem);
	}


	if(save_mphf)
	{
		std::ofstream os(output_filename, std::ios::binary);
		bphf->save(os);
	}

	u_int64_t mphf_value;

	if(check_correctness  && from_disk)
	{
		u_int64_t nb_collision_detected = 0;
		u_int64_t range_problems = 0;
		char * check_table = (char * ) calloc(nelem,sizeof(char));

		auto data_iterator = file_binary("keyfile");

		for (auto const& val: data_iterator)
		{
			mphf_value = bphf->lookup(val);
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

	if(check_correctness  && !from_disk )	//test the mphf
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

	if(bench_lookup && ! from_disk)
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

	if(!from_disk){
		free(data);
	}

	return EXIT_SUCCESS;
}
