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


inline double get_time_usecs() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return double(tv.tv_sec) * 1000000 + double(tv.tv_usec);
}


uint64_t random64 (){
	uint64_t low, high,res;
	low = random();
	high = random();

	res = (high << 32) + low;
	return res;
}


typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
typedef boomphf::mphf<  u_int64_t, hasher_t  > boophf_t;


// iterator from disk file of u_int64_t with buffered read,   todo template
class bfile_iterator : public std::iterator<std::forward_iterator_tag, u_int64_t>{
	public:

	bfile_iterator()
	: _is(nullptr)
	, _pos(0) ,_inbuff (0), _cptread(0)
	{
		_buffsize = 10000;
		_buffer = (u_int64_t *) malloc(_buffsize*sizeof(u_int64_t));
	}

	bfile_iterator(const bfile_iterator& cr)
	{
		_buffsize = cr._buffsize;
		_pos = cr._pos;
		_is = cr._is;
		_buffer = (u_int64_t *) malloc(_buffsize*sizeof(u_int64_t));
		 memcpy(_buffer,cr._buffer,_buffsize*sizeof(u_int64_t) );
		_inbuff = cr._inbuff;
		_cptread = cr._cptread;
		_elem = cr._elem;
	}

	bfile_iterator(FILE* is): _is(is) , _pos(0) ,_inbuff (0), _cptread(0)
	{
		_buffsize = 10000;
		_buffer = (u_int64_t *) malloc(_buffsize*sizeof(u_int64_t));
		int reso = fseek(_is,0,SEEK_SET);
		advance();
	}

	~bfile_iterator()
	{
		if(_buffer!=NULL)
			free(_buffer);
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
		_pos++;

		if(_cptread >= _inbuff)
		{
			int res = fread(_buffer,sizeof(u_int64_t),_buffsize,_is);
			_inbuff = res; _cptread = 0;

			if(res == 0)
			{
				_is = nullptr;
				_pos = 0;
				return;
			}
		}

		_elem = _buffer[_cptread];
		_cptread ++;
	}
	u_int64_t _elem;
	FILE * _is;
	unsigned long _pos;

	u_int64_t * _buffer; // for buffered read
	int _inbuff, _cptread;
	int _buffsize;
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



//simple iterator to generate list of distinct u_int64_t keys (not random, just equally distributed in [0;ULLONG_MAX])
class uint64_iterator : public std::iterator<std::forward_iterator_tag, u_int64_t>{
public:
	
	uint64_iterator() : _nb_elem(0) , _curr(ULLONG_MAX), _step(0),_nb_iterated(0)
	{
	}
	
	uint64_iterator(const uint64_iterator& cr)
	{
		_nb_elem = cr._nb_elem;
		_curr = cr._curr;
		_step = cr._step;
		_stop = cr._stop;
		_nb_iterated = cr._nb_iterated;
	}
	
	uint64_iterator(u_int64_t nb_elem, u_int64_t stop ) : _nb_elem(nb_elem) , _curr(0),_nb_iterated(0), _stop(stop)
	{
		_step =  ULLONG_MAX / _nb_elem;
	}
	
	~uint64_iterator()
	{
	}
	
	u_int64_t const& operator*()  {  return _curr;  }
	
	uint64_iterator& operator++()
	{

		_curr = _curr + _step;
		_nb_iterated++;
		if(_nb_iterated >= _stop) _nb_elem = 0;
		return *this;
	}
	
	friend bool operator==(uint64_iterator const& lhs, uint64_iterator const& rhs)
	{
		if (!lhs._nb_elem || !rhs._nb_elem)  {  if (!lhs._nb_elem && !rhs._nb_elem) {  return true; } else {  return false;  } }
		assert(lhs._nb_elem == rhs._nb_elem);
		return rhs._curr == lhs._curr;
	}
	
	friend bool operator!=(uint64_iterator const& lhs, uint64_iterator const& rhs)  {  return !(lhs == rhs);  }

	u_int64_t _nb_elem;
	u_int64_t _curr;
	u_int64_t _step;
	u_int64_t _stop;
	u_int64_t _nb_iterated;
};

//will generate range over [0;ULLONG_MAX]
//elem equally spaced , separated  by  ULLONG_MAX/nb_elem starting at 0
//will stop after nbiter
class uint64_range{
public:
	uint64_range(u_int64_t nb_elem, u_int64_t nbiter) : _nb_elem(nb_elem), _stop(nbiter)
	{
		
	}
	~uint64_range()
	{
		
	}
	
	uint64_iterator begin() const
	{
		return uint64_iterator(_nb_elem,_stop);
	}
	uint64_iterator end() const
	{
		return uint64_iterator();
	}
	
	
private:
	u_int64_t _nb_elem;
	u_int64_t _stop;

};

//stolen from emphf
struct stats_accumulator {
	stats_accumulator()
	: m_n(0)
	, m_mean(0)
	, m_m2(0)
	{}

	void add(double x)
	{
		m_n += 1;
		auto delta = x - m_mean;
		m_mean += delta / m_n;
		m_m2 += delta * (x - m_mean);
	}

	double mean() const
	{
		return m_mean;
	}

	double variance() const
	{
		return m_m2 / (m_n - 1);
	}

	double relative_stddev() const
	{
		return std::sqrt(variance()) / mean() * 100;
	}

private:
	double m_n;
	double m_mean;
	double m_m2;
};


//PARAMETERS
u_int64_t nelem = 1000*1000;
uint nthreads = 1; //warning must be a divisor of nBuckets
double gammaFactor = 1.0;
bool write_each = false;

u_int64_t nb_in_bench_file;


uint64_t korenXor(uint64_t x){
	x ^= (x << 21);
	x ^= (x >> 35);
	x ^= (x << 4);
	return x;
}


uint nBuckets = 96;
uint nMphfByBucket(96);
vector<FILE*> vFiles(nBuckets);
vector<uint> elinbuckets(nBuckets*nMphfByBucket);
vector<boophf_t> MPHFs(nBuckets*nMphfByBucket);


void multipleMPHF(const vector<vector<u_int64_t>>& datas, uint start, uint n,uint bucketNum){
	for(uint ii(start);ii<start+n;++ii){
		auto data_iterator2 = boomphf::range(static_cast<const u_int64_t*>(&datas[ii][0]), static_cast<const u_int64_t*>(&datas[ii][0]+datas[ii].size()));
		MPHFs[bucketNum*nMphfByBucket+ii]=  boomphf::mphf<u_int64_t,hasher_t>(datas[ii].size(),data_iterator2,1,gammaFactor,write_each,false);
	}
}


void compactBucket(uint start, uint n){
	//foreach bucket

	for(uint i(start);i<start+n;++i){
		auto data_iterator = file_binary(("bucket"+to_string(i)).c_str());
		vector<vector<u_int64_t>> datas(nMphfByBucket);

		for(uint ireserve(0);ireserve<nMphfByBucket;++ireserve){
			datas[ireserve].reserve(elinbuckets[i*nMphfByBucket+ireserve]);
		}

		// we put element in memory
		for(auto it(data_iterator.begin());it!=data_iterator.end();++it){
			datas[(korenXor(*it)%(nMphfByBucket))].push_back(*it);
		}

		vector<thread> threads;
		for(uint tn(0);tn<1;++tn){
			threads.push_back(thread(multipleMPHF,datas,tn*(nMphfByBucket/1),nMphfByBucket/1,i));
		}
		// threads.push_back(thread(multipleMPHF,datas,(nthreads)*(nMphfByBucket/nthreads),nMphfByBucket-(nthreads)*(nMphfByBucket/nthreads),i));

		for(auto &t : threads){t.join();}
	}
}


template <typename phf_t,typename Range>
int check_mphf_correctness (phf_t * bphf, Range const& input_range){
		u_int64_t nb_collision_detected = 0;
		u_int64_t range_problems = 0;
		u_int64_t mphf_value;
			boomphf::bitVector check_table (nelem);
		//auto data_iterator = file_binary("keyfile");

		for (auto const& val: input_range)
		{
			mphf_value = bphf->lookup(val);
			//printf("%llu  mphf_value %llu\n",val,mphf_value);

			if(mphf_value>=nelem)
			{
				range_problems++; continue;
			}
			if(check_table[mphf_value]==0)
			{
				check_table.set(mphf_value);
			}
			else
			{
				//printf("collision for val %lli : \n",mphf_value);
				printf("collision for %llu  mphf_value %llu\n",val,mphf_value);

				nb_collision_detected++;
			}
		}



		if(nb_collision_detected ==  0 && range_problems ==0)
		{
			printf(" --- boophf working correctly --- \n");
			return 0;

		}
		else
		{
			printf("!!! problem, %llu collisions detected; %llu out of range !!!\n",nb_collision_detected,range_problems);
			return 1;

		}
}


template <typename phf_t,typename Range>
void bench_mphf_lookup (phf_t * bphf, Range const& input_range){


	vector<u_int64_t> sample;
	u_int64_t mphf_value;

	//copy sample in ram
	for (auto const& key: input_range) {
		sample.push_back(key);
	}

	printf("bench lookups  sample size %lu \n",sample.size());
	//bench procedure taken from emphf
	stats_accumulator stats;
	double tick = get_time_usecs();
	size_t lookups = 0;
	static const size_t lookups_per_sample = 1 << 16;
	u_int64_t dumb=0;
	double elapsed;
	size_t runs = 10;

	for (size_t run = 0; run < runs; ++run) {
		for (size_t ii = 0; ii < sample.size(); ++ii) {

			mphf_value = bphf->lookup(sample[ii]);
			//do some silly work
			dumb+= mphf_value;

			if (++lookups == lookups_per_sample) {
				elapsed = get_time_usecs() - tick;
				stats.add(elapsed / (double)lookups);
				tick = get_time_usecs();
				lookups = 0;
			}
		}
	}
	printf("BBhash bench lookups average %.2f ns +- stddev  %.2f %%   (fingerprint %llu)  \n", 1000.0*stats.mean(),stats.relative_stddev(),dumb);

}


//#include "bucketing.h"


int main (int argc, char* argv[]){
	//if we want a random seed from timer
//	typedef std::chrono::high_resolution_clock myclock;
//	myclock::time_point beginning = myclock::now();

	bool check_correctness = false;
	bool bench_lookup = false;
	bool save_mphf = false;
	bool load_mphf = false;
	bool buckets = false;
	bool from_disk = true;
	bool bench_lookup_out = false;
	bool on_the_fly= false;
	 write_each = true;
	if(argc <4 ){
		printf("Usage :\n");
		printf("%s <nelem> <nthreads> <gamma>  [options]\n",argv[0]);
		printf("Options:\n");
		printf("\t-check  (check correctness of mphf)\n");
		printf("\t-bench  (bench query time of mphf)\n");
		printf("\t-save\n");
		printf("\t-load\n");
		printf("\t-inram\n");
		printf("\t-nodisk  (do not write each intermediate level on disk)\n");
		printf("\t-buckets\n");
		printf("\t-outquery (bench the fp rate of the mphf)\n");  // bench fp rate
		printf("\t-onthefly (generates key on the fly without storing them on disk or in ram)\n");
	


		return EXIT_FAILURE;
	}

	if(argc >=4 ){
		nelem = strtoul(argv[1], NULL,0);
		nthreads = atoi(argv[2]);
		gammaFactor = atoi(argv[3]);
	}


	
	for (int ii=4; ii<argc; ii++){
		if(!strcmp("-check",argv[ii])) check_correctness= true;
		if(!strcmp("-bench",argv[ii])) bench_lookup= true;
		if(!strcmp("-save",argv[ii])) save_mphf= true;
		if(!strcmp("-load",argv[ii])) load_mphf= true;
		if(!strcmp("-inram",argv[ii])) from_disk= false;
		if(!strcmp("-buckets",argv[ii])) buckets= true;
		if(!strcmp("-outquery",argv[ii])) bench_lookup_out= true;
		if(!strcmp("-onthefly",argv[ii])) on_the_fly= true;
		if(!strcmp("-nodisk",argv[ii])) write_each= false;

	}

	
	if(gammaFactor == 0) {
		fprintf(stderr,"gamma value error\n");
		fprintf(stderr,"Usage should be \n");
		fprintf(stderr,"%s <nelem> <nthreads> <gamma>  [options]\n",argv[0]);
		exit(1);
	}

	
//	///testing terator
//	printf("testing terator :\n");
//	uint64_range terator_in =  uint64_range(10,10);
//	uint64_range terator = uint64_range(terator_in);
//	
//	
//	for (auto const& val: terator) {
//		printf("%llu \n",val);
//
//	}
//	exit(0);
	
	
	FILE * key_file = NULL;
	FILE * bench_file = NULL;

	if(from_disk){
		key_file = fopen("keyfile","w+");
	}

	uint64_t ii, jj;

	/////  generation of keys
	if(!from_disk && !buckets){
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
	}
	else{
		if(!on_the_fly)
		{
			//methode simple pas besoin de  de-dupliquer, mais pas "random"
			u_int64_t step = ULLONG_MAX / nelem;
			// u_int64_t step = 100 / nelem;
			u_int64_t current = 0;
			for (u_int64_t i = 0; i < nelem; i++)
			{
				fwrite(&current, sizeof(u_int64_t), 1, key_file);
				//printf("%llu \n",current);
				current = current + step;
			}
			fclose(key_file);
			printf("key file generated \n");
			
			if(bench_lookup)
			{
				bench_file = fopen("benchfile","w+");
				//create a test file
				//if n < 10 M take all elements, otherwise regular sample to have 10 M elements
				u_int64_t stepb =  nelem  / 10000000;
				if(stepb==0) stepb=1;
				auto data_iterator = file_binary("keyfile");
				u_int64_t cpt = 0;
				nb_in_bench_file=0;
				for (auto const& key: data_iterator) {
					if( (cpt % stepb) == 0)
					{
						fwrite(&key, sizeof(u_int64_t), 1, bench_file);
						nb_in_bench_file++;
					}
					cpt++;
				}
				fclose(bench_file);
			}
		}
	}


	vector<uint> nb_elem_in_previous_buckets (nBuckets*nMphfByBucket);

	
	if(buckets){
		clock_t begin, end;
		double t_begin,t_end; struct timeval timet;
		gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);

		
		for(uint i(0);i<nBuckets;++i){
			vFiles[i]=fopen(("bucket"+to_string(i)).c_str(),"w+");
		}

		printf("splitting keys ..\n");
		
		double tick_split = get_time_usecs();

		
		int buffsize = 10000;
		vector < vector<u_int64_t> > buffers (nBuckets);
		for(int ii=0; ii<nBuckets;ii++)
			buffers[ii].reserve(buffsize);
		
		auto data_iterator = file_binary("keyfile");
		for (auto const& key: data_iterator) {
			u_int64_t hash=(korenXor(key)%(nBuckets*nMphfByBucket)/nMphfByBucket);
			
			if(buffers[hash].size()==buffsize)
			{
				fwrite(buffers[hash].data(), sizeof(u_int64_t), buffers[hash].size(), vFiles[hash]);
				buffers[hash].clear();//hope it keeps capacity intact
			}
			buffers[hash].push_back(key);
			
			++elinbuckets[korenXor(key)%(nBuckets*nMphfByBucket)];
		}
		//flush buffers
		for(int ii=0; ii<nBuckets;ii++)
		{
			fwrite(buffers[ii].data(), sizeof(u_int64_t), buffers[ii].size(), vFiles[ii]);
		}

		for (int ii=0; ii<vFiles.size(); ii++) {
			fclose(vFiles[ii]);
		}

		nb_elem_in_previous_buckets[0] = 0 ;
		for(int ii=1; ii<nBuckets*nMphfByBucket; ii++ ){
			nb_elem_in_previous_buckets[ii] = nb_elem_in_previous_buckets[ii-1] + elinbuckets[ii-1];
		}

		double elapsed_split = get_time_usecs() - tick_split;
		printf("time key split  %.2f s \n", elapsed_split/1000000.0);

		printf("Go compactions !!!\n");

		double integ;

		assert(  modf((double)nBuckets/nthreads ,&integ) == 0  );
		vector<thread> threads;
		for(uint n(0);n<nthreads;++n){
			threads.push_back(thread(compactBucket,n*nBuckets/nthreads,nBuckets/nthreads));
		}
		for(auto &t : threads){t.join();}
		//~ compactBucket(0, nBuckets);

		for(uint i(0);i<nBuckets;++i){
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
			u_int64_t nb_collision_detected(0);
			begin = clock();
			boomphf::bitVector check_table (nelem);
			for (u_int64_t i = 0; i < nelem; i++){

				uint64_t hash=korenXor(current2)%(nBuckets*nMphfByBucket);
				u_int64_t mphf_value = MPHFs[hash].lookup(current2)+  nb_elem_in_previous_buckets [hash];
				if(mphf_value>=nelem){
					range_problems++;
					printf("there is %llu problems \n", range_problems);
				}
				if(check_table[mphf_value]==0)
				{
					check_table.set(mphf_value);
				}
				else
				{
					//printf("collision for val %lli \n",mphf_value);
					nb_collision_detected++;
				}
				current2 += step2;
			}
			printf("there is %llu problems\n", range_problems);
			printf("there is %llu coll\n", nb_collision_detected);

			end = clock();
			//printf("BooPHF %llu lookups in  %.2fs,  approx  %.2f ns per lookup \n", nelem, (double)(end - begin) / CLOCKS_PER_SEC,  ((double)(end - begin) / CLOCKS_PER_SEC)*1000000000/nelem);
		}

		if(bench_lookup)
		{

			auto input_range = file_binary("benchfile");

			vector<u_int64_t> sample;
			u_int64_t mphf_value;
			
			//copy sample in ram
			for (auto const& key: input_range) {
				sample.push_back(key);
			}
			
			printf("bench lookups  sample size %lu \n",sample.size());
			//bench procedure taken from emphf
			stats_accumulator stats;
			double tick = get_time_usecs();
			size_t lookups = 0;
			static const size_t lookups_per_sample = 1 << 16;
			u_int64_t dumb=0;
			double elapsed;
			size_t runs = 10;
			
			for (size_t run = 0; run < runs; ++run) {
				for (size_t ii = 0; ii < sample.size(); ++ii) {
					
					uint64_t hash=korenXor(sample[ii])%(nBuckets*nMphfByBucket);
					mphf_value = MPHFs[hash].lookup(sample[ii]) +  nb_elem_in_previous_buckets [hash];
					dumb+= mphf_value;
					
					//do some silly work
					
					if (++lookups == lookups_per_sample) {
						elapsed = get_time_usecs() - tick;
						stats.add(elapsed / (double)lookups);
						tick = get_time_usecs();
						lookups = 0;
					}
				}
			}
			printf("BBhash buckets bench lookups average %.2f ns +- stddev  %.2f %%   (fingerprint %llu)  \n", 1000.0*stats.mean(),stats.relative_stddev(),dumb);
			
			///
		}



		return EXIT_SUCCESS;
	}



	std::string output_filename;
	output_filename = "saved_mphf";
	boophf_t * bphf = NULL;;


	clock_t begin, end;
	double t_begin,t_end; struct timeval timet;

	if(!load_mphf){
		printf("Construct a BooPHF with  %lli elements  \n",nelem);
		///create the boophf

		gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);

		//MPHF CREATION
		if (on_the_fly)
		{
			auto data_iterator =  uint64_range(nelem,nelem) ;
			bphf = new boomphf::mphf<u_int64_t,hasher_t>(nelem,data_iterator,nthreads,gammaFactor,write_each);
		}
		else if(from_disk)
		{
			auto data_iterator = file_binary("keyfile");
			bphf = new boomphf::mphf<u_int64_t,hasher_t>(nelem,data_iterator,nthreads,gammaFactor,write_each);
		}
		else
		{
			auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(data), static_cast<const u_int64_t*>(data+nelem));
			bphf = new boomphf::mphf<u_int64_t,hasher_t>(nelem,data_iterator,nthreads,gammaFactor,write_each);
		}

		gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);

		double elapsed = t_end - t_begin;


		printf("BooPHF constructed perfect hash for %llu keys in %.2fs\n", nelem,elapsed);
		printf("boophf  bits/elem : %f\n",(float) (bphf->totalBitSize())/nelem);

	}
	else{
		//assumes the mphf was saved before, reload it
		bphf = new boomphf::mphf<u_int64_t,hasher_t>();

		printf("Loading a BooPHF with  %lli elements  \n",nelem);

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


	if(check_correctness && on_the_fly )
	{
		auto data_iterator =  uint64_range(nelem,nelem) ;
		check_mphf_correctness(bphf ,data_iterator);
	}
	else if(check_correctness && from_disk )
	{
		auto data_iterator = file_binary("keyfile");
		check_mphf_correctness(bphf ,data_iterator);
	}
	else if(check_correctness &&  !from_disk )
	{
		auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(data), static_cast<const u_int64_t*>(data+nelem));
		check_mphf_correctness(bphf ,data_iterator);
	}

	
	if(bench_lookup &&  on_the_fly)
	{
		auto data_iterator = uint64_range(nelem,1000000) ;
		
		bench_mphf_lookup(bphf,data_iterator);
		
	}
	else if(bench_lookup &&  from_disk)
	{
		auto data_iterator = file_binary("benchfile");

		bench_mphf_lookup(bphf,data_iterator);

	}
	else if(bench_lookup &&  !from_disk)
	{
		auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(data), static_cast<const u_int64_t*>(data+nelem));

		bench_mphf_lookup(bphf,data_iterator);
	}
	
	if(bench_lookup_out)
	{
		int nrandom = 100000000; //10000000
		static std::mt19937_64 rng;
		rng.seed(std::mt19937_64::default_seed); //default seed
		u_int64_t *  data_random = (u_int64_t * ) calloc(nrandom,sizeof(u_int64_t));
		for (u_int64_t i = 0; i < nrandom; i++){
			data_random[i] = rng();
		//	printf("%llu \n",data_random[i]);
		}
		u_int64_t mphf_value;
		u_int64_t dumb=0;
		u_int64_t nb_fp =0;
		u_int64_t nb_out_of_range =0;

		
		for (size_t ii = 0; ii < nrandom; ++ii) {
			
			mphf_value = bphf->lookup(data_random[ii]);

			//printf("m %llu \n",mphf_value);

			if(mphf_value != ULLONG_MAX)
			{
				nb_fp++;
			}
			
			if((mphf_value != ULLONG_MAX) &&  (mphf_value >= nelem))
			{
				nb_out_of_range++;
			}
			
		}
		
		
		double tick = get_time_usecs();
		
		for (size_t ii = 0; ii < nrandom; ++ii) {
			mphf_value = bphf->lookup(data_random[ii]);
			//do some silly work
			dumb+= mphf_value;

		}
		double elapsed = get_time_usecs() - tick;

		printf("query %i elem  out of set  FP rate %.2f   nb issues %llu    lookup %.2f  ns \n",nrandom, nb_fp/(float)nrandom,nb_out_of_range
			   ,1000.0*elapsed/(double)nrandom );
	}
	

	if(!from_disk){
		free(data);
	}

		delete bphf;

	return EXIT_SUCCESS;
}
