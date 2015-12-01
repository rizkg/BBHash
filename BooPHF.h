// BooPHF library
// intended to be a minimal perfect hash function with fast and low memory construction, at the cost of higher bits/elem than other state of the art libraries once built.
// should work with arbitray large number of elements, based on a cascade of bloom filters

#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include <array>
#include <unordered_map>
#include <vector>
#include <assert.h>
#include <sys/time.h>
#include <string.h>



//#include <algorithm>
//#include <cmath>

namespace boomphf {
	
////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark utils
////////////////////////////////////////////////////////////////

	inline unsigned int popcount_32(unsigned int x)
	{
		unsigned int m1 = 0x55555555;
		unsigned int m2 = 0x33333333;
		unsigned int m4 = 0x0f0f0f0f;
		unsigned int h01 = 0x01010101;
		x -= (x >> 1) & m1;               /* put count of each 2 bits into those 2 bits */
		x = (x & m2) + ((x >> 2) & m2);   /* put count of each 4 bits in */
		x = (x + (x >> 4)) & m4;          /* put count of each 8 bits in partie droite  4bit piece*/
		return (x * h01) >> 24;           /* returns left 8 bits of x + (x<<8) + ... */
	}
	
	
	inline unsigned int popcount_64(uint64_t x)
	{
		
		unsigned int low = x & 0xffffffff ;
		unsigned int high = ( x >> 32) & 0xffffffff ;
		
		return (popcount_32(low) + popcount_32(high));
	}
	
	
	///// progress bar
	class Progress
	{
	public:
		int timer_mode;
		struct timeval timestamp;
		double heure_debut, heure_actuelle ;
		std::string   message;

		uint64_t done;
		uint64_t todo;
		int subdiv ; // progress printed every 1/subdiv of total to do
		double partial;
		double partial_threaded[64];
		uint64_t done_threaded[64];
		
		double steps ; //steps = todo/subidv
		
		void init(uint64_t ntasks, const char * msg)
		{
			message = std::string(msg);
			gettimeofday(&timestamp, NULL);
			heure_debut = timestamp.tv_sec +(timestamp.tv_usec/1000000.0);
			
			//fprintf(stderr,"| %-*s |\n",98,msg);
			
			todo= ntasks;
			done = 0;
			partial =0;
			for (int ii=0; ii<64;ii++) partial_threaded[ii]=0;
			for (int ii=0; ii<64;ii++) done_threaded[ii]=0;
			subdiv= 100;
			steps = (double)todo / (double)subdiv;
			
			if(!timer_mode)
			{
				fprintf(stderr,"[");fflush(stderr);
			}
		}
		
		void finish()
		{
			set(todo);
			if(timer_mode)
				fprintf(stderr,"\n");
			else
				fprintf(stderr,"]\n");
			
			fflush(stderr);
			todo= 0;
			done = 0;
			partial =0;
			
		}
		void finish_threaded()// called by only one of the threads
		{
			done = 0;
			double rem = 0;
			for (int ii=0; ii<64;ii++) done += (done_threaded[ii] );
			for (int ii=0; ii<64;ii++) partial += (partial_threaded[ii] );
			
			finish();
			
		}
		void inc(uint64_t ntasks_done)
		{
			done += ntasks_done;
			partial += ntasks_done;
			
			
			while(partial >= steps)
			{
				if(timer_mode)
				{
					gettimeofday(&timestamp, NULL);
					heure_actuelle = timestamp.tv_sec +(timestamp.tv_usec/1000000.0);
					double elapsed = heure_actuelle - heure_debut;
					double speed = done / elapsed;
					double rem = (todo-done) / speed;
					if(done>todo) rem=0;
					int min_e  = (int)(elapsed / 60) ;
					elapsed -= min_e*60;
					int min_r  = (int)(rem / 60) ;
					rem -= min_r*60;

				fprintf(stderr,"%c[%s]  %-5.3g%%   elapsed: %3i min %-2.0f sec   remaining: %3i min %-2.0f sec",13,
						message.c_str(),
						100*(double)done/todo,
						min_e,elapsed,min_r,rem);

				}
				else
				{
					fprintf(stderr,"-");fflush(stderr);
				}
				partial -= steps;
			}
			
			
		}
		
		void inc(uint64_t ntasks_done, int tid) //threads collaborate to this same progress bar
		{
			partial_threaded[tid] += ntasks_done;
			done_threaded[tid] += ntasks_done;
			while(partial_threaded[tid] >= steps)
			{
				if(timer_mode)
				{
					struct timeval timet;
					double now;
					gettimeofday(&timet, NULL);
					now = timet.tv_sec +(timet.tv_usec/1000000.0);
					uint64_t total_done  = 0;
					for (int ii=0; ii<64;ii++) total_done += (done_threaded[ii] );
					double elapsed = now - heure_debut;
					double speed = total_done / elapsed;
					double rem = (todo-total_done) / speed;
					if(total_done > todo) rem =0;
					int min_e  =  (int)(elapsed / 60) ;
					elapsed -= min_e*60;
					int min_r  =  (int)(rem / 60) ;
					rem -= min_r*60;
					
					fprintf(stderr,"%c[%s]  %-5.3g%%   elapsed: %3i min %-2.0f sec   remaining: %3i min %-2.0f sec",13,
							message.c_str(),
							100*(double)done/todo,
							min_e,elapsed,min_r,rem);
				}
				else
				{
					fprintf(stderr,"-");fflush(stderr);
				}
				partial_threaded[tid] -= steps;
				
			}
			
		}
		
		void set(uint64_t ntasks_done)
		{
			if(ntasks_done > done)
				inc(ntasks_done-done);
		}
		Progress () :     timer_mode(0) {}
		//include timer, to print ETA ?
	};
	
	
	
	
////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark hasher
////////////////////////////////////////////////////////////////
	
	typedef std::array<uint64_t,7> hash_set_t;
	
	
	
	template <typename Item> class HashFunctors
	{
		
		
	public:
		
		/** Constructor.
		 * \param[in] nbFct : number of hash functions to be used
		 * \param[in] seed : some initialization code for defining the hash functions. */
		HashFunctors ()
		{
			_nbFct = 7; // use 7 hash func
			_user_seed = 0;
			generate_hash_seed ();
		}
		
		//return one hash
		uint64_t operator ()  (const Item& key, size_t idx)  {  return hash64 (key, _seed_tab[idx]);  }
		
		//this one returns all the 7 hashes
		//maybe use xorshift instead, for faster hash compute
		hash_set_t operator ()  (const Item& key)
		{
			hash_set_t	 hset;
			
			for(size_t ii=0;ii<7; ii++)
			{
				hset[ii] =  hash64 (key, _seed_tab[ii]);
			}
			return hset;
		}
		
	private:
		
		
		inline static uint64_t hash64 (Item key, uint64_t seed)
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
		
		/* */
		void generate_hash_seed ()
		{
			static const uint64_t rbase[MAXNBFUNC] =
			{
				0xAAAAAAAA55555555ULL,  0x33333333CCCCCCCCULL,  0x6666666699999999ULL,  0xB5B5B5B54B4B4B4BULL,
				0xAA55AA5555335533ULL,  0x33CC33CCCC66CC66ULL,  0x6699669999B599B5ULL,  0xB54BB54B4BAA4BAAULL,
				0xAA33AA3355CC55CCULL,  0x33663366CC99CC99ULL
			};
			
			for (size_t i=0; i<MAXNBFUNC; ++i)  {  _seed_tab[i] = rbase[i];  }
			for (size_t i=0; i<MAXNBFUNC; ++i)  {  _seed_tab[i] = _seed_tab[i] * _seed_tab[(i+3) % MAXNBFUNC] + _user_seed ;  }
		}
		
		size_t _nbFct;
		
		static const size_t MAXNBFUNC = 10;
		uint64_t _seed_tab[MAXNBFUNC];
		uint64_t _user_seed;
	};
	
	
////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark iterators
////////////////////////////////////////////////////////////////
	
	template <typename Iterator>
	struct iter_range
	{
		iter_range(Iterator b, Iterator e)
		: m_begin(b)
		, m_end(e)
		{}
		
		Iterator begin() const
		{ return m_begin; }
		
		Iterator end() const
		{ return m_end; }
		
		Iterator m_begin, m_end;
	};
	
	template <typename Iterator>
	iter_range<Iterator> range(Iterator begin, Iterator end)
	{
		return iter_range<Iterator>(begin, end);
	}
	
////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark BitVector
////////////////////////////////////////////////////////////////
	
	class bitVector {
		
	public:
		
		bitVector() : _size(0)
		{}
		
		bitVector(uint64_t n) : _size(n)
		{
			_nchar  = (1+n/64LL);
			_bitArray = (uint64_t *) calloc (_nchar,sizeof(uint64_t));
		}
		
		~bitVector()
		{
			free(_bitArray);
		}
		
		void resize(uint64_t newsize)
		{
			//printf("bitvector resize from  %llu bits to %llu \n",_size,newsize);
			_nchar  = (1+newsize/64LL);
			_bitArray = (uint64_t *) realloc(_bitArray,_nchar*sizeof(uint64_t));
			_size = newsize;
		}
		
		size_t size() const
		{
			return _size;
		}
		
		uint64_t bitSize() const {return (_nchar*64LL + _ranks.capacity()*64 );}
		
		//clear whole array
		void clear()
		{
			memset(_bitArray,0,_nchar*sizeof(uint64_t));
		}
		
		//clear interval, only works with start and size multiple of 64
		void clear(uint64_t start, size_t size)
		{
			assert( (start & 63) ==0);
			assert( (size & 63) ==0);
			memset(_bitArray + (start/64),0,(size/64)*sizeof(uint64_t));
		}
		
		//for debug purposes
		void print() const
		{
			printf("bit array of size %lli: \n",_size);
			for(uint64_t ii = 0; ii< _size; ii++)
			{
				if(ii%10==0)
					printf(" (%llu) ",ii);
				int val = (_bitArray[ii >> 6] >> (ii & 63 ) ) & 1;
				printf("%i",val);
			}
			printf("\n");
			
			printf("rank array : size %lu \n",_ranks.size());
			for (uint64_t ii = 0; ii< _ranks.size(); ii++)
			{
				printf("%llu :  %lli,  ",ii,_ranks[ii]);
			}
			printf("\n");
		}
		
		//return value at pos
		uint64_t operator[](uint64_t pos) const
		{
			return (_bitArray[pos >> 6] >> (pos & 63 ) ) & 1;
			//return (_bitArray[pos >> 6] & (1ULL << (pos & 63) )  ) ;//renvoit nul ou non nul (pas 1 )
		}
		
		
		//set bit pos to 1
		void set(uint64_t pos)
		{
			if(pos>=_size)
			{
				resize( _size * 1.1);
			}
			assert(pos<_size);
			//_bitArray [pos >> 6] |=   (1ULL << (pos & 63) ) ;
			__sync_fetch_and_or (_bitArray + (pos >> 6), (1ULL << (pos & 63)) );

		}
		
		//set bit pos to 0
		void reset(uint64_t pos)
		{
			//_bitArray [pos >> 6] &=   ~(1ULL << (pos & 63) ) ;
			__sync_fetch_and_and (_bitArray + (pos >> 6), ~(1ULL << (pos & 63) ));

		}
		
		void build_ranks()
		{
			
			_ranks.reserve(2+ _size/_nb_bits_per_rank_sample);
			
			uint64_t curent_rank = 0;
			for (size_t ii = 0; ii < _nchar; ii++) {
				if (((ii*64)  % _nb_bits_per_rank_sample) == 0) {
					_ranks.push_back(curent_rank);
				}
				curent_rank +=  popcount_64(_bitArray[ii]);
			}
		}
		
		uint64_t rank(uint64_t pos) const
		{
			uint64_t word_idx = pos / 64;
			uint64_t word_offset = pos % 64;
			uint64_t block = pos / _nb_bits_per_rank_sample;
			uint64_t r = _ranks[block];
			for (uint64_t w = block * _nb_bits_per_rank_sample / 64; w < word_idx; ++w) {
				r += popcount_64( _bitArray[w] );
			}
			uint64_t mask = (uint64_t(1) << word_offset ) - 1;
			r += popcount_64( _bitArray[word_idx] & mask);
			
			return r;
		}
		
	protected:
		uint64_t* _bitArray;
		uint64_t _size;
		uint64_t _nchar;
		
		 // epsilon =  64 / _nb_bits_per_rank_sample   bits
		// additional size for rank is epsilon * _size
		static const uint64_t _nb_bits_per_rank_sample = 256; //256
		std::vector<uint64_t> _ranks;
	};
	
////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark bloom
////////////////////////////////////////////////////////////////
	
	//simple blocked bloom class, does not compute hashes internally, must be given the N hashes as a parameter
	//thus does not need to be templated
	//(allows reuse of hashes for other things externally)
	
	
	u_int8_t bit_mask [] = {
		0x01,  //00000001
		0x02,  //00000010
		0x04,  //00000100
		0x08,  //00001000
		0x10,  //00010000
		0x20,  //00100000
		0x40,  //01000000
		0x80   //10000000
	};
	
	
	//blocked bloom
	class bbloom
	{
	public:
		
		bbloom (uint64_t tai_bloom, size_t nbHash = 7, size_t block_nbits = 12)
		:  _n_hash_func(nbHash), _blooma(0), _tai(tai_bloom+2*(1<<block_nbits)), _nchar(0), _nbits_BlockSize(block_nbits)
		{
			_nchar  = (1+_tai/8LL);
			_blooma = (unsigned char *) malloc (_nchar*sizeof(unsigned char));
			memset (_blooma, 0, _nchar*sizeof(unsigned char));
			
			_mask_block = (1<<_nbits_BlockSize) - 1;
			_reduced_tai = _tai -  2*(1<<_nbits_BlockSize) ;//2* for neighbor coherent
		}
		virtual ~bbloom ()  { free (_blooma); }
		
		size_t getNbHash () const { return _n_hash_func; }
		uint64_t bitSize() const {return _nchar*8;}

		bool contains (hash_set_t & hashes)
		{
			uint64_t h0 =  hashes[0] % _reduced_tai;
			
			if ((_blooma[h0 >> 3 ] & bit_mask[h0 & 7]) == 0 )  {  return false;  }
			
			for (size_t i=1; i<_n_hash_func; i++)
			{
				uint64_t h1 = h0  + (hashes[i] & _mask_block ) ;
				if ((_blooma[h1 >> 3 ] & bit_mask[h1 & 7]) == 0)  {  return false;  }
			}
			return true;
		}
		
		void insert (hash_set_t & hashes)
		{
			uint64_t h0 = hashes[0] % _reduced_tai;
			__sync_fetch_and_or (_blooma + (h0 >> 3), bit_mask[h0 & 7]);
			
			for (size_t i=1; i< _n_hash_func; i++)
			{
				uint64_t h1 = h0  + (hashes[i] & _mask_block )   ;
				__sync_fetch_and_or (_blooma + (h1 >> 3), bit_mask[h1 & 7]);
			}
		}
		
	private:
		
		size_t _n_hash_func;
		u_int8_t* _blooma;
		uint64_t _tai;
		uint64_t _nchar;
		uint64_t _mask_block;
		size_t    _nbits_BlockSize;
		uint64_t _reduced_tai;
	};
	

////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark level
////////////////////////////////////////////////////////////////
	
	class level{
	public:
		level(){ }
		
		~level() { delete bloom;}
		
		uint64_t idx_begin;
		uint64_t hash_domain;
		bbloom * bloom;

	};
	

////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark mphf
////////////////////////////////////////////////////////////////



	template <typename elem_t, typename Hasher_t>
	class mphf {
		
		
		
	public:
		mphf()
		{}
		
		~mphf()
		{
			for(int ii=0; ii<_nb_levels; ii++)
			{
				delete _levels[ii];
			}
			delete _bitset;
		}
		

		template <typename Range>
		mphf( size_t n, Range const& input_range, double gamma = 2.4) :
		_gamma(gamma), _hash_domain(size_t(ceil(double(n) * gamma))), _nelem(n)
		{
			
			setup();
			
			_progressBar.timer_mode=1;
			_progressBar.init( _nelem * _nb_levels ,"Building BooPHF");
			
			initLevel0(input_range);
			
			for(int ii = 0; ii< _nb_levels-1; ii++)
			{
				processLevel(input_range,ii);
			}
			
			_progressBar.finish();

			
			_bitset->build_ranks();
			
		}
		
		
		uint64_t lookup(elem_t elem)
		{
			auto hashes = _hasher(elem);
			uint64_t non_minimal_hp;
			
			int level =  getLevel(hashes);
			
			
			if( level == (_nb_levels-1))
			{
				non_minimal_hp = _final_hash[elem];
			}
			else
			{
				non_minimal_hp = _levels[level]->idx_begin + ( hashes[level] %  _levels[level]->hash_domain);
			}
			
			uint64_t minimal_hp = _bitset->rank(non_minimal_hp);
			
			return minimal_hp;
		}
		
		uint64_t totalBitSize()
		{
			uint64_t bloomsizes = 0;
			for (int ii=0; ii< _nb_levels-1; ii++)
			{
				bloomsizes+= _levels[ii]->bloom->bitSize();
			}
			uint64_t totalsize = _bitset->bitSize() + bloomsizes + _final_hash.size()*42*8 ;  // unordered map takes approx 42B per elem [personal test]
			
			
			printf("Bitarray    %12llu  bits (%.2f %%)   (array + ranks )\n",
				   _bitset->bitSize(), 100*(float)_bitset->bitSize()/totalsize);
			printf("Blooms      %12llu  bits (%.2f %%)\n",
				   bloomsizes, 100*(float)bloomsizes/totalsize);
			printf("final hash  %12lu  bits (%.2f %%) (nb in final hash %lu)\n",
				   _final_hash.size()*42*8, 100*(float)(_final_hash.size()*42*8)/totalsize,
				   _final_hash.size() );
			return totalsize;
		}
		

		
		private :
		
		void setup()
		{
			double proba_collision = 1.0 - ( (1.0 -  pow( 1.0 - 1.0/(_gamma*_nelem), _nelem ) )  *_gamma) ;
			double sum_geom = proba_collision / (1.0 - proba_collision);
			//printf("proba collision %f  sum_geom  %f  bitvector size %lli \n",proba_collision,sum_geom,(uint64_t) (_hash_domain + (2.0*proba_collision*_hash_domain)));
			
			int bloom_bit_per_elem = 12;//12 / 5
			
			_nb_levels = 6;
			
			_levels = (level **) malloc(_nb_levels * sizeof(level *) );
			
			//build levels
			uint64_t previous_idx =0;
			for(int ii=0; ii<_nb_levels; ii++)
			{
				_levels[ii] = new level();
				
				_levels[ii]->idx_begin = previous_idx;
				
				// round size to nearet superior multiple of 64, makes it easier to clear a level
				_levels[ii]->hash_domain =  (( (uint64_t) (_hash_domain * pow(proba_collision,ii)) + 63) / 64 ) * 64;
				previous_idx += _levels[ii]->hash_domain;
				
				//printf("build level %i bit array : start %llu, size %llu, bloom %llu \n",ii,_levels[ii]->idx_begin,_levels[ii]->hash_domain,(uint64_t)( pow(proba_collision,ii+1) * _nelem * bloom_bit_per_elem) );
				
				if(ii<(_nb_levels-1))
					_levels[ii]->bloom = new bbloom( pow(proba_collision,ii+1) * _nelem * bloom_bit_per_elem);
				else
					_levels[ii]->bloom  = NULL; // last level has no bloom
				
			}
			
			_bitset = new bitVector(previous_idx);  // approx :   sum_geom * _hash_domain
			
		}
		
		//query the level, up to given maxlevel (useful during construction)
		int getLevel(hash_set_t & hashes, int maxlevel = 100)
		{
			int level = 0;
			for (int ii=0; ii<(_nb_levels-1) &&  ii < maxlevel ; ii++ )
			{
				if( ! _levels[ii]->bloom->contains(hashes))
				{
					break;
				}
				level++;
			}
			return level;
		}
		
		//insert into bitarray ot into  bloom if collision in bitarray
		void insertIntoLevel(hash_set_t &  hashes, int i)
		{
			uint64_t hashi = _levels[i]->idx_begin +    hashes[i] % _levels[i]->hash_domain;
			
			if( ! (*_bitset)[hashi])
			{
				_bitset->set(hashi);
			}
			else
			{
				_levels[i]->bloom->insert(hashes);
			}
			
		}
		
		
		//first loop to fill level 0
		template <typename Range>
		void initLevel0(Range const& input_range)
		{
			uint64_t nb_done =0;
			for (auto const &val: input_range) {
				
				auto hashes = _hasher(val);
				
				insertIntoLevel(hashes,0);
				
				nb_done++;
				if((nb_done&1023) ==0 ) {_progressBar.inc(nb_done);nb_done=0; }
			}
			
		}
		
		
		//loop combining debloom of level i and insert into level i+1
		template <typename Range>
		void processLevel(Range const& input_range,int i)
		{
			
			uint64_t hashi;
			uint64_t nb_done =0;
			
			//clear level before deblooming
			_bitset->clear(_levels[i]->idx_begin, _levels[i]->hash_domain);
			
			uint64_t hashidx = _levels[i+1]->idx_begin;
			
			////////// loop to detect FP of bloom and insert into next level
			for (auto const& val: input_range) {
				
				auto hashes = _hasher(val);
				int level = getLevel(hashes, i+1);
				
				if(level == i+1)
				{
					//insert to level i+1 : either next level of the cascade or final hash if last level reached
					
					if(i+1== _nb_levels-1) //stop cascade here, insert into exact hash
					{
						_bitset->set(hashidx);
						_final_hash[val] = hashidx;
						hashidx++;
					}
					else
					{
						insertIntoLevel(hashes,i+1);
					}
				}
				else if (level == i)
				{
					//insert to table level i
					hashi = _levels[i]->idx_begin + hashes[i] % _levels[i]->hash_domain;
					//printf("nc %llu --> %llu\n",val,hashi);
					_bitset->set(hashi);
				}
				
				nb_done++;
				if((nb_done&1023) ==0 ) {_progressBar.inc(nb_done);nb_done=0; }
			}
		}
	
	private:
		level ** _levels;
		int _nb_levels;
		Hasher_t _hasher;
		bitVector * _bitset;
		double _gamma;
		uint64_t _hash_domain;
		uint64_t _nelem;
		std::unordered_map<elem_t,uint64_t> _final_hash;
		Progress _progressBar;

	};

	
}