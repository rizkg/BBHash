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
#include <memory> // for make_shared


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
		unsigned int high = ( x >> 32LL) & 0xffffffff ;
		
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
			subdiv= 1000;
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
							100*(double)total_done/todo,
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
        uint64_t operator ()  (const Item& key, size_t idx)  const {  return hash64 (key, _seed_tab[idx]);  }
       
        uint64_t hashWithSeed(const Item& key, uint64_t seed)  const {  return hash64 (key, seed);  }
		
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
	
/* alternative hash functor based on xorshift, taking a single hash functor as input.
we need this 2-functors scheme because HashFunctors won't work with unordered_map.
(rayan) 
*/

    // wrapper around HashFunctors to return only one value instead of 7
    template <typename Item> class SingleHashFunctor 
	{
	public:
		uint64_t operator ()  (const Item& key, uint64_t seed=0xAAAAAAAA55555555ULL) const  {  return hashFunctors.hashWithSeed(key, seed);  }
		
	private:
		HashFunctors<Item> hashFunctors;
		
	};

	
	
	// the SingleHasher_t must have  operator()(elem_t key, uint64_t seed)
	//this class simply generates a list of seeds
	template <typename Item, class SingleHasher_t> class IndepHashFunctors
	{
		
	public:

		IndepHashFunctors ()
		{
			generate_hash_seed ();
		}
		hash_set_t operator ()  (const Item& key)
		{
			hash_set_t	 hset;
			
			for(size_t ii=0;ii<7; ii++)
			{
				hset[ii] =  singleHasher (key, _seed_tab[ii]);
			}
			return hset;
		}
		
	private:
		void generate_hash_seed ()
		{
			static const uint64_t rbase[MAXNBFUNC] =
			{
				0xAAAAAAAA55555555ULL,  0x33333333CCCCCCCCULL,  0x6666666699999999ULL,  0xB5B5B5B54B4B4B4BULL,
				0xAA55AA5555335533ULL,  0x33CC33CCCC66CC66ULL,  0x6699669999B599B5ULL,  0xB54BB54B4BAA4BAAULL,
				0xAA33AA3355CC55CCULL,  0x33663366CC99CC99ULL
			};
			
			for (size_t i=0; i<MAXNBFUNC; ++i)  {  _seed_tab[i] = rbase[i];  }
			for (size_t i=0; i<MAXNBFUNC; ++i)  {  _seed_tab[i] = _seed_tab[i] * _seed_tab[(i+3) % MAXNBFUNC]  ;  }
		}
		
		static const size_t MAXNBFUNC = 10;
		uint64_t _seed_tab[MAXNBFUNC];
		SingleHasher_t singleHasher;
	};
	
    template <typename Item, class SingleHasher_t> class XorshiftHashFunctors
    {
        /*  Xorshift128*
            Written in 2014 by Sebastiano Vigna (vigna@acm.org)

            To the extent possible under law, the author has dedicated all copyright
            and related and neighboring rights to this software to the public domain
            worldwide. This software is distributed without any warranty.

            See <http://creativecommons.org/publicdomain/zero/1.0/>. */
        /* This is the fastest generator passing BigCrush without
           systematic failures, but due to the relatively short period it is
           acceptable only for applications with a mild amount of parallelism;
           otherwise, use a xorshift1024* generator.

           The state must be seeded so that it is not everywhere zero. If you have
           a nonzero 64-bit seed, we suggest to pass it twice through
           MurmurHash3's avalanching function. */

      //  uint64_t s[ 2 ];

        uint64_t next(uint64_t * s) {
            uint64_t s1 = s[ 0 ];
            const uint64_t s0 = s[ 1 ];
            s[ 0 ] = s0;
            s1 ^= s1 << 23; // a
            return ( s[ 1 ] = ( s1 ^ s0 ^ ( s1 >> 17 ) ^ ( s0 >> 26 ) ) ) + s0; // b, c
        }

        public:
        //this one returns all the 7 hashes
        hash_set_t operator ()  (const Item& key)
        {
			uint64_t s[ 2 ];

            hash_set_t   hset;
            
            hset[0] =  singleHasher (key, 0xAAAAAAAA55555555ULL); 
            hset[1] =  singleHasher (key, 0x33333333CCCCCCCCULL); 
            
            s[0] = hset[0];
            s[1] = hset[1];

            for(size_t ii=2;ii< 7 /* it's much better have a constant here, for inlining; this loop is super performance critical*/; ii++)
            {
                hset[ii] = next(s);
            }

            return hset;
        }
    private:
        SingleHasher_t singleHasher;
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
			_nchar  = (1ULL+n/64ULL);
			_bitArray = (uint64_t *) calloc (_nchar,sizeof(uint64_t));
		}
		
		~bitVector()
		{
			free(_bitArray);
		}
		
		void resize(uint64_t newsize)
		{
			//printf("bitvector resize from  %llu bits to %llu \n",_size,newsize);
			_nchar  = (1ULL+newsize/64ULL);
			_bitArray = (uint64_t *) realloc(_bitArray,_nchar*sizeof(uint64_t));
			_size = newsize;
		}
		
		size_t size() const
		{
			return _size;
		}
		
		uint64_t bitSize() const {return (_nchar*64ULL + _ranks.capacity()*64ULL );}
		
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
			memset(_bitArray + (start/64ULL),0,(size/64ULL)*sizeof(uint64_t));
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
			return (_bitArray[pos >> 6ULL] >> (pos & 63 ) ) & 1;
		}
		
		//atomically   return old val and set to 1
		uint64_t atomic_test_and_set(uint64_t pos)
		{
			uint64_t oldval = 	__sync_fetch_and_or (_bitArray + (pos >> 6), (uint64_t) (1ULL << (pos & 63)) );
			
			return  ( oldval >> (pos & 63 ) ) & 1;
		}
		
		
		//set bit pos to 1
		void set(uint64_t pos)
		{
			assert(pos<_size);
			//_bitArray [pos >> 6] |=   (1ULL << (pos & 63) ) ;
			__sync_fetch_and_or (_bitArray + (pos >> 6ULL), (1ULL << (pos & 63)) );

		}
		
		//set bit pos to 0
		void reset(uint64_t pos)
		{
			//_bitArray [pos >> 6] &=   ~(1ULL << (pos & 63) ) ;
			__sync_fetch_and_and (_bitArray + (pos >> 6ULL), ~(1ULL << (pos & 63) ));

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
			uint64_t word_idx = pos / 64ULL;
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
		static const uint64_t _nb_bits_per_rank_sample = 512; //512 seems ok
		std::vector<uint64_t> _ranks;
	};
	
////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark bloom
////////////////////////////////////////////////////////////////
	
	//simple blocked bloom class, does not compute hashes internally, must be given the N hashes as a parameter
	//thus does not need to be templated
	//(allows reuse of hashes for other things externally)
	
	
    static u_int8_t bit_mask [] = {
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
			_nchar  = (1ULL+_tai/8ULL);
			_blooma = (unsigned char *) malloc (_nchar*sizeof(unsigned char));
			memset (_blooma, 0, _nchar*sizeof(unsigned char));
			
			_mask_block = (1ULL<<_nbits_BlockSize) - 1ULL;
			_reduced_tai = _tai -  2ULL*(1ULL<<_nbits_BlockSize) ;//2* for neighbor coherent
		}
		virtual ~bbloom ()  { free (_blooma); }
		
		size_t getNbHash () const { return _n_hash_func; }
		uint64_t bitSize() const {return _nchar*8;}

		bool contains (hash_set_t & hashes)
		{
			uint64_t h0 =  hashes[0] % _reduced_tai;
			
			if ((_blooma[h0 >> 3ULL ] & bit_mask[h0 & 7]) == 0 )  {  return false;  }
			
			for (size_t i=1; i<_n_hash_func; i++)
			{
				uint64_t h1 = h0  + (hashes[i] & _mask_block ) ;
				if ((_blooma[h1 >> 3ULL ] & bit_mask[h1 & 7]) == 0)  {  return false;  }
			}
			return true;
		}
		
		void insert (hash_set_t & hashes)
		{
			uint64_t h0 = hashes[0] % _reduced_tai;
			__sync_fetch_and_or (_blooma + (h0 >> 3ULL), bit_mask[h0 & 7]);
			
			for (size_t i=1; i< _n_hash_func; i++)
			{
				uint64_t h1 = h0  + (hashes[i] & _mask_block )   ;
				__sync_fetch_and_or (_blooma + (h1 >> 3ULL), bit_mask[h1 & 7]);
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

	
#define NBBUFF 10000
	
	template<typename Range,typename Iterator>
	struct thread_args
	{
		void * boophf;
		Range const * range;
		std::shared_ptr<void> it_p; /* used to be "Iterator it" but because of fastmode, iterator is polymorphic; TODO: think about whether it should be a unique_ptr actually */
		std::shared_ptr<void> until_p; /* to cache the "until" variable */
		int level;
	};
	
	//forward declaration
    template <typename elem_t, typename Hasher_t, typename Range, typename it_type>
	void * thread_initLevel0(void * args);
	
    template <typename elem_t, typename Hasher_t, typename Range, typename it_type>
	void * thread_processLevel(void * args);
	

    /* Hasher_t returns a single hash when operator()(elem_t key) is called.
       if used with XorshiftHashFunctors, it must have the following operator: operator()(elem_t key, uint64_t seed) */
    template <typename elem_t, typename Hasher_t>
	class mphf {
		
        /* this mechanisms gets 7 hashes (for the Bloom filters) out of Hasher_t */
        typedef XorshiftHashFunctors<elem_t,Hasher_t> MultiHasher_t ;
       // typedef HashFunctors<elem_t> MultiHasher_t; // original code (but only works for int64 keys)  (seems to be as fast as the current xorshift)
		//typedef IndepHashFunctors<elem_t,Hasher_t> MultiHasher_t; //faster than xorshift
		
	public:
		mphf()
		{}
		
		~mphf()
		{
			pthread_mutex_destroy(&_mutex);

			for(int ii=0; ii<_nb_levels; ii++)
			{
				delete _levels[ii];
			}
			delete _bitset;
		}
		

		template <typename Range>
		mphf( size_t n, Range const& input_range,int num_thread = 1, bool fastmode = true,  double gamma = 2.5) :
		_gamma(gamma), _hash_domain(size_t(ceil(double(n) * gamma))), _nelem(n), _num_thread(num_thread), _fastmode (fastmode)
		{
			
			setup();
			
			_progressBar.timer_mode=1;
			if(_fastmode)
				_progressBar.init( _nelem * 3 +  ( _nelem * _proba_collision * _proba_collision) * (_nb_levels-3)    ,"Building BooPHF");
			else
				_progressBar.init( _nelem * _nb_levels ,"Building BooPHF");
			
			initLevel0(input_range);
			
			for(int ii = 0; ii< _nb_levels-1; ii++)
			{
				processLevel(input_range,ii);
			}
			
			_progressBar.finish_threaded();

			
			_bitset->build_ranks();
			
			_lastbitsetrank = _bitset->rank( _bitset->size() -1);
			
			//printf("used temp ram for construction : %lli MB \n",setLevel2.capacity()* sizeof(elem_t) /1024ULL/1024ULL);
			
			std::vector<elem_t>().swap(setLevel2);   // clear setLevel2 reallocating
			
		}
		
		
		uint64_t lookup(elem_t elem)
		{
			auto hashes = _hasher(elem);
			uint64_t non_minimal_hp,minimal_hp;
			
			int level =  getLevel(hashes);
			
			
			if( level == (_nb_levels-1))
			{
				minimal_hp = _final_hash[elem] + _lastbitsetrank;
				return minimal_hp;
			}
			else
			{
				non_minimal_hp = _levels[level]->idx_begin + ( hashes[level] %  _levels[level]->hash_domain);
			}
			
			 minimal_hp = _bitset->rank(non_minimal_hp);
			
			return minimal_hp;
		}
		
		uint64_t nbKeys()
		{
            return _nelem;
        }

		uint64_t totalBitSize()
		{
			uint64_t bloomsizes = 0;
			for (int ii=0; ii< _nb_levels-1; ii++)
			{
				bloomsizes+= _levels[ii]->bloom->bitSize();
			}
			uint64_t totalsize = _bitset->bitSize() + bloomsizes + _final_hash.size()*42*8 ;  // unordered map takes approx 42B per elem [personal test] (42B with uint64_t key, would be larger for other type of elem)
			
			
			printf("Bitarray    %12llu  bits (%.2f %%)   (array + ranks )\n",
				   _bitset->bitSize(), 100*(float)_bitset->bitSize()/totalsize);
			printf("Blooms      %12llu  bits (%.2f %%)\n",
				   bloomsizes, 100*(float)bloomsizes/totalsize);
			printf("final hash  %12lu  bits (%.2f %%) (nb in final hash %lu)\n",
				   _final_hash.size()*42*8, 100*(float)(_final_hash.size()*42*8)/totalsize,
				   _final_hash.size() );
			return totalsize;
		}
		
		template <typename Range,typename Iterator>
        void pthread_init0(elem_t * buffer, Range input_range, std::shared_ptr<Iterator> shared_it, std::shared_ptr<Iterator> until_p)
		{
			uint64_t nb_done =0;
			
			int tid =  __sync_fetch_and_add (&_nb_living, 1);
			
			auto until = *until_p;
			
			
			uint64_t inbuff =0;
			
			for (bool isRunning=true;  isRunning ; )
			{
				
				pthread_mutex_lock(&_mutex);
				
				//copy n items into buffer
                for(; inbuff<NBBUFF && (*shared_it)!=until;  ++(*shared_it) /* subtle: if it was sahared++, we would need to implenent operator++(int) for our iterator */)
				{
                    buffer[inbuff]= *(*shared_it); inbuff++;
				}
				
                if((*shared_it)==until) isRunning =false;
				
				pthread_mutex_unlock(&_mutex);
				
				
				//do work on the n elems of the buffer
                for(uint64_t ii=0; ii<inbuff ; ii++)
				{
					elem_t val = buffer[ii];
					
					auto hashes = _hasher(val);
					
					insertIntoLevel(hashes,0);
					
					nb_done++;
					if((nb_done&1023) ==0 ) {_progressBar.inc(nb_done,tid);nb_done=0; }
				}
				
				inbuff = 0;
			}
			
		}
	
		template <typename Range,typename Iterator>
        void pthread_processLevel(elem_t * buffer, Range input_range, std::shared_ptr<Iterator> shared_it, std::shared_ptr<Iterator> until_p, int i)
		{
			uint64_t nb_done =0;
			int tid =  __sync_fetch_and_add (&_nb_living, 1);
			
			auto until = *until_p;
			uint64_t inbuff =0;

			for (bool isRunning=true;  isRunning ; )
			{
				
				//safely copy n items into buffer
				pthread_mutex_lock(&_mutex);
                for(; inbuff<NBBUFF && (*shared_it)!=until;  ++(*shared_it))
				{
                    buffer[inbuff]= *(*shared_it); inbuff++;
				}
                if((*shared_it)==until) isRunning =false;
				pthread_mutex_unlock(&_mutex);
				
				
				//do work on the n elems of the buffer
                for(uint64_t ii=0; ii<inbuff ; ii++)
				{
					elem_t val = buffer[ii];
					
					auto hashes = _hasher(val);
					int level = getLevel(hashes, i+1); //should be safe
					
					if(level == i+1)
					{
						if(i == 1 && _fastmode)
						{
							uint64_t idxl2 = __sync_fetch_and_add(& _idxLevelSetLevel2,1);
							//si depasse taille attendue pour setLevel2, fall back sur slow mode mais devrait pas arriver si hash ok et proba avec nous
							if(idxl2> setLevel2.size())
								_fastmode = false;
							else
								setLevel2[idxl2] = val; // create set E2
						}
						
						//insert to level i+1 : either next level of the cascade or final hash if last level reached
						if(i+1== _nb_levels-1) //stop cascade here, insert into exact hash
						{
							uint64_t hashidx =  __sync_fetch_and_add (& _hashidx, 1);
							
							pthread_mutex_lock(&_mutex); //see later if possible to avoid this, mais pas bcp item vont la
							//_bitset->set(hashidx);
							// calc rank de fin  precedent level qq part, puis init hashidx avec ce rank, direct minimal, pas besoin inser ds bitset et rank
							_final_hash[val] = hashidx;
							pthread_mutex_unlock(&_mutex);
						}
						else
						{
							insertIntoLevel(hashes,i+1); //should be safe
						}
					}
					else if (level == i)
					{
						//insert to table level i
						uint64_t hashi = _levels[i]->idx_begin + hashes[i] % _levels[i]->hash_domain;
						_bitset->set(hashi);
					}
					
					nb_done++;
					if((nb_done&1023) ==0 ) {_progressBar.inc(nb_done,tid);nb_done=0; }

				}
				
				inbuff = 0;
			}
			
		}
		
		private :
		
		void setup()
		{
			pthread_mutex_init(&_mutex, NULL);

			if(_fastmode)
				setLevel2.resize(0.035 * (double)_nelem );
			
			 _proba_collision = 1.0 - ( (1.0 -  pow( 1.0 - 1.0/(_gamma*_nelem), _nelem ) )  *_gamma) ;
			double sum_geom = _proba_collision / (1.0 - _proba_collision);
			//printf("proba collision %f  sum_geom  %f  bitvector size %lli \n",proba_collision,sum_geom,(uint64_t) (_hash_domain + (2.0*proba_collision*_hash_domain)));
			
			int bloom_bit_per_elem = 12;//12 / 5
			
			if(_fastmode)
				_nb_levels = 7; // in fast mode we can afford extra level
			else
				_nb_levels = 6;

			
			_levels = (level **) malloc(_nb_levels * sizeof(level *) );
			
			//build levels
			uint64_t previous_idx =0;
			for(int ii=0; ii<_nb_levels; ii++)
			{
				_levels[ii] = new level();
				
				_levels[ii]->idx_begin = previous_idx;
				
				// round size to nearet superior multiple of 64, makes it easier to clear a level
				_levels[ii]->hash_domain =  (( (uint64_t) (_hash_domain * pow(_proba_collision,ii)) + 63) / 64 ) * 64;
				previous_idx += _levels[ii]->hash_domain;
				
				//printf("build level %i bit array : start %12llu, size %12llu, bloom %12llu \n",ii,_levels[ii]->idx_begin,_levels[ii]->hash_domain,(uint64_t)( pow(_proba_collision,ii+1) * _nelem * bloom_bit_per_elem) );
				
				if(ii<(_nb_levels-1))
					_levels[ii]->bloom = new bbloom( pow(_proba_collision,ii+1) * _nelem * bloom_bit_per_elem);
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
			
			if( _bitset->atomic_test_and_set(hashi) )
			{
				_levels[i]->bloom->insert(hashes);

			}

		}
		
		
		//first loop to fill level 0
		template <typename Range>
		void initLevel0(Range const& input_range)
		{

			_nb_living =0;
			//create  threads
			pthread_t *tab_threads= new pthread_t [_num_thread];
			
			typedef decltype(input_range.begin()) it_type;
			
			thread_args<Range, it_type> t_arg; // meme arg pour tous
			t_arg.boophf = this;
			t_arg.range = &input_range;
			t_arg.it_p =  std::static_pointer_cast<void>(std::make_shared<it_type>(input_range.begin()));
			t_arg.until_p =  std::static_pointer_cast<void>(std::make_shared<it_type>(input_range.end()));
			
			for(int ii=0;ii<_num_thread;ii++)
			{
                pthread_create (&tab_threads[ii], NULL,  thread_initLevel0<elem_t, Hasher_t, Range, it_type>, &t_arg); //&t_arg[ii]
			}
			
			//joining
			for(int ii=0;ii<_num_thread;ii++)
			{
				pthread_join(tab_threads[ii], NULL);
			}
		}
		
		
		//loop combining debloom of level i and insert into level i+1
		template <typename Range>
		void processLevel(Range const& input_range,int i)
		{
			//clear level before deblooming
			_bitset->clear(_levels[i]->idx_begin, _levels[i]->hash_domain);

			_hashidx = 0;
			_idxLevelSetLevel2 =0;
			_nb_living =0;
			//create  threads
			pthread_t *tab_threads= new pthread_t [_num_thread];
			typedef decltype(input_range.begin()) it_type;
			thread_args<Range, it_type> t_arg; // meme arg pour tous
			t_arg.boophf = this;
			t_arg.range = &input_range;
			t_arg.it_p =  std::static_pointer_cast<void>(std::make_shared<it_type>(input_range.begin()));
			t_arg.until_p =  std::static_pointer_cast<void>(std::make_shared<it_type>(input_range.end()));
			t_arg.level = i;
			
			if(i >= 2 && _fastmode)
			{
				auto data_iterator = boomphf::range(static_cast<const elem_t*>( &setLevel2[0]), static_cast<const elem_t*>( (&setLevel2[0]) +setLevel2.size()));
                typedef decltype(data_iterator.begin()) fastmode_it_type;
				t_arg.it_p =  std::static_pointer_cast<void>(std::make_shared<fastmode_it_type>(data_iterator.begin()));
				t_arg.until_p =  std::static_pointer_cast<void>(std::make_shared<fastmode_it_type>(data_iterator.end()));
                
                /* we'd like to do t_arg.it = data_iterator.begin() but types are different;
                    so, casting to (void*) because of that; and we remember the type in the template */
			 
                for(int ii=0;ii<_num_thread;ii++)
                    pthread_create (&tab_threads[ii], NULL,  thread_processLevel<elem_t, Hasher_t, Range, fastmode_it_type>, &t_arg); //&t_arg[ii]
			}
			else
			{
			    for(int ii=0;ii<_num_thread;ii++)
                    pthread_create (&tab_threads[ii], NULL,  thread_processLevel<elem_t, Hasher_t, Range, decltype(input_range.begin())>, &t_arg); //&t_arg[ii]
			}
			//joining
			for(int ii=0;ii<_num_thread;ii++)
			{
				pthread_join(tab_threads[ii], NULL);
			}
		}
	
	private:
		level ** _levels;
		int _nb_levels;
        MultiHasher_t _hasher;
		bitVector * _bitset;
		double _gamma;
		uint64_t _hash_domain;
		uint64_t _nelem;
        std::unordered_map<elem_t,uint64_t,Hasher_t> _final_hash;
		Progress _progressBar;
		int _nb_living;
		int _num_thread;
		uint64_t _hashidx;
		double _proba_collision;
		
		uint64_t _lastbitsetrank;
		uint64_t _idxLevelSetLevel2;
		
		bool _fastmode; // fast build mode , requires that 3.5 % elem are loaded in ram
		std::vector< elem_t > setLevel2; // this set should represent approx  3.5 % of total nb of elements, greatly speed up levels 3 and later

	public:
		pthread_mutex_t _mutex;
	};

////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark threading
////////////////////////////////////////////////////////////////
	
	//bon sang de template
    template <typename elem_t, typename Hasher_t, typename Range, typename it_type>
	void * thread_initLevel0(void * args)
	{
		if(args ==NULL) return NULL;
		
		thread_args<Range,it_type> *targ = (thread_args<Range,it_type>*) args;
		mphf<elem_t, Hasher_t>  * obw = (mphf<elem_t, Hasher_t > *) targ->boophf;
		elem_t * buffer =  (elem_t *)  malloc(NBBUFF*sizeof(elem_t));
		pthread_mutex_t * mutex =  & obw->_mutex;
		
		//get starting iterator for this thread, must be protected (must not be currently used by other thread to copy elems in buff)
		pthread_mutex_lock(mutex);
        std::shared_ptr<it_type> startit = std::static_pointer_cast<it_type>(targ->it_p);
        std::shared_ptr<it_type> until = std::static_pointer_cast<it_type>(targ->until_p);
		pthread_mutex_unlock(mutex);

		obw->pthread_init0(buffer, *(targ->range), startit, until);
		
		free(buffer);
		return NULL;
	}
	
    template <typename elem_t, typename Hasher_t, typename Range, typename it_type>
	void * thread_processLevel(void * args)
	{
		if(args ==NULL) return NULL;
		
		thread_args<Range,it_type> *targ = (thread_args<Range,it_type>*) args;
		
		mphf<elem_t, Hasher_t>  * obw = (mphf<elem_t, Hasher_t > *) targ->boophf;
		int level = targ->level;
		elem_t * buffer =  (elem_t *)  malloc(NBBUFF*sizeof(elem_t));
		pthread_mutex_t * mutex =  & obw->_mutex;
		
		pthread_mutex_lock(mutex); // from comment above: "//get starting iterator for this thread, must be protected (must not be currently used by other thread to copy elems in buff)"
        std::shared_ptr<it_type> startit = std::static_pointer_cast<it_type>(targ->it_p);
        std::shared_ptr<it_type> until_p = std::static_pointer_cast<it_type>(targ->until_p);
		pthread_mutex_unlock(mutex);
		
		obw->pthread_processLevel(buffer, *(targ->range), startit, until_p, level); // i
		
		free(buffer);
		return NULL;
	}
}
