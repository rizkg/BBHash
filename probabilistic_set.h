//
//  probabilistic_set.h
//  bank_request
//
//  Created by Pierre Peterlongo on 01/03/16.
//
//

#ifndef probabilistic_set_h
#define probabilistic_set_h
using namespace std;

#include <iostream>
#include "native_bit_vector_array.h"

using namespace std;



class probabilisticSet {
public:
    probabilisticSet(const uint64_t nb_elements, const int fingerprint_size): _nb_elements(nb_elements), _fingerprint_size(fingerprint_size){
        _bas = bitArraySet(nb_elements,fingerprint_size);
        _fingerprint_range = (uint64_t)1<<_fingerprint_size;
    }
    
    void add(const uint64_t i, const uint64_t key){
        uint64_t fingerprint = korenXor(key)%_fingerprint_range;
        _bas.set_i(i,fingerprint);
    }
    
    bool exists(const uint64_t i, const uint64_t key){
        uint64_t fingerprint = korenXor(key)%_fingerprint_range;
        uint64_t stored_fingerprint = _bas.get_i(i);
        if (fingerprint == stored_fingerprint) return true;
        return false;
    }
    
private:
    
    inline uint64_t korenXor(uint64_t x){
        x ^= (x << 21);
        x ^= (x >> 35);
        x ^= (x << 4);
        return x;
    }
    
    bitArraySet _bas;
    int _fingerprint_size;
    uint64_t _fingerprint_range;
    uint64_t _nb_elements;
};



#endif
