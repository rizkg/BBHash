//
//  bit_vector_array.h
//  bank_request
//
//  Created by Pierre Peterlongo on 01/03/16.
//
//

#ifndef bank_request_bit_vector_array_h
#define bank_request_bit_vector_array_h
using namespace std;
#include <boost/dynamic_bitset.hpp>
typedef boost::dynamic_bitset<> _bitset;


class bitArraySet {
public:
    
    bitArraySet() : _nb_elements(0), _nb_bit_per_element(0)
    {
    }
    
    bitArraySet(const uint64_t nb_elements, const int nb_bit_per_element) : _nb_elements(nb_elements), _nb_bit_per_element(nb_bit_per_element)
    {
        _table = _bitset (_nb_elements*_nb_bit_per_element);
//        set_mask();
    }
    
    void set_i(const uint64_t i, const _bitset toset){
        const uint64_t start = i*_nb_bit_per_element;
        for (int p=0; p<_nb_bit_per_element; p++) {
            _table[start+p] = toset[p];
        }
    }
    
     _bitset get_i(const uint64_t i){
//        return (_table >> ((_nb_elements-i-1)*_nb_bit_per_element)) ; // I'm afraid this code creates a copy of the original bitset, that can be huge.
         const uint64_t start = i*_nb_bit_per_element;
         _bitset result (_nb_bit_per_element);
         for (int p=0; p<_nb_bit_per_element; p++) {
             result[p]=_table[start+p];
         }
         return result;
    }
    
    _bitset get_set(){
        return _table;
    }
    
private:
    
//    void set_mask(){
//        _mask = _bitset (_nb_bit_per_element,( 2 << _nb_bit_per_element )-1);
//    }
//    
  
    _bitset _table;
//    _bitset _mask;
    uint64_t _nb_elements;
    int _nb_bit_per_element;
};



#endif
