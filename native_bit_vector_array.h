//
//  bit_vector_array.h
//  bank_request
//
//  Created by Pierre Peterlongo on 01/03/16.
//
//

#ifndef bank_request_native_bit_vector_array_h
#define bank_request_native_bit_vector_array_h
using namespace std;

class bitArraySet {
public:
    bitArraySet() : _nb_elements(0), _nb_bit_per_element(0), _nb_bit_per_unit(64)
    {
    }
    
    bitArraySet(const uint64_t nb_elements, const int nb_bit_per_element) : _nb_elements(nb_elements), _nb_bit_per_element(nb_bit_per_element),  _nb_bit_per_unit(64)
    {
        if (_nb_bit_per_element>_nb_bit_per_unit){
            cerr<<"Cannot create a bitArraySet for elements bigger than "<<_nb_bit_per_unit<<" bits"<<endl;
            exit(0);
        }
        
//        cout<<_nb_bit_per_element/(float)_nb_bit_per_unit<<endl;
        _nb_unit = _get_starting_unit_indice(_nb_elements)+1;
        
        set_mask();
        _table = (uint64_t *)malloc(sizeof(uint64_t)*_nb_unit);
        if(_table == NULL){
            cerr<<"Cannot allocate bitArraySet, exit"<<endl;
            exit(0);
        }
//        cout << "_nb_unit = "<<_nb_unit<<endl;
    }
    
    
    
    void set_i(const uint64_t indice_element, uint64_t to_set){
        if (to_set>_mask_unit) {
            cerr<<"Cannot add elements bigger than "<<_mask_unit<<" ("<<_nb_bit_per_unit<<" bits), element "<<to_set<<" to big"<<endl;
            exit(1);
        }
        const uint64_t starting_unit_indice = _get_starting_unit_indice(indice_element);
        const int starting_position_in_the_unit = _get_starting_position_in_the_unit(indice_element, starting_unit_indice);
        const int excluded_ending_position_in_the_unit=starting_position_in_the_unit+_nb_bit_per_element;
        /** case 1: the whole element is contained in the unit */
//        cout<<"excluded_ending_position_in_the_unit "<<excluded_ending_position_in_the_unit<<endl;
        if (excluded_ending_position_in_the_unit<=_nb_bit_per_unit) {
            // 012345670123456701234567 (positions in each unit)
            // 000000001111111122222222 (3 units, each of 8 bits)
            //          ----            (element starts unit 1, position 1, excluded_ending_position_in_the_unit = 1+4 = 5)
            // 11111111 = 11101010      (bits contained in element 1)
            // 1011                     (bit to put in element 1, position 1) --> element 1 becomes 1.1011.010
            // idea:
            // 1/ create an eraser_mask 10000111 (for erasing the old value)
            // 2/ res = old_value & mask
            // 3/ shift the new value by 3 bits (3=8-5=_nb_bit_per_unit-excluded_ending_position_in_the_unit)
            // 4/ res = res | shifted
            const uint64_t eraser_mask = ~(get_mask(_nb_bit_per_element)<<(_nb_bit_per_unit-excluded_ending_position_in_the_unit));
//            cout<<get_mask(_nb_bit_per_element)<<" ~eraser_mask "<<(get_mask(_nb_bit_per_element)<<(_nb_bit_per_unit-excluded_ending_position_in_the_unit))<<endl;
//            cout<<" erased table = "<<(_table[starting_unit_indice] & eraser_mask)<<endl;
            _table[starting_unit_indice] = (_table[starting_unit_indice] & eraser_mask) | (to_set<<(_nb_bit_per_unit-excluded_ending_position_in_the_unit));
//            cout<<"table1 "<<_table[starting_unit_indice]<<endl;
        }
        
        /** case 2: the element spans 2 units */
        else{
//            cout<<"table1avt "<<_table[starting_unit_indice]<<endl;
            // 012345670123456701234567 (positions in each unit)
            // 000000001111111122222222 (3 units, each of 8 bits)
            //                ----       (element starts unit 1, position 7, excluded_ending_position_in_the_unit = 7+4 = 11>nb_bit_per_unit)
            const int excluded_ending_position_in_the_next_unit = excluded_ending_position_in_the_unit - _nb_bit_per_unit;
            // excluded_ending_position_in_the_next_unit = 11-8=3
            
            //---------- starting unit ------------
            // idea:
            // 1/ create an eraser_mask 11111110 that erases last bits corresponding to new value (1 bit to zero = _nb_bit_per_unit-starting_position_in_the_unit)
            // 2/ res = old_value & eraser_mask
            // 3/ res = res | (new value >> 3) (3 = excluded_ending_position_in_the_next_unit)
            // --> the starting using is changed with res
            uint64_t eraser_mask = _mask_unit<<(_nb_bit_per_unit-starting_position_in_the_unit);
            _table[starting_unit_indice] = (_table[starting_unit_indice] & eraser_mask) | (to_set>>excluded_ending_position_in_the_next_unit);
            
            //---------- next affected unit ------
            // idea:
            // 1/ create an eraser_mask 00011111 (3=excluded_ending_position_in_the_next_unit) get_mask(8-3)
            // 2/ res = old_value & eraser_mask
            // 3/ shift the new value by 5 bits to the left
            // 4/ res = res | shifted
            eraser_mask = get_mask(_nb_bit_per_unit-excluded_ending_position_in_the_next_unit);
            _table[starting_unit_indice+1] = (_table[starting_unit_indice+1] & eraser_mask) | (to_set<<(_nb_bit_per_unit-excluded_ending_position_in_the_next_unit));
//            cout<<"table1 "<<_table[starting_unit_indice]<<endl;
//            cout<<"table2 "<<_table[starting_unit_indice+1]<<endl;
        }
    }
    
    uint64_t get_i(const uint64_t indice_element){
        const uint64_t starting_unit_indice = _get_starting_unit_indice(indice_element);
        const int starting_position_in_the_unit = _get_starting_position_in_the_unit(indice_element, starting_unit_indice);
//        cout<<"element "<<indice_element<<" in array composed of "<<_nb_unit<<" units of size "<<_nb_bit_per_unit<<endl;
//        cout<<"starting unit "<<starting_unit_indice<<" position "<<starting_position_in_the_unit<<endl;
        
        const int excluded_ending_position_in_the_unit=starting_position_in_the_unit+_nb_bit_per_element;
//        cout <<"excluded_ending_position_in_the_unit="<<excluded_ending_position_in_the_unit<<endl;
        /** case 1: the whole element is contained in the unit */
        
        // 012345670123456701234567 (positions in each unit)
        // 000000001111111122222222 (3 units, each of 8 bits)
        //          ----            (element starts unit 1, position 1, excluded_ending_position_in_the_unit = 1+4 = 5)
        // 11111111 = 11101010      (bits contained in element 1)
        //             ----
        // (01101010 >> 8-5) & 00001111 = 00001101
//        cout<<"res table1 "<<_table[starting_unit_indice]<<endl;
        if (excluded_ending_position_in_the_unit<=_nb_bit_per_unit) {
            return (_table[starting_unit_indice]>>(_nb_bit_per_unit-excluded_ending_position_in_the_unit)) & _mask_element;
        }
        
//        cout<<"res table2 "<<_table[starting_unit_indice+1]<<endl;
        /** case 2: the element spans 2 units */
        // 012345670123456701234567 (positions in each unit)
        // 000000001111111122222222 (3 units, each of 8 bits)
        //                ----      (element starts unit 1, position 7, excluded_ending_position_in_the_unit = 7+4 = 11>nb_bit_per_unit)
        const int excluded_ending_position_in_the_next_unit = excluded_ending_position_in_the_unit - _nb_bit_per_unit;
        // excluded_ending_position_in_the_next_unit = 11-8=3
        uint64_t result=get_mask(_nb_bit_per_unit-starting_position_in_the_unit);           // set a mask for the last bits present in the starting unit
        result &=_table[starting_unit_indice];                                              // get the last bits of the starting unit
        result = result << excluded_ending_position_in_the_next_unit;                       // give space to the incomming bits from the next unit
        uint64_t starting_bits_of_next_unit =                                               // get the first bits of the next unit that correspond to the end of the element.
            _table[starting_unit_indice+1] >> (_nb_bit_per_unit-excluded_ending_position_in_the_next_unit);
//        cout<<"result "<<result<<endl;
//        cout<<"starting_bits_of_next_unit "<<starting_bits_of_next_unit<<endl;
        result = result | starting_bits_of_next_unit;
        
        return result;
        
    }
    
    
private:
    
    /**
     * Given an element_indice, returns the starting unit_indice
     * BUG WITH MORE THAN ON BILLION ELEMENTS
     */
    inline uint64_t _get_starting_unit_indice(const uint64_t indice_element){
        return ((indice_element*_nb_bit_per_element)/_nb_bit_per_unit);
    }
    
    
    inline int _get_starting_position_in_the_unit(const uint64_t element_indice, const uint64_t starting_unit_indice){
        return element_indice*_nb_bit_per_element - starting_unit_indice*_nb_bit_per_unit;
    }
    
    uint64_t get_mask(const int nb_bit){
        uint64_t res=~0;
        return res>>(_nb_bit_per_unit-nb_bit);
    }
    
    void set_mask(){
        _mask_unit = get_mask(_nb_bit_per_unit);
        _mask_element = get_mask(_nb_bit_per_element);
    }
    
    int  _nb_bit_per_unit; // 64 or better 128
    uint64_t   _nb_unit;      // an unit = an entry of the array. A unit stores 1 or more element(s)
    uint64_t * _table;
    uint64_t   _mask_unit;
    uint64_t   _mask_element;
    
    uint64_t   _nb_elements;     // an element is a stored stuff whose size is <= _nb_bit_per_unit
    int        _nb_bit_per_element;
    
};



#endif
