//
//  test.cpp
//  bank_request
//
//  Created by Pierre Peterlongo on 01/03/16.
//
//
#include <iostream>
#include "bit_vector_array.h"
#include "probabilistic_set.h"

using namespace std;

typedef unsigned long long u64b;
u64b _hash(char *key, u64b len)
{
	u64b result = 0;
	while (len)
	{
		result += *key;
		key++;
		len--;
	}
	
	return result;
}
int main (int argc, char* argv[]){
    
    if(argc<3){
		printf("Usage :\n");
		printf("%s <nelem> <nbitperelem>  [options]\n",argv[0]);
		return EXIT_FAILURE;
        
    }
    cout<<"####################################"<<endl;
    cout<<"TESTING BISET"<<endl;
    cout<<"####################################"<<endl;
    
    int nelem = atoi(argv[1]);
    int nbitperelem =atoi(argv[2]);
    bitArraySet bas(nelem,nbitperelem);
    
    cout<<"bas created"<<endl;
    
    _bitset toset(nbitperelem);
    for(int i=0;i<nbitperelem;i++)
        toset[i]=i%2;

    cout<<"insert "<<toset<<" postion 2"<<endl;
    
//    cout<<bas.get_set()<<endl;
    
    bas.set_i(2,toset);
    
    //    cout<<bas.get_set()<<endl;
    
    cout<<"validation of position 2:"<<bas.get_i(2)<<endl;

    
    
    
    cout<<"####################################"<<endl;
    cout<<"TESTING PROBABILISTIC SET"<<endl;
    cout<<"####################################"<<endl;
    
    probabilisticSet ps (nelem,nbitperelem);
    
//    ps.add(test_element,test_element);
    
//    cout<<ps.exists(test_element,test_element)<<endl;
    
    // Fill the full vector array (with key = index for now)
    for(uint64_t i=0;i<nelem;i++) {
        ps.add(i,i);
    }
    
    
    // Assert each element is present (no FN)
    for(uint64_t i=0;i<nelem;i++) {
        if(!ps.exists(i,i)){
            cout<<"ERROR: FN "<<i<<endl;
            
            return EXIT_FAILURE;
        }
    }

    // count the number of FP among one million random queries
    std::srand(std::time(0)); // use current time as seed for random generator
    uint64_t nb_FP=0;
    uint64_t nb_queries = 10000000;
    for(uint64_t i=0;i<nb_queries;i++) {
        if(ps.exists(i%nelem,std::rand())) nb_FP++;
    }

    cout<<nb_FP<<" FP among "<<nb_queries<<" queries "<<(100*nb_FP)/(float)(nb_queries)<<" %"<<endl;
    return EXIT_SUCCESS;
}

