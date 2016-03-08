//
//  test.cpp
//  bank_request
//
//  Created by Pierre Peterlongo on 01/03/16.
//
//
#include <iostream>
#include "native_bit_vector_array.h"
#include "probabilistic_set.h"
#include <chrono>

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
//    cout<<"####################################"<<endl;
//    cout<<"TESTING BITSET"<<endl;
//    cout<<"####################################"<<endl;
    
    int nelem = atoi(argv[1]);
    int nbitperelem =atoi(argv[2]);
 
    bitArraySet bas (nelem,nbitperelem);
    
//    bas.set_i(0,298);
    
//    cout<<bas.get_i(0)<<endl;
    
    
//    uint64_t max=(((uint64_t)1)<<nbitperelem)-1;
//    if (max<nelem) nelem=max;
//    for (int i=0; i<nelem; i++) {
//        bas.set_i(i, i+1);
//        
//        
//        for (int j=0; j<=i; j++) {
//            //        cout<<"-------------"<<i<<","<<j<<endl;
//            uint64_t res=bas.get_i(j);
//            if(res != j+1){
//                
//                cerr<<"FILL PB indice"<<j<<" should be "<<j+1<<" while it is "<<res<<" i = "<<i<<endl;
//                return EXIT_FAILURE;
//            }
//            
//            
//        }
//    }
    
    
    
//
//    cout<<"####################################"<<endl;
//    cout<<"TESTING PROBABILISTIC SET"<<endl;
//    cout<<"####################################"<<endl;
//    
    probabilisticSet ps (nelem,nbitperelem);
    
//    ps.add(test_element,test_element);
    
//    cout<<ps.exists(test_element,test_element)<<endl;
    
    // Fill the full vector array (with key = index for now)
    auto start = std::chrono::high_resolution_clock::now();
    for(uint64_t i=0;i<nelem;i++) {
        ps.add(i,i);
    }
    auto finish = std::chrono::high_resolution_clock::now();
    
    auto time_fill = std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count()/nelem;
//    std::cout<<"FILL: " << std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count()/nelem<<" ns per element"<<endl;
    
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
    uint64_t nb_queries = 100000;
    
    

    start = std::chrono::high_resolution_clock::now();
    for(uint64_t i=0;i<nb_queries;i++) {
        if(ps.exists(i%nelem,std::rand())) nb_FP++;
    }

//    cout<<nb_FP<<" FP among "<<nb_queries<<" queries "<<(100*nb_FP)/(float)(nb_queries)<<" %"<<endl;
    finish = std::chrono::high_resolution_clock::now();
//    std::cout <<"QUERY: "<<std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count()/nb_queries<<" ns per query"<<endl;
    auto time_query = std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count()/nb_queries;
    auto FP_rate = (100*nb_FP)/(float)(nb_queries);
    
    cout<<nbitperelem<<" "<<time_fill<<" "<<time_query<<" "<<FP_rate<<endl;
    return EXIT_SUCCESS;
}

