//
//  test.cpp
//  bank_request
//
//  Created by Pierre Peterlongo on 01/03/16.
//
//
#include <iostream>
#include "bit_vector_array.h"

using namespace std;

int main (int argc, char* argv[]){
    
    if(argc<3){
		printf("Usage :\n");
		printf("%s <nelem> <nbitperelem>  [options]\n",argv[0]);
		return EXIT_FAILURE;
        
    }
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
    return EXIT_SUCCESS;

    
}

