# BBHash
BBHash is a simple library for building minimal perfect hash function.
It is designed to handle large scale datasets. The function is just a little bit larger than other state-of-the-art libraries, it takes approximately 3 bits / elements (compared to 2.62 bits/elem for the emphf lib), but construction is faster and does not require additional memory. 

It is easy to include in other projects (just include a single .h file) and has no dependencies.


# How to run test

A sample usage is provided in file bootest.cpp, compile and run with : ( params are nb_elements nb_threads)

    make
    ./Bootest 100000000 1
    
Use -check to check correctness of the hash function, and -bench to benchmark lookup performance.
    
Here is a sample output :
    
    $./Bootest 100000000 8 -check -bench
    key file generated 
    Construct a BooPHF with  100000000 elements  
    [Building BooPHF]  100  %   elapsed:   0 min 31 sec   remaining:   0 min 0  sec
    BooPHF constructed perfect hash for 100000000 keys in 30.95s
    Bitarray       305808384  bits (99.49 %)   (array + ranks )
    final hash       1557024  bits (0.51 %) (nb in final hash 4634)
    boophf  bits/elem : 3.073654
     --- boophf working correctly --- 
    bench lookups  sample size 9999872 
    BBhash bench lookups average 258.06 ns +- stddev  5.55 %   (fingerprint 5000111281850410) 



# Author
Guillaume Rizk guillaume.rizk@algorizk.com
