# BooPHF
BooPHF is a simple library for building minimal perfect hash function.
It is designed to handle large scale datasets. The function is larger than other state-of-the-art libraries, it takes approximately 6.5 bits / elements (compared to 2.62 bits/elem for the emphf lib), but construction is faster and does not require additional memory. 

It is easy to include in other projects (just include a single .h file) and has no dependencies.


# How to run test

A sample usage is provided in file bootest.cpp, compile and run with :

    make
    ./Bootest 100000000
    
Use -check to check correctness of the hash function, and -bench to benchmark lookup performance.
    
Here is a sample output :
    
    $./Bootest 100000000  -check -bench
    Construct a BooPHF with  100000000 elements (762 MB for holding elems in ram)
    [Building BooPHF]  100  %   elapsed:   1 min 50 sec   remaining:   0 min 0  sec
    BooPHF constructed perfect hash for 100000000 keys in 109.84s
    Bitarray       330133248  bits (54.24 %)   (array + ranks )
    Blooms         267298408  bits (43.91 %)
    final hash      11244576  bits (1.85 %) (nb in final hash 33466)
    boophf  bits/elem : 6.086762
     --- boophf working correctly --- 
    BooPHF 100000000 lookups in  26.39s,  approx  263.90 ns per lookup   (fingerprint 4999999950000000) 
    


Construction is currently sequential,  parallel version is under work.

# Author
Guillaume Rizk guillaume.rizk@algorizk.com
