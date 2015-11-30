# BooPHF
BooPHF is a simple library for building minimal perfect hash function.
It is designed to handle large scale datasets. The function is larger than other state-of-the-art libraries, it takes approximately 6.5 bits / elements (compared to 2.62 bits/elem for the emphf lib), but construction is faster and does not require additional memory. 

It is easy to include in other projects (just include a single .h file) and has no dependencies.


# How to run test

A sample usage is provided in file bootest.cpp, compiled with :
    make
    
    ./Bootest 100000000



Construction is currently sequential,  parallel version is under work.

# Author
Guillaume Rizk guillaume.rizk@algorizk.com
