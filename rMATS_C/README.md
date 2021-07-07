Bioinformatics Code Transplantation
===================================

Dependencies
------------

- GSL - GNU Scientific Library

It should be convenient to install GSL on most of the platform. Please follow the installation instruction from gsl's page.

You can try these instruction if you're using Ubuntu,

    sudo apt-get install libgsl0-dev
    sudo apt-get install libgsl0ldbl
    
or using Mac,

    brew install gsl

Build
-----

    make

Usage
-----

Suppose we have an input file ./test.txt:

    ./rMATSexe -i ./test.txt -t 1 -o ./test_output.txt -c 0.001
    
- -t the number of thread.
- -o output filename.
- -i input filename.
- -c cutoff.

Detail
------------

- Simplify some formula in order to avoid round-off error and speed up the program.
- Using GSL to manipulate vector.
- GSL provide a collection of numerical routine for scientific computation.
