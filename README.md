# invertor
For Blockwise Matrix Inversion

By treating the matrix in block form, we perform the inversion using different ways of memory handling.
The programs are in C language.
An algorithm for parallel processing using blockwise inversion for large block partitioned matrix is developed and implemented using OpenMP API.  This parallel code is also available in this repository.

Summary of programs included.
File 1: 'test_invertor.c' - Sample program to test the written other invertor functions.
File 2: 'sampleoutput_inplace_by_a.out' - Sample output generated by running the 'test_invertor.e' after including 'test_invertor_inplace_by_a.c' file.
File 3: 'invertor_by_a.c' - Program performs inversion for partitioned matrix where block A and its Schur complements are invertible.
File 4: 'invertor_inplace_by_a.c' - Program performs inplace inversion for partitioned matrix where block A and its Schur complements are invertible 
File 5: 'invertor_by_ad.c' - Program performs inversion for partitioned matrix where block A, D and their Schur complements are invertible
File 6: 'invertor_by_prll.c' - Program performs inversion for large partitioned block matrix where diagonal blocks and their Schur complements are invertible.

Partitioning Scheme of matrix for File 3, 4 and 5:

Input matrix X = 	A B
			C D
			
============================================================================
Detailed Instruction for running the sample program:
============================================================================
File1: 
'test_invertor.c' - Sample program to test the written invertor functions.

This program is written to test the different functions written in other files.  The user has to uncomment in the 'test_invertor.c' file and include any one of the following list of include files to test the corresponding function.
//#include "invertor_by_a.c"
//#include "invertor_inplace_by_a.c"
//#include "invertor_by_ad.c"
//#include "invertor_by_prll.c"

The 'test_invertor.c' generates a random matrix and inverts it using the different methods.  The order of input matrix can be varied manually editing the 'test_invertor.c' file where the value of n is declared (First line inside the main function).
For eg. 	
	int n=23;

The compilation can be performed using gcc as follows.
	gcc -o test_invertor.e test_invertor.c -lm
	
For the case of running using with OpenMp, include the invertor_by_prll.c file and comment out the remaining files.
For compilation,
	gcc -o test_invertor.e test_invertor.c -lm -fopenmp
	
-------------------------------------------------------------------------------

File 2: 'sampleoutput_inplace_by_a.out' - Sample output generated by running the 'test_invertor.e' after including 'test_invertor_inplace_by_a.c' file.
This output is generate for random generated matrix of order (23 * 23).  The output file contains, the generated input matrix, inverted matrix and amount of time taken for computing the inverse.  It also shows the verification that input matrix * inverted matrix as identity matrix.

--------------------------------------------------------------------------------

Further details are available in the webpage:
https://www.imsc.res.in/~rtsenthil/invertor.html
