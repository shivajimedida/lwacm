lwacm
=====

Link-wise Artificial Compressibility Method C implementation

by Yifan Yang with supervisor Thomas Zeiser @ FAU   01.12.2014

The source code is free, do whatever you want.

________________________________________________________________________________________________

>> project description
   
   the goal of this project is to implement the lwacm in C and optmize it 

>> file description
   
   1, lwacm.c : original source code, after the program teminates, the performance date will be 
      written into file 'log', with the format as: [domain_size]  [run_time]  [MLUps]
      
      e.g  '500000.000000    0.398274    1.255417'
      
   2, lwacm_raw : soure code without domain definition
   
   3, code_gen.c : the code generator that will add the domain size to 'lwacm_raw', and generate
      a new source file named 'lwacm_run_[domain_size].c' ([domain_size] is a integer).
      After that, a 'run.sh' will be generated to run all excutables one by one.

>> run guide:

        [command]          [description]
        
   1,   make               compile code_gen
   2,   ./code_gen         generate all source code for different domain size and run.sh
   3,   make               compile all source codes
   4,   chmod +x run.sh    make the script executable
   5,   ./run.sh           run all programs
   
   
   wait until all program terminates, and all performance data can be found in 'log'.
   
   Note: data in file 'log' are not ordered. to order it, we can use:
   
   6,   sort -n +0 -1 log > log2
   7,   make clean         this will delete all generated files except 'log' 
   
________________________________________________________________________________________________

change log:

v0.4
1, add the code generator to generate the source code for different domain size
2, modify the make file to compile all source codes
3, add run.sh script to excute all excutables
4, add MLUps measurement for the code

v0.3
1, fix the boundary condition code, now the sun of p is stable all the time


v0.23
1, get rid of the store step in mesh loop, which is redundant since we use toggle flag

v0.22
1, use toggle flag to speed up process, and merge the stroe loop into calculation loop
2, add time measurement
3, add MLUps messurement function


v0.21
1, fix some typo
2, add test() function to show sum of rho() and f()

v0.2
1, fix the boundary problem by adding 1 to all x, y and z index for p[] and u[].
2, seperate calculation and store loop to fix the bug of modifing p[0][] while it is still in use.
3, delete f_e[] and f_e_o[] array, use one variable instead.
4, hard code compute_p() and compute_u() functions.

v0.1
1, get rid of all if statements.
2, instead of using alpha loop, hard coded call back functions to speed up process.


Old questioins:

1, equatioin 11, what is f(a,e,o)(x,t) and how to calculate it?
   it is the current value of E.q 10

2, why do we separate the loop for mesh point x into two different loops?
   we can also merge it

3, do we store all the intermediata data into the stime step t ?
   not all of them just 2, current time step and previoud time step
   e.g. double p_store[t][N_X][N_Y][N_Z]  double u[t][N_X][N_Y][N_Z][3];
   
4, it is a good idea to declare variables as globle?
   you can do it

5, with the initial values of p and u provided, why the values of array f dropped dramatically after the 3rd time step ?
   because you have illegal access at index "x-xi_alpha"

6, seem that we do not need 2 arraies to store the values of p and u, we can just use e.g. double p_store[N_X][N_Y][N_Z]
   because the calculation of p and u does not depend on the previous values ?
   No, because the calculation of f[] still depends on p[0][]

7, is it a good practice to merge the formulas to speed up calculation? e.g change "3.0*1.0/36.0" into "1.0/12.0"
   No worries, the compiler will do this instead

8, if we need to copy data from p[1][] to p[0][], how about just use *memcpy() *directly?
   memcpy is evil, too. It causes the same amount of data (memory) traffic as plain C statements.

9, add time measurement before and after the iteration loop and print the number of lattice site
   updates per second, i.e. (N_X*N_Y*N_Z*T_MAX) / Delta_t. A nice graph also would be
   performance as function of the domain size.



