lwacm
=====

Link-wise Artificial Compressibility Method

by Yifan Yang with supervisor Thomas Zeiser @ FAU
                                       01.12.2014

The source code is free, do whatever you want.



New questions:

1, is it a good practice to merge the formulas to speed up calculation?
   e.g change "3.0*1.0/36.0" into "1.0/12.0"




________________________________________________________________________________________________

change log:

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






