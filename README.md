lwacm
=====

Link-wise Artificial Compressibility Method

update:

1, get rid of all if statements.
2, instead of using alpha loop, I used hard coded call back functions to speed up process.


questions:

1, with the initial values of p and u provided, why the values of array f dropped dramatically after the 3rd time step ?

2, seem that we do not need 2 arraies to store the values of p and u, we can just use e.g. double p_store[N_X][N_Y][N_Z]
   because the calculation of p and u does not depend on the previous values ?



old questioins:

1, equatioin 11, what is f(a,e,o)(x,t) and how to calculate it: it is the current value of E.q 10

2, why do we separate the loop for mesh point x into two different loops? we can also merge it

3, do we store all the intermediata data into the stime step t ?           not all of them just 2, current time step and previoud time step
   e.g. double p_store[t][N_X][N_Y][N_Z]  double u[t][N_X][N_Y][N_Z][3];
   
4, it is a good idea to declare variables as globle? you can do it








