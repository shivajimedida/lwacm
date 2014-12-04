lwacm
=====

Link-wise Artificial Compressibility Method


questioins:

1, equatioin 11, what is f(a,e,o)(x,t) and how to calculate it: it is the current value of E.q 10

2, why do we separate the loop for mesh point x into two different loops? we can also merge it

3, do we store all the intermediata data into the stime step t ?           not all of them just 2, current time step and previoud time step
   e.g. double p_store[t][N_X][N_Y][N_Z]  double u[t][N_X][N_Y][N_Z][3];
   
4, it is a good idea to declare variables as globle? you can do it




f[0] = 1.0/3.0 * p(x-xi0) * ( 1 + 3 )







