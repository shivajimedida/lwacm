/*
 *  Link-wise Artificial Compressibility Method
 *  by Yifan Yang with supervisor Thomas Zeiser @ FAU
 *                                             01.12.2014
 *   xi velocities in D3Q19 stencil, 19 total
 *   
 *   xi(0)  {  0,  0,  0, }
 *   
 *   xi(1)  {  1,  0,  0, }
 *   xi(2)  { -1,  0,  0, }
 *   xi(3)  {  0,  1,  0, }
 *   xi(4)  {  0, -1,  0, }
 *   xi(5)  {  0,  0,  1, }
 *   xi(6)  {  0,  0, -1, }
 *   
 *   xi(7)  {  1,  1,  0, }
 *   xi(8)  { -1,  1,  0, }
 *   xi(9)  {  1, -1,  0, }
 *   xi(10) { -1, -1,  0, }
 *   
 *   xi(11) {  1,  0,  1, }
 *   xi(12) { -1,  0,  1, }
 *   xi(13) {  1,  0, -1, }
 *   xi(14) { -1,  0, -1, }
 *   
 *   xi(15) {  0,  1,  1, }
 *   xi(16) {  0, -1,  1, }
 *   xi(17) {  0,  1, -1, }
 *   xi(18) {  0, -1, -1, }
 */

#include <stdio.h>
#include <time.h>
//#include <unistd.h>

#define T_MAX     200
#define N_X       10
#define N_Y       10
#define N_Z       10
#define ALPHA_MAX 10

//collision frequency
const double omega = 1.5;

// time step and index variables
int t = 0;

int x = 0;
int y = 0;
int z = 0;

int i = 0;

// toggle flag
int t_now = 0;
int t_next = 0;

//time parameter
clock_t start_t, loop_t, end_t, total_t;

// array to store value of p, t=0 is t, t=1 is t+1
// use x+1, y+1, z+1 to avoid illegal access like p(x-xi_alpha), same for u[]
// 
//       t |  x  |  y   |  z
double p[2][N_X+1][N_Y+1][N_Z+1];

// array to store value u
//       t |  x  |  y   |  z    |u_xyz, t=0 is t, t=1 is t+1
double u[2][N_X+1][N_Y+1][N_Z+1][3];

// array to store f
double f[19]     = { 0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0 };

// intermediate variables for calculating f[]
double u_xi = 0;    // u(x-xi) * xi
double u_x_xi = 0;  // u(x) * xi
double u_2 = 0;     // u square

double f_e = 0;     // f(e)
double f_e_o = 0;   // f(e, o)(x-xi_alpha, t)
double f_x_t = 0;   // f(e, o)(x, t)

// for load p and u in alpha callback function
double p_load = 0;
double u_load[3] = {0, 0, 0};

// this is a debug function
void test()
{
    printf(">> time step = %d \n", t);
    
    loop_t = clock();
    
    total_t = (loop_t - start_t) / CLOCKS_PER_SEC;
    
    printf("   time taken by CPU: %f seconds\n\n", (double)total_t  );
    
    // part to sum up rho globally
    double result = 0;
      
    for( x = 1; x < N_X; x++)
    {
      for( y = 1; y < N_Y; y++)
      {
        for( z = 1; z < N_Z; z++)
        {
            // step 13, store p(x, t+1)
            result += p[t_now][x][y][z];
            
        }
      }
    }
    
    printf("   sum of all rho[x][y][z] = %e \n", result );
    
    // print out f[]
    for( i = 0; i < 19; i++)
    {
        printf("   f(%d) = %e, ", i, f[i] );
    }
    
    printf("\n_______________________________________________________________________________________________\n");
}

// call back functions to calculate f for each alpha from 0 to 18
void alpha_0_call()
{      
      // alpha = 0 calculate u * xi and u square  xi(0)  {  0,  0,  0, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x][y][z];
             
      u_load[0] = u[t_now][x][y][z][0];
      u_load[1] = u[t_now][x][y][z][1];
      u_load[2] = u[t_now][x][y][z][2];
      
      u_xi = 0;
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = 0;
      
      //  the follwoing code can be abbreviated
      // step 5, compute f(e)(x-xi[0], t)   
      f_e = 1.0/3.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/3.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/3.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9
      f[0] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_1_call()
{
      // alpha = 1 calculate u * xi and u square  xi(1)  {  1,  0,  0, }
      
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x-1][y][z];
      
      u_load[0] = u[t_now][x-1][y][z][0];
      u_load[1] = u[t_now][x-1][y][z][1];
      u_load[2] = u[t_now][x-1][y][z][2];
      
      u_xi = u_load[0];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0];
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e = 1.0/18.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/18.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/18.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[1] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_2_call()
{
      // alpha = 0 calculate u * xi and u square  xi(2)  { -1,  0,  0, }
      
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x+1][y][z];
      
      u_load[0] = u[t_now][x+1][y][z][0];
      u_load[1] = u[t_now][x+1][y][z][1];
      u_load[2] = u[t_now][x+1][y][z][2];
      
      u_xi = u_load[0]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e = 1.0/18.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/18.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/18.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[2] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_3_call()
{
      // alpha = 0 calculate u * xi and u square  xi(3)  {  0,  1,  0, }
      
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x][y-1][z];
      
      u_load[0] = u[t_now][x][y-1][z][0];
      u_load[1] = u[t_now][x][y-1][z][1];
      u_load[2] = u[t_now][x][y-1][z][2];
      
      u_xi = u_load[1];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][1];
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e = 1.0/18.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/18.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/18.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[3] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_4_call()
{
      // alpha = 0 calculate u * xi and u square  xi(4)  {  0, -1,  0, }
      
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x][y+1][z];
      
      u_load[0] = u[t_now][x][y+1][z][0];
      u_load[1] = u[t_now][x][y+1][z][1];
      u_load[2] = u[t_now][x][y+1][z][2];
      
      u_xi = u_load[1]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][1] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e = 1.0/18.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/18.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/18.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[4] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_5_call()
{
      // alpha = 0 calculate u * xi and u square  xi(5)  {  0,  0,  1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x][y][z-1];
      
      u_load[0] = u[t_now][x][y][z-1][0];
      u_load[1] = u[t_now][x][y][z-1][1];
      u_load[2] = u[t_now][x][y][z-1][2];
      
      u_xi = u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][2];
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e = 1.0/18.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/18.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/18.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[5] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_6_call()
{
      // alpha = 0 calculate u * xi and u square  xi(6)  {  0,  0, -1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x][y][z+1];
      
      u_load[0] = u[t_now][x][y][z+1][0];
      u_load[1] = u[t_now][x][y][z+1][1];
      u_load[2] = u[t_now][x][y][z+1][2];
      
      u_xi = u_load[2]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][2] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e = 1.0/18.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/18.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/18.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[6] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_7_call()
{
      // alpha = 0 calculate u * xi and u square  xi(7)  {  1,  1,  0, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x-1][y-1][z];
      
      u_load[0] = u[t_now][x-1][y-1][z][0];
      u_load[1] = u[t_now][x-1][y-1][z][1];
      u_load[2] = u[t_now][x-1][y-1][z][2];
      
      u_xi = u_load[0] + u_load[1];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] + u[t_now][x][y][z][1];
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[7] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_8_call()
{
      // alpha = 0 calculate u * xi and u square  xi(8)  { -1,  1,  0, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x+1][y-1][z];
      
      u_load[0] = u[t_now][x+1][y-1][z][0];
      u_load[1] = u[t_now][x+1][y-1][z][1];
      u_load[2] = u[t_now][x+1][y-1][z][2];
      
      u_xi = u_load[0]*(-1) + u_load[1];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] * (-1) + u[t_now][x][y][z][1];
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[8] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_9_call()
{
      // alpha = 0 calculate u * xi and u square  xi(9)  {  1, -1,  0, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x-1][y+1][z];
      
      u_load[0] = u[t_now][x-1][y+1][z][0];
      u_load[1] = u[t_now][x-1][y+1][z][1];
      u_load[2] = u[t_now][x-1][y+1][z][2];
      
      u_xi = u_load[0] + u_load[1]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] + u[t_now][x][y][z][1] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[9] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_10_call()
{
      // alpha = 0 calculate u * xi and u square  xi(10) { -1, -1,  0, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x+1][y+1][z];
      
      u_load[0] = u[t_now][x+1][y+1][z][0];
      u_load[1] = u[t_now][x+1][y+1][z][1];
      u_load[2] = u[t_now][x+1][y+1][z][2];
      
      u_xi = u_load[0]*(-1) + u_load[1]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] * (-1) + u[t_now][x][y][z][1] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[10] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_11_call()
{
      // alpha = 0 calculate u * xi and u square  xi(11) {  1,  0,  1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x-1][y][z-1];
      
      u_load[0] = u[t_now][x-1][y][z-1][0];
      u_load[1] = u[t_now][x-1][y][z-1][1];
      u_load[2] = u[t_now][x-1][y][z-1][2];
      
      u_xi = u_load[0] + u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] + u[t_now][x][y][z][2];
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[11] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_12_call()
{
      // alpha = 0 calculate u * xi and u square  xi(12) { -1,  0,  1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x+1][y][z-1];
      
      u_load[0] = u[t_now][x+1][y][z-1][0];
      u_load[1] = u[t_now][x+1][y][z-1][1];
      u_load[2] = u[t_now][x+1][y][z-1][2];
      
      u_xi = u_load[0]*(-1) + u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] * (-1) + u[t_now][x][y][z][2];
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[12] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_13_call()
{
      // alpha = 0 calculate u * xi and u square  xi(13) {  1,  0, -1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x-1][y][z+1];
      
      u_load[0] = u[t_now][x-1][y][z+1][0];
      u_load[1] = u[t_now][x-1][y][z+1][1];
      u_load[2] = u[t_now][x-1][y][z+1][2];
      
      u_xi = u_load[0] + u_load[2]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] + u[t_now][x][y][z][2] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[13] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_14_call()
{
      // alpha = 0 calculate u * xi and u square  xi(14) { -1,  0, -1, }
      
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x+1][y][z+1];
      
      u_load[0] = u[t_now][x+1][y][z+1][0];
      u_load[1] = u[t_now][x+1][y][z+1][1];
      u_load[2] = u[t_now][x+1][y][z+1][2];
      
      u_xi = u_load[0]*(-1) + u_load[2]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][0] * (-1) + u[t_now][x][y][z][2] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[14] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_15_call()
{
      // alpha = 0 calculate u * xi and u square  xi(15) {  0,  1,  1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x][y-1][z-1];
      
      u_load[0] = u[t_now][x][y-1][z-1][0];
      u_load[1] = u[t_now][x][y-1][z-1][1];
      u_load[2] = u[t_now][x][y-1][z-1][2];
      
      u_xi = u_load[1] + u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][1] + u[t_now][x][y][z][2];
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[15] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_16_call()
{
      // alpha = 0 calculate u * xi and u square  xi(16) {  0, -1,  1, }
      
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x][y+1][z-1];
      
      u_load[0] = u[t_now][x][y+1][z-1][0];
      u_load[1] = u[t_now][x][y+1][z-1][1];
      u_load[2] = u[t_now][x][y+1][z-1][2];
      
      u_xi = u_load[1]*(-1) + u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][1] * (-1) + u[t_now][x][y][z][2];
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[16] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_17_call()
{
      // alpha = 0 calculate u * xi and u square  xi(17) {  0,  1, -1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x][y-1][z+1];
      
      u_load[0] = u[t_now][x][y-1][z+1][0];
      u_load[1] = u[t_now][x][y-1][z+1][1];
      u_load[2] = u[t_now][x][y-1][z+1][2];
      
      u_xi = u_load[1] + u_load[2] * (-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][1] + u[t_now][x][y][z][2] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[17] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void alpha_18_call()
{
      // alpha = 0 calculate u * xi and u square  xi(18) {  0, -1, -1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x][y+1][z+1];
      
      u_load[0] = u[t_now][x][y+1][z+1][0];
      u_load[1] = u[t_now][x][y+1][z+1][1];
      u_load[2] = u[t_now][x][y+1][z+1][2];
      
      u_xi = u_load[1] * (-1) + u_load[2] * (-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][z][1] * (-1) + u[t_now][x][y][z][2] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[18] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}


int main()
{
    printf("\n>> lwacm start...\n");
    printf("\n   domain size : N_X = %d, N_Y = %d, N_Z = %d\n", N_X, N_Y, N_Z );
    printf("\n   timestep required : T_MAX = %d\n\n", T_MAX );
    
    
    
    // initialize p and u
    for( x = 0; x < N_X+1; x++)
    {
      for( y = 0; y < N_Y+1; y++)
      {
        for( z = 0; z < N_Z+1; z++)
        {
            p[0][x][y][z] = 1.0f; // set it to 1.0
            p[1][x][y][z] = 1.0f; // set it to 1.0
            
            u[0][x][y][z][0] = 0;
            u[0][x][y][z][1] = 0;
            u[0][x][y][z][2] = 0;
            
            u[1][x][y][z][0] = 0;
            u[1][x][y][z][1] = 0;
            u[1][x][y][z][2] = 0;
        }
      }
    }
    
    
    // initialize toggle flag
    t_now = 0;
    t_next = 1;
    
    
    start_t = clock();
    
    // for all time step t
    for( t = 0; t < T_MAX; t++)
    {
        
        // loop 1, calculation
        // for all mesh point x, calculate p[t_next][] and u[t_next][]
        for( x = 1; x < N_X; x++)
        {
          for( y = 1; y < N_Y; y++)
          {
            for( z = 1; z < N_Z; z++)
            {
                
                // for all alpha from 0 to 18, calculate f[]
                alpha_0_call();
                alpha_1_call();
                alpha_2_call();
                alpha_3_call();
                alpha_4_call();
                alpha_5_call();
                alpha_6_call();
                alpha_7_call();
                alpha_8_call();
                alpha_9_call();
                alpha_10_call();
                alpha_11_call();
                alpha_12_call();
                alpha_13_call();
                alpha_14_call();
                alpha_15_call();
                alpha_16_call();
                alpha_17_call();
                alpha_18_call();
                
                // step 11, compute p(x, t+1)
                p_load = f[0]+f[1]+f[2]+f[3]+f[4]+f[5]+f[6]+f[7]+f[8]+f[9]+f[10]+f[11]+f[12]+f[13]+f[14]+f[15]+f[16]+f[17]+f[18];
                
                p[t_next][x][y][z] = p_load;
                
                // step 12, compute u(x, t+1)
                u[t_next][x][y][z][0] = ( f[1]-f[2]+f[7]-f[8]+f[9]-f[10]+f[11]-f[12]+f[13]-f[14] )/p_load;
                u[t_next][x][y][z][1] = ( f[3]-f[4]+f[7]+f[8]-f[9]-f[10]+f[15]-f[16]+f[17]-f[18] )/p_load;
                u[t_next][x][y][z][2] = ( f[5]-f[6]+f[11]+f[12]-f[13]-f[14]+f[15]+f[16]-f[17]-f[18] )/p_load;
                
                // step 13, store p(x, t+1)
                p[t_now][x][y][z] = p[t_next][x][y][z];
                
                // step 13, store u(x, t+1)
                u[t_now][x][y][z][0] = u[t_next][x][y][z][0];
                u[t_now][x][y][z][1] = u[t_next][x][y][z][1];
                u[t_now][x][y][z][2] = u[t_next][x][y][z][2];
                
                
            }
          }
        }
        
        // toggle t_now and t_next
        t_now  = 1-t_now;
        t_next = 1-t_next; 
        
        test();
        
        //usleep(50000);
        
        
    }
    
    
    end_t = clock();
    
    total_t = (end_t - start_t) / CLOCKS_PER_SEC;
    
    printf("\n   Total loop time taken by CPU: %f seconds\n", (double)total_t  );
    printf("   MLUps =  %f \n\n", (N_X * N_Y * N_Z * T_MAX)/(double)total_t/1000000.0  );
    
    return 0;
}




