/*  v 0.4
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
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

// initilize domain size and time steps
int T_MAX = 0;
int N_X = 0;
int N_Y = 0;
int N_Z = 0;

//collision frequency
const double omega = 1.5;

// time step and index variables
int t = 0;

int x = 0;
int y = 0;
int z = 0;

int i = 0;

long long bytes = 0;

// toggle flag
int t_now = 0;
int t_next = 0;

//time parameter
struct timeval tm_start, tm_end;


//time_t begin, end;
//clock_t start_t = 0;
//clock_t end_t = 0;

double sec_elapsed = 0;

//performance parameter
double domain_size = 0;
double total_lattice_update = 0;
double mlups = 0;


// array to store value of p, t=0 is t, t=1 is t+1
// use x+2, y+2, z+2 to avoid illegal access like p(x-xi_alpha), same for u[]
// 
// array format
//          t |  x  |  y   |  z
// double p[2][N_X+2][N_Y+2][N_Z+2];

// 4 dimensional array pointer for dynamic memory allocation
double ****p = NULL;


// array to store value u
//         t |  x  |  y   |  z    |u_xyz, t=0 is t, t=1 is t+1
// double u[2][N_X+2][N_Y+2][3][N_Z+2];

// 5 dimensional array pointer for dynamic memory allocation
double *****u = NULL;


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
void sum_up_p()
{
    fprintf(stderr,"\n   > step = %d ", t);
    
    // part to sum up rho globally
    double result = 0;
      
    for( x = 1; x < N_X+1; x++)
    {
        for( y = 1; y < N_Y+1; y++)
        {
            for( z = 1; z < N_Z+1; z++)
            {
                // calculate sum of p[]
                result += p[t_next][x][y][z];
            }
        }
    }
    
    fprintf(stderr,"   sum of all rho[x][y][z] = %e ", result );
}

// small function to show progress bar if you are tired of waiting
void progress_bar()
{   
    int j = 0;
    
    i = t*100/T_MAX;
    
    fprintf(stderr, "\r");
    
    if(i < 10) { fprintf(stderr, "   [00%d %%] ", i); }    
    else if(i < 100){ fprintf(stderr, "   [0%d %%] ", i); }
    else { fprintf(stderr, "   [%d %%] ", i); }
        
    for(j = 0; j < i; ++j)
    {
        if(j == i-1){ fprintf(stderr, ">"); }
        else{ fprintf(stderr, "="); }
    }
    
    fflush(stdout);
}

// call back functions to calculate f for each alpha from 0 to 18
void alpha_0_call()
{      
      // alpha = 0 calculate u * xi and u square  xi(0)  {  0,  0,  0, }
      
      //load p(x-xia) and u(x-xia)
      p_load = p[t_now][x][y][z];
      
      u_load[0] = u[t_now][x][y][0][z];
      u_load[1] = u[t_now][x][y][1][z];
      u_load[2] = u[t_now][x][y][2][z];
      
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
      
      u_load[0] = u[t_now][x-1][y][0][z];
      u_load[1] = u[t_now][x-1][y][1][z];
      u_load[2] = u[t_now][x-1][y][2][z];
      
      u_xi = u_load[0];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][0][z];
      
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
      
      u_load[0] = u[t_now][x+1][y][0][z];
      u_load[1] = u[t_now][x+1][y][1][z];
      u_load[2] = u[t_now][x+1][y][2][z];
      
      u_xi = u_load[0]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][0][z] * (-1);
      
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
      
      u_load[0] = u[t_now][x][y-1][0][z];
      u_load[1] = u[t_now][x][y-1][1][z];
      u_load[2] = u[t_now][x][y-1][2][z];
      
      u_xi = u_load[1];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][1][z];
      
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
      
      u_load[0] = u[t_now][x][y+1][0][z];
      u_load[1] = u[t_now][x][y+1][1][z];
      u_load[2] = u[t_now][x][y+1][2][z];
      
      u_xi = u_load[1]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][1][z] * (-1);
      
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
      
      u_load[0] = u[t_now][x][y][0][z-1];
      u_load[1] = u[t_now][x][y][1][z-1];
      u_load[2] = u[t_now][x][y][2][z-1];
      
      u_xi = u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][2][z];
      
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
      
      u_load[0] = u[t_now][x][y][0][z+1];
      u_load[1] = u[t_now][x][y][1][z+1];
      u_load[2] = u[t_now][x][y][2][z+1];
      
      u_xi = u_load[2]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][2][z] * (-1);
      
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
      
      u_load[0] = u[t_now][x-1][y-1][0][z];
      u_load[1] = u[t_now][x-1][y-1][1][z];
      u_load[2] = u[t_now][x-1][y-1][2][z];
      
      u_xi = u_load[0] + u_load[1];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][0][z] + u[t_now][x][y][1][z];
      
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
      
      u_load[0] = u[t_now][x+1][y-1][0][z];
      u_load[1] = u[t_now][x+1][y-1][1][z];
      u_load[2] = u[t_now][x+1][y-1][2][z];
      
      u_xi = u_load[0]*(-1) + u_load[1];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][0][z] * (-1) + u[t_now][x][y][1][z];
      
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
      
      u_load[0] = u[t_now][x-1][y+1][0][z];
      u_load[1] = u[t_now][x-1][y+1][1][z];
      u_load[2] = u[t_now][x-1][y+1][2][z];
      
      u_xi = u_load[0] + u_load[1]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][0][z] + u[t_now][x][y][1][z] * (-1);
      
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
      
      u_load[0] = u[t_now][x+1][y+1][0][z];
      u_load[1] = u[t_now][x+1][y+1][1][z];
      u_load[2] = u[t_now][x+1][y+1][2][z];
      
      u_xi = u_load[0]*(-1) + u_load[1]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][0][z] * (-1) + u[t_now][x][y][1][z] * (-1);
      
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
      
      u_load[0] = u[t_now][x-1][y][0][z-1];
      u_load[1] = u[t_now][x-1][y][1][z-1];
      u_load[2] = u[t_now][x-1][y][2][z-1];
      
      u_xi = u_load[0] + u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][0][z] + u[t_now][x][y][2][z];
      
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
      
      u_load[0] = u[t_now][x+1][y][0][z-1];
      u_load[1] = u[t_now][x+1][y][1][z-1];
      u_load[2] = u[t_now][x+1][y][2][z-1];
      
      u_xi = u_load[0]*(-1) + u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][0][z] * (-1) + u[t_now][x][y][2][z];
      
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
      
      u_load[0] = u[t_now][x-1][y][0][z+1];
      u_load[1] = u[t_now][x-1][y][1][z+1];
      u_load[2] = u[t_now][x-1][y][2][z+1];
      
      u_xi = u_load[0] + u_load[2]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][0][z] + u[t_now][x][y][2][z] * (-1);
      
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
      
      u_load[0] = u[t_now][x+1][y][0][z+1];
      u_load[1] = u[t_now][x+1][y][1][z+1];
      u_load[2] = u[t_now][x+1][y][2][z+1];
      
      u_xi = u_load[0]*(-1) + u_load[2]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][0][z] * (-1) + u[t_now][x][y][2][z] * (-1);
      
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
      
      u_load[0] = u[t_now][x][y-1][0][z-1];
      u_load[1] = u[t_now][x][y-1][1][z-1];
      u_load[2] = u[t_now][x][y-1][2][z-1];
      
      u_xi = u_load[1] + u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][1][z] + u[t_now][x][y][2][z];
      
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
      
      u_load[0] = u[t_now][x][y+1][0][z-1];
      u_load[1] = u[t_now][x][y+1][1][z-1];
      u_load[2] = u[t_now][x][y+1][2][z-1];
      
      u_xi = u_load[1]*(-1) + u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][1][z] * (-1) + u[t_now][x][y][2][z];
      
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
      
      u_load[0] = u[t_now][x][y-1][0][z+1];
      u_load[1] = u[t_now][x][y-1][1][z+1];
      u_load[2] = u[t_now][x][y-1][2][z+1];
      
      u_xi = u_load[1] + u_load[2] * (-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][1][z] + u[t_now][x][y][2][z] * (-1);
      
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
      
      u_load[0] = u[t_now][x][y+1][0][z+1];
      u_load[1] = u[t_now][x][y+1][1][z+1];
      u_load[2] = u[t_now][x][y+1][2][z+1];
      
      u_xi = u_load[1] * (-1) + u_load[2] * (-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[t_now][x][y][1][z] * (-1) + u[t_now][x][y][2][z] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[t_now][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[18] =  f_e + 2*( omega-1/omega )*( f_x_t - f_e_o ) ;  // Eq.11
}

void free_array_p()
{
    fprintf(stderr, "\n   > free array p...    ");
    
    // free mem for z
    for(t = 0; t < 2; t++)
    {
        for(x = 0; x < N_X+2; x++)
        {
            for(y = 0; y < N_Y+2; y++)
            {
                free( p[t][x][y] );
            }  
        }
    }
    
    // free mem for y
    for(t = 0; t < 2; t++)
    {
        for(x = 0; x < N_X+2; x++)
        {
            free( p[t][x] );
        }
    }
    
    // free mem for x
    for(t = 0; t < 2; t++)
    {
        free( p[t] );
    }
    
    // free data pointer
    free(p);
    
    fprintf(stderr, "done\n");
    
}

void free_array_u()
{
    fprintf(stderr, "\n   > free array u...    ");
    
    // free mem for z
    for(t = 0; t < 2; t++)
    {
        for(x = 0; x < N_X+2; x++)
        {
            for(y = 0; y < N_Y+2; y++)
            {
                for(i = 0; i < 3; i++)
                {
                    free( u[t][x][y][i] );
                }
            }
        }
    }
    
    // free mem for u_xyz
    for(t = 0; t < 2; t++)
    {
        for(x = 0; x < N_X+2; x++)
        {
            for(y = 0; y < N_Y+2; y++)
            {
                free( u[t][x][y] );
            }  
        }
    }
    
    // free mem for y
    for(t = 0; t < 2; t++)
    {
        for(x = 0; x < N_X+2; x++)
        {
            free( u[t][x] );
        }
    }
    
    // free mem for x
    for(t = 0; t < 2; t++)
    {
        free( u[t] );
    }
    
    // free data pointer
    free(u);
    
    fprintf(stderr, "done\n");
    
}



void useage()
{
    fprintf(stderr, "\n>> lwacm, a benchmark tool for LWACM kernel\n");
    fprintf(stderr, "   usage: lwacm <flag> <parameter>\n");
    fprintf(stderr, "          flag: s <parameter:int>, specify the domain size\n");
    fprintf(stderr, "          flag: t <parameter:int>, specify the max time step\n");
    fprintf(stderr, "\n          e.g ./lwacm s 100 t 20\n\n");
}


int main( int argc, char *argv[] )
{
    FILE *fileout = NULL;
    fileout = fopen("log", "a+");
    
    // argument parsing
    if(argc < 5) {
        fprintf(stderr, "\n>> parameter not specified, please check\n");
        useage();
        return 1;
    }
    
    else
    {
        for (i = 1; i < argc; ++i)
        {
            if(  *argv[i] == 'h'  )
            {
                useage();
                return 0;
            }
            
            if(  *argv[i] == 's'  )
            {
                N_X = atoi(argv[i + 1]);
                N_Y = N_X;
                N_Z = N_X;
            }
                
            if(  *argv[i] == 't'  )
            {
                T_MAX = atoi(argv[i + 1]);
            }
        }
    }
    
    if(N_X == 0 || N_Y == 0 || N_Z == 0 || T_MAX == 0) {
        fprintf(stderr, "\n>> domain size or T_MAX cannot be set to 0, please check\n");
        useage();
        return 1;
    }
    
    
    fprintf(stderr, "\n>> program start\n\n   domain size : N_X = %d, N_Y = %d, N_Z = %d\n", N_X, N_Y, N_Z );
    fprintf(stderr, "\n   timestep required : T_MAX = %d\n", T_MAX );
    
    // calculate memory consumption for p
    bytes = 2 * (N_X+2) * (N_Y+2) * (N_Z+2) * sizeof(double);
    
    if( bytes > 1024*1024 ) {
        fprintf(stderr, "\n   size of array p : %lld MBytes\n", bytes/1024/1024);
    }
    
    else if( bytes > 1024 ) {
        fprintf(stderr, "\n   size of array p : %lld KBytes\n", bytes/1024);
    }
    
    else {
        fprintf(stderr, "\n   size of array p :  %lld Bytes\n", bytes);
    }
    
    
    // calculate memory consumption for p
    bytes = bytes * 3;
    
    if( bytes > 1024*1024 ) {
        fprintf(stderr, "\n   size of array u : %lld MBytes\n", bytes/1024/1024);
    }
    
    else if( bytes > 1024 ) {
        fprintf(stderr, "\n   size of array u : %lld KBytes\n", bytes/1024);
    }
    
    else {
        fprintf(stderr, "\n   size of array u : %lld Bytes\n", bytes);
    }
    
    
    
    fprintf(stderr, "\n   > start allocating memory for array p...    ");
    
    // allocate mem for array pointer
    p = malloc( 2 * sizeof(double ***) );
    
    if(p == NULL)
    {
        fprintf(stderr, "\n   error allocating memory for array pointer (size = 2)\n\n");
        return 1;
    }
    
    // allocate mem for t_now and t_next
    for(t = 0; t < 2; t++)
    {
        p[t] = malloc( (N_X+2) * sizeof(double**) );
        
        if(p[t] == NULL)
        {
            fprintf(stderr, "\n   error allocating memory for t_now/t_next (size = %d)\n\n", N_X);
            return 1;
        }
        //else { printf("\n   mem [%d] ok\n", i); }
        
        // allocate mem for N_X
        for(x = 0; x < N_X+2; x++)
        {
            p[t][x] = malloc( (N_Y+2) * sizeof(double*) );
        
            if(p[t][x] == NULL)
            {
                fprintf(stderr, "\n   error allocating memory for N_X (size = %d)\n\n", N_Y);
                return 1;
            }
            //else { printf("\n   mem [%d][%d] ok\n", i, x); }
            
            // allocate mem for N_X(actuall data)
            for(y = 0; y < N_Y+2; y++)
            {
                p[t][x][y] = malloc( (N_Z+2) * sizeof(double) );
        
                if(p[t][x][y] == NULL)
                {
                    fprintf(stderr, "\n   error allocating memory for N_Y (size =%d)\n\n", N_Z);
                    return 1;
                }
                //else { printf("\n   mem [%d][%d][%d] ok\n", i, x, y); }
            }
        }
    }
    
    
    fprintf(stderr, "done\n\n   > start allocating memory for array u...    ");
    
    // allocate mem for array pointer
    u = malloc( 2 * sizeof(double ****) );
    
    if(p == NULL)
    {
        fprintf(stderr, "\n   error allocating memory for array pointer(size = 2)\n\n");
        return 1;
    }
    
    // allocate mem for t_now and t_next
    for(t = 0; t < 2; t++)
    {
        u[t] = malloc( (N_X+2) * sizeof(double***) );
        
        if(u[t] == NULL)
        {
            fprintf(stderr, "\n   error allocating memory for t_now/t_next(size = %d)\n\n", N_X);
            return 1;
        }
        //else { printf("\n   mem [%d] ok\n", i); }
        
        // allocate mem for N_X
        for(x = 0; x < N_X+2; x++)
        {
            u[t][x] = malloc( (N_Y+2) * sizeof(double**) );
        
            if(u[t][x] == NULL)
            {
                fprintf(stderr, "\n   error allocating memory for N_X(size = %d)\n\n", N_Y);
                return 1;
            }
            //else { printf("\n   mem [%d][%d] ok\n", i, x); }
            
            // allocate mem for N_X
            for(y = 0; y < N_Y+2; y++)
            {
                u[t][x][y] = malloc( (N_Z+2) * sizeof(double*) );
        
                if(u[t][x][y] == NULL)
                {
                    fprintf(stderr, "\n   error allocating memory for N_Y(size =%d)\n\n", N_Z);
                    return 1;
                }
                //else { printf("\n   mem [%d][%d][%d] ok\n", i, x, y); }
                
                // allocate mem for N_Z
                for(i = 0; i < 3; i++)
                {
                    u[t][x][y][i] = malloc( (N_Z+2) * sizeof(double) );
        
                    if(u[t][x][y][i] == NULL)
                    {
                        fprintf(stderr, "\n   error allocating memory for N_Z(size = 3)\n\n");
                        return 1;
                    }
                }
            }
        }
    }

    fprintf(stderr, "done\n\n   > initializing array p and u...    ");
    
    // initialize p and u
    for( x = 0; x < N_X+2; x++)
    {
      for( y = 0; y < N_Y+2; y++)
      {
        for( z = 0; z < N_Z+2; z++)
        {
            p[0][x][y][z] = 1.0f; // set it to 1.0
            p[1][x][y][z] = 1.0f; // set it to 1.0
            
            u[0][x][y][0][z] = 0;
            u[0][x][y][1][z] = 0;
            u[0][x][y][2][z] = 0;
            
            u[1][x][y][0][z] = 0;
            u[1][x][y][1][z] = 0;
            u[1][x][y][2][z] = 0;
        }
      }
    }
    
    // initialize toggle flag
    t_now = 0;
    t_next = 1;
    
    fprintf(stderr, "done\n\n   > updating cell...    ");
    
    // record the start time
    //time(&begin);
    //start_t = clock();
    
    gettimeofday( &tm_start, NULL );
    
    // for all time step t
    for( t = 0; t < T_MAX; t++)
    {
        // for all mesh point x, calculate p[t_next][] and u[t_next][]
        for( x = 1; x < N_X+1; x++)
        {
          for( y = 1; y < N_Y+1; y++)
          {
            for( z = 1; z < N_Z+1; z++)
            {
                // for all alpha from 0 to 18, calculate f[] for each node
                // test function for alpha
                //for( a = 0; a < 19; a++)
                //{
                //    alpha_call(a);
                //}
                
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
                u[t_next][x][y][0][z] = ( f[1]-f[2]+f[7]-f[8]+f[9]-f[10]+f[11]-f[12]+f[13]-f[14] )/p_load;
                u[t_next][x][y][1][z] = ( f[3]-f[4]+f[7]+f[8]-f[9]-f[10]+f[15]-f[16]+f[17]-f[18] )/p_load;
                u[t_next][x][y][2][z] = ( f[5]-f[6]+f[11]+f[12]-f[13]-f[14]+f[15]+f[16]-f[17]-f[18] )/p_load;
            }
          }
        }
        
        //progress_bar();
        //sum_up_p();
        
        //######################################################################################################################
        //##################################################### start ##########################################################
        // explicit peridic boundary conditions
        // surfaces: xy
        for( x = 1; x < N_X+1; x++)
        {
            for( y = 1; y < N_Y+1; y++)
            {
                p[t_next][ x ][ y ][ 0 ] = p[t_next][ x ][ y ][N_Z];
                u[t_next][ x ][ y ][ 0 ][ 0 ] = u[t_next][ x ][ y ][ 0 ][N_Z];
                u[t_next][ x ][ y ][ 1 ][ 0 ] = u[t_next][ x ][ y ][ 1 ][N_Z];
                u[t_next][ x ][ y ][ 2 ][ 0 ] = u[t_next][ x ][ y ][ 2 ][N_Z];
          
                p[t_next][ x ][ y ][N_Z+1] = p[t_next][ x ][ y ][1];
                u[t_next][ x ][ y ][ 0 ][N_Z+1] = u[t_next][ x ][ y ][ 0 ][ 1 ];
                u[t_next][ x ][ y ][ 1 ][N_Z+1] = u[t_next][ x ][ y ][ 1 ][ 1 ];
                u[t_next][ x ][ y ][ 2 ][N_Z+1] = u[t_next][ x ][ y ][ 2 ][ 1 ];
          
            }
        }
        
        // surfaces: xz
        for( x = 1; x < N_X+1; x++)
        {
            for( z = 1; z < N_Z+1; z++)
            {
                p[t_next][ x ][ 0 ][ z ] = p[t_next][ x ][N_Y][ z ];
                u[t_next][ x ][ 0 ][ 0 ][ z ] = u[t_next][ x ][N_Y][ 0 ][ z ];
                u[t_next][ x ][ 0 ][ 1 ][ z ] = u[t_next][ x ][N_Y][ 1 ][ z ];
                u[t_next][ x ][ 0 ][ 2 ][ z ] = u[t_next][ x ][N_Y][ 2 ][ z ];
          
                p[t_next][ x ][N_Y+1][ z ] = p[t_next][ x ][ 1 ][ z ];
                u[t_next][ x ][N_Y+1][ 0 ][ z ] = u[t_next][ x ][ 1 ][ 0 ][ z ];
                u[t_next][ x ][N_Y+1][ 1 ][ z ] = u[t_next][ x ][ 1 ][ 1 ][ z ];
                u[t_next][ x ][N_Y+1][ 2 ][ z ] = u[t_next][ x ][ 1 ][ 2 ][ z ];
            }
        }
        
        // surfaces: yz
        for( y = 1; y < N_Y+1; y++)
        {
            for( z = 1; z < N_Z+1; z++)
            {
                p[t_next][ 0 ][ y ][ z ] = p[t_next][N_X][ y ][ z ];
                u[t_next][ 0 ][ y ][ 0 ][ z ] = u[t_next][N_X][ y ][ 0 ][ z ];
                u[t_next][ 0 ][ y ][ 1 ][ z ] = u[t_next][N_X][ y ][ 1 ][ z ];
                u[t_next][ 0 ][ y ][ 2 ][ z ] = u[t_next][N_X][ y ][ 2 ][ z ];
          
                p[t_next][N_X+1][ y ][ z ] = p[t_next][ 1 ][ y ][ z ];
                u[t_next][N_X+1][ y ][ 0 ][ z ] = u[t_next][ 1 ][ y ][ 0 ][ z ];
                u[t_next][N_X+1][ y ][ 1 ][ z ] = u[t_next][ 1 ][ y ][ 1 ][ z ];
                u[t_next][N_X+1][ y ][ 2 ][ z ] = u[t_next][ 1 ][ y ][ 2 ][ z ];
            }
        }
      
      
        // lines x: 4
        for( x = 1; x < N_X+1; x++)
        {
            // (x, 0, 0) 
            p[t_next][ x ][ 0 ][ 0 ] = p[t_next][ x ][N_Y][N_Z];
            u[t_next][ x ][ 0 ][ 0 ][ 0 ] = u[t_next][ x ][N_Y][ 0 ][N_Z];
            u[t_next][ x ][ 0 ][ 1 ][ 0 ] = u[t_next][ x ][N_Y][ 1 ][N_Z];
            u[t_next][ x ][ 0 ][ 2 ][ 0 ] = u[t_next][ x ][N_Y][ 2 ][N_Z];
          
            // (x, N_Y+1, N_Z+1)
            p[t_next][ x ][N_Y+1][N_Z+1] = p[t_next][ x ][ 1 ][ 1 ];
            u[t_next][ x ][N_Y+1][ 0 ][N_Z+1] = u[t_next][ x ][ 1 ][ 0 ][ 1 ];
            u[t_next][ x ][N_Y+1][ 1 ][N_Z+1] = u[t_next][ x ][ 1 ][ 1 ][ 1 ];
            u[t_next][ x ][N_Y+1][ 2 ][N_Z+1] = u[t_next][ x ][ 1 ][ 2 ][ 1 ];
          
            // (x, 0, N_Z+1)
            p[t_next][ x ][ 0 ][N_Z+1] = p[t_next][ x ][N_Y][ 1 ];
            u[t_next][ x ][ 0 ][ 0 ][N_Z+1] = u[t_next][ x ][N_Y][ 0 ][ 1 ];
            u[t_next][ x ][ 0 ][ 1 ][N_Z+1] = u[t_next][ x ][N_Y][ 1 ][ 1 ];
            u[t_next][ x ][ 0 ][ 2 ][N_Z+1] = u[t_next][ x ][N_Y][ 2 ][ 1 ];
          
            // (x, N_Y+1, 0)
            p[t_next][ x ][N_Y+1][ 0 ] = p[t_next][ x ][ 1 ][N_Z];
            u[t_next][ x ][N_Y+1][ 0 ][ 0 ] = u[t_next][ x ][ 1 ][ 0 ][N_Z];
            u[t_next][ x ][N_Y+1][ 1 ][ 0 ] = u[t_next][ x ][ 1 ][ 1 ][N_Z];
            u[t_next][ x ][N_Y+1][ 2 ][ 0 ] = u[t_next][ x ][ 1 ][ 2 ][N_Z];
        }
        
        // lines y: 4
        for( y = 1; y < N_Y+1; y++)
        {
            // (0, y, 0) 
            p[t_next][ 0 ][ y ][ 0 ] = p[t_next][N_X][y][N_Z];
            u[t_next][ 0 ][ y ][ 0 ][ 0 ] = u[t_next][N_X][ y ][ 0 ][N_Z];
            u[t_next][ 0 ][ y ][ 1 ][ 0 ] = u[t_next][N_X][ y ][ 1 ][N_Z];
            u[t_next][ 0 ][ y ][ 2 ][ 0 ] = u[t_next][N_X][ y ][ 2 ][N_Z];
          
            // (N_X+1, y, N_Z+1)
            p[t_next][N_X+1][ y ][N_Z+1] = p[t_next][ 1 ][ y ][ 1 ];
            u[t_next][N_X+1][ y ][ 0 ][N_Z+1] = u[t_next][ 1 ][ y ][ 0 ][ 1 ];
            u[t_next][N_X+1][ y ][ 1 ][N_Z+1] = u[t_next][ 1 ][ y ][ 1 ][ 1 ];
            u[t_next][N_X+1][ y ][ 2 ][N_Z+1] = u[t_next][ 1 ][ y ][ 2 ][ 1 ];
          
            // (0, y, N_Z+1)
            p[t_next][ 0 ][ y ][N_Z+1] = p[t_next][N_X][ y ][ 1 ];
            u[t_next][ 0 ][ y ][ 0 ][N_Z+1] = u[t_next][N_X][ y ][ 0 ][ 1 ];
            u[t_next][ 0 ][ y ][ 1 ][N_Z+1] = u[t_next][N_X][ y ][ 1 ][ 1 ];
            u[t_next][ 0 ][ y ][ 2 ][N_Z+1] = u[t_next][N_X][ y ][ 2 ][ 1 ];
          
            // (N_X+1, y, 0)
            p[t_next][N_X+1][ y ][ 0 ] = p[t_next][ 1 ][ y ][N_Z];
            u[t_next][N_X+1][ y ][ 0 ][ 0 ] = u[t_next][ 1 ][ y ][ 0 ][N_Z];
            u[t_next][N_X+1][ y ][ 1 ][ 0 ] = u[t_next][ 1 ][ y ][ 1 ][N_Z];
            u[t_next][N_X+1][ y ][ 2 ][ 0 ] = u[t_next][ 1 ][ y ][ 2 ][N_Z];
        }
        
        // lines z: 4
        for( z = 1; z < N_Z+1; z++)
        {
            // (0, 0, z) 
            p[t_next][ 0 ][ 0 ][ z ] = p[t_next][N_X][N_Y][ z ];
            u[t_next][ 0 ][ 0 ][ 0 ][ z ] = u[t_next][N_X][N_Y][ 0 ][ z ];
            u[t_next][ 0 ][ 0 ][ 1 ][ z ] = u[t_next][N_X][N_Y][ 1 ][ z ];
            u[t_next][ 0 ][ 0 ][ 2 ][ z ] = u[t_next][N_X][N_Y][ 2 ][ z ];
          
            // (N_X+1, N_Y+1, z)
            p[t_next][N_X+1][N_Y+1][ z ] = p[t_next][ 1 ][ 1 ][ z ];
            u[t_next][N_X+1][N_Y+1][ 0 ][ z ] = u[t_next][ 1 ][ 1 ][ 0 ][ z ];
            u[t_next][N_X+1][N_Y+1][ 1 ][ z ] = u[t_next][ 1 ][ 1 ][ 1 ][ z ];
            u[t_next][N_X+1][N_Y+1][ 1 ][ z ] = u[t_next][ 1 ][ 1 ][ 2 ][ z ];
          
            // (0, N_Y+1, z)
            p[t_next][ 0 ][N_Y+1][ z ] = p[t_next][N_X][ 1 ][ z ];
            u[t_next][ 0 ][N_Y+1][ 0 ][ z ] = u[t_next][N_X][ 1 ][ 0 ][ z ];
            u[t_next][ 0 ][N_Y+1][ 1 ][ z ] = u[t_next][N_X][ 1 ][ 1 ][ z ];
            u[t_next][ 0 ][N_Y+1][ 2 ][ z ] = u[t_next][N_X][ 1 ][ 2 ][ z ];
          
            // (N_X+1, 0, z)
            p[t_next][N_X+1][ 0 ][ z ] = p[t_next][ 1 ][N_Y][ z ];
            u[t_next][N_X+1][ 0 ][ 0 ][ z ] = u[t_next][ 1 ][N_Y][ 0 ][ z ];
            u[t_next][N_X+1][ 0 ][ 1 ][ z ] = u[t_next][ 1 ][N_Y][ 1 ][ z ];
            u[t_next][N_X+1][ 0 ][ 2 ][ z ] = u[t_next][ 1 ][N_Y][ 2 ][ z ];
        }
        
        //######################################################## end #########################################################
        //######################################################################################################################
        
        // toggle t_now and t_next
        t_now  = 1-t_now;
        t_next = 1-t_next;
    }
    
    gettimeofday( &tm_end, NULL );
    
    fprintf(stderr, "done\n");
    
    // record the end time
    
    sec_elapsed = (tm_end.tv_sec-tm_start.tv_sec) + (tm_end.tv_usec-tm_start.tv_usec)/1000000.0;
    
    //time(&end);
    //sec_real = difftime(end, begin);
    
    //end_t = clock();
    //sec_cpu = (double)(end_t-start_t) / CLOCKS_PER_SEC;
    
    // free dynamically allocated array p and u
    free_array_p();
    free_array_u();
    
    domain_size = N_X * N_Y * N_Z;
    total_lattice_update = N_X * N_Y * N_Z * T_MAX;
    mlups = total_lattice_update/sec_elapsed/1000000.0;
    
    fprintf(stderr, "\n   Wall Clock:  %f seconds\n", sec_elapsed );
    fprintf(stderr, "   domain size: %f\n", domain_size );
    fprintf(stderr, "   MLUps(Wall): %f\n\n", mlups );
    
    fprintf(fileout, "%d    %f    %f    %f\n", N_X, domain_size, sec_elapsed, mlups);
    fclose(fileout);
    
    return 0;
}


