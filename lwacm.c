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

#define T_MAX     10
#define N_X       10
#define N_Y       10
#define N_Z       10
#define ALPHA_MAX 18

const double omega = 1.5;

// time and index variables
int t = 0;
int x = 0;
int y = 0;
int z = 0;

// array to store value of p, t=0 is t, t=1 is t+1 
//       t | x  | y  |  z
double p[2][N_X][N_Y][N_Z];

// array to store value u
//     t| x | y | z |u_xyz, t=0 is t, t=1 is t+1
double u[2][N_X][N_Y][N_Z][3];

// array to store f, f(e), and f(e,o)
double f[19]     = { 0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0 };
double f_e[19]   = { 0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0 };
double f_e_o[19] = { 0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0 };

// temperary variable for calculation
double u_xi = 0;   // u(x-xi) * xi
double u_x_xi = 0;   // u(x) * xi
double u_2 = 0;    // u square
double f_x_t = 0;  // f(e, o)(x, t)

// for load p and u in alpha callback function
double p_load = 0;
double u_load[3];


// compute new rho()
double compute_p()
{
    return f[0]+f[1]+f[2]+f[3]+f[4]+f[5]+f[6]+f[7]+f[8]+f[9]+f[10]+f[11]+f[12]+f[13]+f[14]+f[15]+f[16]+f[17]+f[18];
}


// compute new u for coordinate x,y,z
void compute_u( int x, int y, int z, double p)
{    
    u[1][x][y][z][0] = ( f[1]-f[2]+f[7]-f[8]+f[9]-f[10]+f[11]-f[12]+f[13]-f[14] )/p;
    u[1][x][y][z][1] = ( f[3]-f[4]+f[7]+f[8]-f[9]-f[10]+f[15]-f[16]+f[17]-f[18] )/p;
    u[1][x][y][z][2] = ( f[5]-f[6]+f[11]+f[12]-f[13]-f[14]+f[15]+f[16]-f[17]-f[18] )/p;
}


// call back functions to calculate f for each alpha from 0 to 18
void alpha_0_call()
{      
      // alpha = 0 calculate u * xi and u square  xi(0)  {  0,  0,  0, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x][y][z];
             
      u_load[0] = u[0][x][y][z][0];
      u_load[1] = u[0][x][y][z][1];
      u_load[2] = u[0][x][y][z][2];
             
      u_xi = 0;
      u_2  = u_load[0] * u_load[0] + u_load[0] * u_load[0] + u_load[0] * u_load[0];
      u_x_xi = 0;
      
      //  the follwoing code can be abbreviated
      // step 5, compute f(e)(x-xi[0], t)   
      f_e[0] = 1.0/3.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[0] = 3.0 * 1.0/3.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/3.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9
      f[0] =  f_e[0] + 2*( omega-1/omega )*( f_x_t - f_e_o[0] ) ;  // Eq.11
}

void alpha_1_call()
{
      // alpha = 1 calculate u * xi and u square  xi(1)  {  1,  0,  0, }
      
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x-1][y][z];
      
      u_load[0] = u[0][x-1][y][z][0];
      u_load[1] = u[0][x-1][y][z][1];
      u_load[2] = u[0][x-1][y][z][2];
      
      u_xi = u_load[0];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][0];
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e[1] = 1.0/18.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[1] = 3.0 * 1.0/18.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/18.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[1] =  f_e[1] + 2*( omega-1/omega )*( f_x_t - f_e_o[1] ) ;  // Eq.11
}

void alpha_2_call()
{
      // alpha = 0 calculate u * xi and u square  xi(2)  { -1,  0,  0, }
      
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x-1][y][z];
      
      u_load[0] = u[0][x-1][y][z][0];
      u_load[1] = u[0][x-1][y][z][1];
      u_load[2] = u[0][x-1][y][z][2];
      
      u_xi = u_load[0]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][0] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e[2] = 1.0/18.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[2] = 3.0 * 1.0/18.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/18.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[2] =  f_e[2] + 2*( omega-1/omega )*( f_x_t - f_e_o[2] ) ;  // Eq.11
}

void alpha_3_call()
{
      // alpha = 0 calculate u * xi and u square  xi(3)  {  0,  1,  0, }
      
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x][y-1][z];
      
      u_load[0] = u[0][x][y-1][z][0];
      u_load[1] = u[0][x][y-1][z][1];
      u_load[2] = u[0][x][y-1][z][2];
      
      u_xi = u_load[1];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][1];
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e[3] = 1.0/18.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[3] = 3.0 * 1.0/18.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/18.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[3] =  f_e[3] + 2*( omega-1/omega )*( f_x_t - f_e_o[3] ) ;  // Eq.11
}

void alpha_4_call()
{
      // alpha = 0 calculate u * xi and u square  xi(4)  {  0, -1,  0, }
      
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x][y+1][z];
      
      u_load[0] = u[0][x][y+1][z][0];
      u_load[1] = u[0][x][y+1][z][1];
      u_load[2] = u[0][x][y+1][z][2];
      
      u_xi = u_load[1]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][1] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e[4] = 1.0/18.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[4] = 3.0 * 1.0/18.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/18.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[4] =  f_e[4] + 2*( omega-1/omega )*( f_x_t - f_e_o[4] ) ;  // Eq.11
}

void alpha_5_call()
{
      // alpha = 0 calculate u * xi and u square  xi(5)  {  0,  0,  1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x][y][z-1];
      
      u_load[0] = u[0][x][y][z-1][0];
      u_load[1] = u[0][x][y][z-1][1];
      u_load[2] = u[0][x][y][z-1][2];
      
      u_xi = u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][2];
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e[5] = 1.0/18.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[5] = 3.0 * 1.0/18.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/18.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[5] =  f_e[5] + 2*( omega-1/omega )*( f_x_t - f_e_o[5] ) ;  // Eq.11
}

void alpha_6_call()
{
      // alpha = 0 calculate u * xi and u square  xi(6)  {  0,  0, -1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x][y][z+1];
      
      u_load[0] = u[0][x][y][z+1][0];
      u_load[1] = u[0][x][y][z+1][1];
      u_load[2] = u[0][x][y][z+1][2];
      
      u_xi = u_load[2]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][2] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e[6] = 1.0/18.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[6] = 3.0 * 1.0/18.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/18.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[6] =  f_e[6] + 2*( omega-1/omega )*( f_x_t - f_e_o[6] ) ;  // Eq.11
}

void alpha_7_call()
{
      // alpha = 0 calculate u * xi and u square  xi(7)  {  1,  1,  0, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x-1][y-1][z];
      
      u_load[0] = u[0][x-1][y-1][z][0];
      u_load[1] = u[0][x-1][y-1][z][1];
      u_load[2] = u[0][x-1][y-1][z][2];
      
      u_xi = u_load[0] + u_load[1];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][0] + u[0][x][y][z][1];
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e[7] = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[7] = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[7] =  f_e[7] + 2*( omega-1/omega )*( f_x_t - f_e_o[7] ) ;  // Eq.11
}

void alpha_8_call()
{
      // alpha = 0 calculate u * xi and u square  xi(8)  { -1,  1,  0, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x+1][y-1][z];
      
      u_load[0] = u[0][x+1][y-1][z][0];
      u_load[1] = u[0][x+1][y-1][z][1];
      u_load[2] = u[0][x+1][y-1][z][2];
      
      u_xi = u_load[0]*(-1) + u_load[1];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][0] * (-1) + u[0][x][y][z][1];
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e[8] = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[8] = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[8] =  f_e[8] + 2*( omega-1/omega )*( f_x_t - f_e_o[8] ) ;  // Eq.11
}

void alpha_9_call()
{
      // alpha = 0 calculate u * xi and u square  xi(9)  {  1, -1,  0, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x-1][y+1][z];
      
      u_load[0] = u[0][x-1][y+1][z][0];
      u_load[1] = u[0][x-1][y+1][z][1];
      u_load[2] = u[0][x-1][y+1][z][2];
      
      u_xi = u_load[0] + u_load[1]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][0] + u[0][x][y][z][1] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e[9] = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[9] = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[9] =  f_e[9] + 2*( omega-1/omega )*( f_x_t - f_e_o[9] ) ;  // Eq.11
}

void alpha_10_call()
{
      // alpha = 0 calculate u * xi and u square  xi(10) { -1, -1,  0, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x+1][y+1][z];
      
      u_load[0] = u[0][x+1][y+1][z][0];
      u_load[1] = u[0][x+1][y+1][z][1];
      u_load[2] = u[0][x+1][y+1][z][2];
      
      u_xi = u_load[0]*(-1) + u_load[1]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][0] * (-1) + u[0][x][y][z][1] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e[10] = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[10] = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[10] =  f_e[10] + 2*( omega-1/omega )*( f_x_t - f_e_o[10] ) ;  // Eq.11
}

void alpha_11_call()
{
      // alpha = 0 calculate u * xi and u square  xi(11) {  1,  0,  1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x-1][y][z-1];
      
      u_load[0] = u[0][x-1][y][z-1][0];
      u_load[1] = u[0][x-1][y][z-1][1];
      u_load[2] = u[0][x-1][y][z-1][2];
      
      u_xi = u_load[0] + u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][0] + u[0][x][y][z][2];
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e[11] = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[11] = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[11] =  f_e[11] + 2*( omega-1/omega )*( f_x_t - f_e_o[11] ) ;  // Eq.11
}

void alpha_12_call()
{
      // alpha = 0 calculate u * xi and u square  xi(12) { -1,  0,  1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x+1][y][z-1];
      
      u_load[0] = u[0][x+1][y][z-1][0];
      u_load[1] = u[0][x+1][y][z-1][1];
      u_load[2] = u[0][x+1][y][z-1][2];
      
      u_xi = u_load[0]*(-1) + u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][0] * (-1) + u[0][x][y][z][2];
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e[12] = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[12] = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[12] =  f_e[2] + 2*( omega-1/omega )*( f_x_t - f_e_o[12] ) ;  // Eq.11
}

void alpha_13_call()
{
      // alpha = 0 calculate u * xi and u square  xi(13) {  1,  0, -1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x-1][y][z+1];
      
      u_load[0] = u[0][x-1][y][z+1][0];
      u_load[1] = u[0][x-1][y][z+1][1];
      u_load[2] = u[0][x-1][y][z+1][2];
      
      u_xi = u_load[0] + u_load[2]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][0] + u[0][x][y][z][2] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e[13] = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[13] = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[13] =  f_e[13] + 2*( omega-1/omega )*( f_x_t - f_e_o[13] ) ;  // Eq.11
}

void alpha_14_call()
{
      // alpha = 0 calculate u * xi and u square  xi(14) { -1,  0, -1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x+1][y][z+1];
      
      u_load[0] = u[0][x+1][y][z+1][0];
      u_load[1] = u[0][x+1][y][z+1][1];
      u_load[2] = u[0][x+1][y][z+1][2];
      
      u_xi = u_load[0]*(-1) + u_load[2]*(-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][0] * (-1) + u[0][x][y][z][2] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e[14] = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[14] = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[14] =  f_e[14] + 2*( omega-1/omega )*( f_x_t - f_e_o[14] ) ;  // Eq.11
}

void alpha_15_call()
{
      // alpha = 0 calculate u * xi and u square  xi(15) {  0,  1,  1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x][y-1][z-1];
      
      u_load[0] = u[0][x][y-1][z-1][0];
      u_load[1] = u[0][x][y-1][z-1][1];
      u_load[2] = u[0][x][y-1][z-1][2];
      
      u_xi = u_load[1] + u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][1] + u[0][x][y][z][2];
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e[15] = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[15] = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[15] =  f_e[15] + 2*( omega-1/omega )*( f_x_t - f_e_o[15] ) ;  // Eq.11
}

void alpha_16_call()
{
      // alpha = 0 calculate u * xi and u square  xi(16) {  0, -1,  1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x][y+1][z-1];
      
      u_load[0] = u[0][x][y+1][z-1][0];
      u_load[1] = u[0][x][y+1][z-1][1];
      u_load[2] = u[0][x][y+1][z-1][2];
      
      u_xi = u_load[1]*(-1) + u_load[2];
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][1] * (-1) + u[0][x][y][z][2];
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e[16] = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[16] = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[16] =  f_e[16] + 2*( omega-1/omega )*( f_x_t - f_e_o[16] ) ;  // Eq.11
}

void alpha_17_call()
{
      // alpha = 0 calculate u * xi and u square  xi(17) {  0,  1, -1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x][y-1][z+1];
      
      u_load[0] = u[0][x][y-1][z+1][0];
      u_load[1] = u[0][x][y-1][z+1][1];
      u_load[2] = u[0][x][y-1][z+1][2];
      
      u_xi = u_load[1] + u_load[2] * (-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][1] + u[0][x][y][z][2] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)
      f_e[17] = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[17] = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[17] =  f_e[17] + 2*( omega-1/omega )*( f_x_t - f_e_o[17] ) ;  // Eq.11
}

void alpha_18_call()
{
      // alpha = 0 calculate u * xi and u square  xi(18) {  0, -1, -1, }
            
      //load p(x-xia) and u(x-xia)
      p_load = p[0][x][y+1][z+1];
      
      u_load[0] = u[0][x][y+1][z+1][0];
      u_load[1] = u[0][x][y+1][z+1][1];
      u_load[2] = u[0][x][y+1][z+1][2];
      
      u_xi = u_load[1] * (-1) + u_load[2] * (-1);
      u_2  = u_load[0] * u_load[0] + u_load[1] * u_load[1] + u_load[2] * u_load[2];
      u_x_xi = u[0][x][y][z][1] * (-1) + u[0][x][y][z][2] * (-1);
      
      // step 5, compute f(e)(x-xi[0], t)   
      f_e[18] = 1.0/36.0 * p_load * ( 1 + 3*u_xi + 4.5*u_xi*u_xi - 1.5*u_2 );  // Eq.8
      
      // step 6, compute f(e, o)(x-xi[0], t)
      f_e_o[18] = 3.0 * 1.0/36.0 * p_load * u_xi;  // Eq.10
      
      // step 8.5, compute f(e, o)(x, t)
      f_x_t = 3.0 * 1.0/36.0 * p[0][x][y][z] * u_x_xi;
      
      // step 9, compute f(x, t+1)
      f[18] =  f_e[18] + 2*( omega-1/omega )*( f_x_t - f_e_o[18] ) ;  // Eq.11
}


int main()
{
    printf("\nlwacm start...\n");
    
    int i = 0;
    

    // initialize p and u
    for( x = 0; x < N_X; x++)
    {
      for( y = 0; y < N_Y; y++)
      {
        for( z = 0; z < N_Z; z++)
        {
            p[0][x][y][z] = 1.0f; // set it to 1.0
            
            u[0][x][y][z][0] = 0;
            u[0][x][y][z][1] = 0;
            u[0][x][y][z][2] = 0;
        }
      }
    }
    
    // for all time step t
    for( t = 0; t < T_MAX; t++)
    {
        
        printf("time step = %d \n", t);
        
        // for all mesh point x
        for( x = 0; x < N_X; x++)
        {
          for( y = 0; y < N_Y; y++)
          {
            for( z = 0; z < N_Z; z++)
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
                p[1][x][y][z] = compute_p();
      
                // step 12, compute u(x, t+1)
                compute_u(x, y, z, p[1][x][y][z]);
      
                // step 13, store p(x, t+1)
                p[0][x][y][z] = p[1][x][y][z];
      
                // step 13, store u(x, t+1)
                u[0][x][y][z][0] = u[1][x][y][z][0];
                u[0][x][y][z][1] = u[1][x][y][z][1];
                u[0][x][y][z][2] = u[1][x][y][z][2];
            }
          }
        }
        
        for( i = 0; i < 19; i++)
        {
            printf("f(%d) = %e \n", i, f[i] );
        }
        
    }
    
    return 0;
}


