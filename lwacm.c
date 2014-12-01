/*
 *  Link-wise Artificial Compressibility Method
 *  by Yifan Yang with supervisor Thomas Zeiser @ FAU
 *                                         01.12.2014
 */

#include <stdio.h>

#define T_MAX     10
#define N_X       10
#define N_Y       10
#define N_Z       10
#define ALPHA_MAX 18

// variables
int alpha = 0;
int t = 0;

int x = 0;
int y = 0;
int z = 0;


const float omega = 1.5;

// array to store value of p
float p_store[N_X][N_Y][N_Z];


// array to store value u
//    | x | y | z |u_xyz
float u[10][10][10][3];


// array to store f, f(e), and f(e,o)
float f[19]     = { 0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0 };
float f_e[19]   = { 0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0 };
float f_e_o[19] = { 0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0 };

// u square
float u_square( int x, int y, int z )
{
    float result = u[x-1][y-1][z-1][0]*u[x-1][y-1][z-1][0];
    result += u[x-1][y-1][z-1][1]*u[x-1][y-1][z-1][1];
    result += u[x-1][y-1][z-1][2]*u[x-1][y-1][z-1][2];
    
    return  result;
}

// multiply funtion for u and xi stencil
float u_x_xi( int x, int y, int z, int alpha )
{
    int xix = 0;
    int xiy = 0;
    int xiz = 0;
    
    if( 0 == alpha ) { xix = 0; xiy = 0; xiz = 0; }
    
    else if( 1 == alpha )  { xix =  1; xiy =  0; xiz =  0; }
    else if( 2 == alpha )  { xix = -1; xiy =  0; xiz =  0; }
    else if( 3 == alpha )  { xix =  0; xiy =  1; xiz =  0; }
    else if( 4 == alpha )  { xix =  0; xiy = -1; xiz =  0; }
    else if( 5 == alpha )  { xix =  0; xiy =  0; xiz =  1; }
    else if( 6 == alpha )  { xix =  0; xiy =  0; xiz = -1; }
    
    else if( 7 == alpha )  { xix =  1; xiy =  1; xiz = 0; }
    else if( 8 == alpha )  { xix = -1; xiy =  1; xiz = 0; }
    else if( 9 == alpha )  { xix =  1; xiy = -1; xiz = 0; }
    else if( 10 == alpha ) { xix = -1; xiy = -1; xiz = 0; }
    
    else if( 11 == alpha ) { xix =  1; xiy = 0; xiz =  1; }
    else if( 12 == alpha ) { xix = -1; xiy = 0; xiz =  1; }
    else if( 13 == alpha ) { xix =  1; xiy = 0; xiz = -1; }
    else if( 14 == alpha ) { xix = -1; xiy = 0; xiz = -1; }
    
    else if( 15 == alpha ) { xix = 0; xiy =  1; xiz =  1; }
    else if( 16 == alpha ) { xix = 0; xiy = -1; xiz =  1; }
    else if( 17 == alpha ) { xix = 0; xiy =  1; xiz = -1; }
    else if( 18 == alpha ) { xix = 0; xiy = -1; xiz = -1; }
    
    else
    {
        printf("illegal alpha !\n");
    }
    
    return u[x-1][y-1][z-1][0]*xix + u[x-1][y-1][z-1][1]*xiy + u[x-1][y-1][z-1][2]*xiz;
}

// function of rho()
float compute_p( int x, int y, int z, int time_step, float *store)
{    
    int alpha = 0;
    
    // for all index alpha
    for( alpha = 0; alpha <= ALPHA_MAX; alpha++)
    {
    }
    
    return 1.0;
}

// function of w(alpha)
float w(int alpha)
{
    if( 0 == alpha ) { return 1.0/3.0; }
    
    else if( (alpha >= 1)||(alpha <= 6) ) { return 1.0/18.0; }
    
    else if( (alpha >= 7)||(alpha <= 18) ) { return 1.0/36.0; }
    
    else
    {
        printf("invalid alpha !!\n");
        return 0;
    }
}

int main()
{
    printf("\nprogram start\n");
    
    // initial value of p
    p_store[0][0][0] = 1.0f;
    
    // initial value of u
    u[0][0][0][0] = 0;
    u[0][0][0][1] = 0;
    u[0][0][0][2] = 0;
    
    float u_temp = 0;
    float u_2_temp = 0;
    
    float f_x_t = 0;
    
    
    // for all time step t
    for( t = 0; t <= T_MAX; t++)
    {
        
        // for all mesh point x
        for( x = 1; x <= N_X; x++)
        {
          for( y = 1; y <= N_Y; y++)
          {
            for( z = 1; z <= N_Z; z++)
            {
                
                // for all index alpha
                for( alpha = 0; alpha <= ALPHA_MAX; alpha++)
                {
                    // load p and u
                    
                    u_temp = u_x_xi( x, y, z, alpha );    // u times xi(alpha)
                    u_2_temp = u_square( x, y, z );       // u square
                    
                    // step 5, compute f(e)(x-xi, t)
                    f_e[alpha]   =  w(alpha) * p_store[x][y][z] * ( 1 + 3*u_temp + 9*u_temp*u_temp/2 - 3*u_2_temp/2 );  // Eq.8
                    
                    // step 6, compute f(e, o)(x-xi, t)
                    f_e_o[alpha] =  w(alpha) * p_store[x][y][z] * ( 1 + 9*u_temp*u_temp/2 - 3*u_2_temp/2 );  // Eq.10
                    
                }
                
                
                // for all index alpha
                for( alpha = 0; alpha <= ALPHA_MAX; alpha++)
                {
                    // step 9, compute f(a)(x, t+1)
                    u_2_temp = u_square( x, y, z );       // u square
                    
                    f_x_t = w(alpha) * p_store[x][y][z] * ( 1 - 3*u_2_temp/2 );
                    
                    f[alpha] =  f_e[alpha] + 2*( omega-1/omega )*( f_x_t - f_e_o[alpha] ) ;  // Eq.11
                }
                
                
                // step 11, compute p(x, t+1)
                
                // step 12, compute u(x, t+1)
                
                // step 13, store p(x, t+1) and u(x, t+1)
                
                printf("x = %d, y = %d, z = %d \n", x, y, z  );
                 
            }
          }
        }
    }
    
    return 0;
}


