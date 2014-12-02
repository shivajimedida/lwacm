/*
 *  Link-wise Artificial Compressibility Method
 *  by Yifan Yang with supervisor Thomas Zeiser @ FAU
 *                                         01.12.2014
 */

#include <stdio.h>

#define T_MAX     1000
#define N_X       10
#define N_Y       10
#define N_Z       10
#define ALPHA_MAX 18

const double omega = 1.5;

// array to store value of p, t=0 is t, t=1 is t+1 
//            t | x  | y  |  z
double p_store[2][N_X][N_Y][N_Z];

// array to store value u
//     t| x | y | z |u_xyz, t=0 is t, t=1 is t+1
double u[2][N_X][N_Y][N_Z][3];

// array to store f, f(e), and f(e,o)
double f[19]     = { 0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0 };
double f_e[19]   = { 0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0 };
double f_e_o[19] = { 0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0 };

// u square
double u_square( int x, int y, int z )
{
    double result = u[0][x][y][z][0] * u[0][x][y][z][0];
    result += u[0][x][y][z][1] * u[0][x][y][z][1];
    result += u[0][x][y][z][2] * u[0][x][y][z][2];
    
    return  result;
}

// multiply funtion for u and xi stencil
double u_x_xi( int x, int y, int z, int alpha )
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
    
    return u[0][x][y][z][0]*xix + u[0][x][y][z][1]*xiy + u[0][x][y][z][2]*xiz;
}

// compute new u for coordinate x,y,z
void compute_u( int x, int y, int z, double p)
{    
    u[1][x][y][z][0] = ( f[1]-f[2]+f[7]-f[8]+f[9]-f[10]+f[11]-f[12]+f[13]-f[14] )/p;
    u[1][x][y][z][1] = ( f[3]-f[4]+f[7]+f[8]-f[9]-f[10]+f[15]-f[16]+f[17]-f[18] )/p;
    u[1][x][y][z][2] = ( f[5]-f[6]+f[11]+f[12]-f[13]-f[14]+f[15]+f[16]-f[17]-f[18] )/p;
}

// compute new rho()
double compute_p()
{    
    int i = 0;
    double result = 0;
    
    // sum up all the f(a)
    for( i = 0; i < 19; i++)
    {
        result += f[i];
    }
    
    return result;
}

// put the data in t+1 into t
void store_p(int x, int y, int z )
{
    p_store[0][x][y][z] = p_store[1][x][y][z];
}

// compute new u for coordinate x,y,z
void store_u( int x, int y, int z )
{    
    u[0][x][y][z][0] = u[1][x][y][z][0];
    u[0][x][y][z][1] = u[1][x][y][z][1];
    u[0][x][y][z][2] = u[1][x][y][z][2];
}


// function of w(alpha)
double w(int alpha)
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
    
    int i = 0;
    
    // variables
    int alpha = 0;
    int t = 0;

    int x = 0;
    int y = 0;
    int z = 0;
    
    // initialize p and u
    for( x = 0; x < N_X; x++)
    {
      for( y = 0; y < N_Y; y++)
      {
        for( z = 0; z < N_Z; z++)
        {
            p_store[0][x][y][z] = 1.0f; // set it to 1.0
            
            u[0][x][y][z][0] = 0;
            u[0][x][y][z][1] = 0;
            u[0][x][y][z][2] = 0;
        }
      }
    }
    
    // temperary variable for calculation
    double u_temp = 0;
    double u_2_temp = 0;
    double f_x_t = 0;
    
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
                
                // for all index alpha
                for( alpha = 0; alpha <= ALPHA_MAX; alpha++)
                {
                    // load p and u
                    
                    u_temp = u_x_xi( x, y, z, alpha );    // u times xi(alpha)
                    //printf("u_temp = %f \n", u_temp);
                    
                    u_2_temp = u_square( x, y, z );       // u square
                    
                    // step 5, compute f(e)(x-xi, t)
                    f_e[alpha]   =  w(alpha) * p_store[0][x][y][z] * ( 1 + 3*u_temp + 9*u_temp*u_temp/2 - 3*u_2_temp/2 );  // Eq.8
                    
                    // step 6, compute f(e, o)(x-xi, t)
                    f_e_o[alpha] =  w(alpha) * p_store[0][x][y][z] * ( 1 + 9*u_temp*u_temp/2 - 3*u_2_temp/2 );  // Eq.10
                    
                }
                
                // for all index alpha
                for( alpha = 0; alpha <= ALPHA_MAX; alpha++)
                {
                    // step 9, compute f(a)(x, t+1)
                    u_2_temp = u_square( x, y, z );       // u square
                    
                    f_x_t = w(alpha) * p_store[0][x][y][z] * ( 1 - 3*u_2_temp/2 );
                    
                    f[alpha] =  f_e[alpha] + 2*( omega-1/omega )*( f_x_t - f_e_o[alpha] ) ;  // Eq.11
                }
                
                // step 11, compute p(x, t+1)
                p_store[1][x][y][z] = compute_p();
                
                // step 12, compute u(x, t+1)
                compute_u(x, y, z, p_store[1][x][y][z]);
                
                // step 13, store p(x, t+1) and u(x, t+1)
                
                store_p(x, y, z);
                store_u(x, y, z);
                
                //printf("x = %d, y = %d, z = %d \n", x, y, z );
                
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


