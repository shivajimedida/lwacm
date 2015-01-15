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
#include <unistd.h>

#define T_MAX     20
#define N_X       10
#define N_Y       10
#define N_Z       10


int x = 0;
int y = 0;
int z = 0;

int xx = 0;
int yy = 0;
int zz = 0;

int i = 0;
int t_next = 0;


// array to store value of p, t=0 is t, t=1 is t+1
// use x+1, y+1, z+1 to avoid illegal access like p(x-xi_alpha), same for u[]
// 
//       t |  x  |  y   |  z
double p[2][N_X+2][N_Y+2][N_Z+2];

// array to store value u
//       t |  x  |  y   |  z    |u_xyz, t=0 is t, t=1 is t+1
double u[2][N_X+2][N_Y+2][N_Z+2][3];

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


int main()
{
    printf("\n>> program start...\n");
    printf("\n   domain size : N_X = %d, N_Y = %d, N_Z = %d\n", N_X, N_Y, N_Z );
    
    
    
    // initialize p
    for( x = 0; x < N_X+2; x++)
    {
        for( y = 0; y < N_Y+2; y++)
	{
	    for( z = 0; z < N_Z+2; z++)
	    {
	        p[0][x][y][z] = 0; // set all to 0
	    }
	}
    }
    
    for( x = 1; x < N_X+1; x++)
    {
        for( y = 1; y < N_Y+1; y++)
	{
	    for( z = 1; z < N_Z+1; z++)
	    {
	        p[0][x][y][z] = 1.0; // set inner to 1
	    }
	}
    }
    
    
    printf("\n______________________________________________________" );
    printf("\n>> before setting bc " );
    
    // print p
    for( x = 0; x < N_X+2; x++)
    {
        for( y = 0; y < N_Y+2; y++)
	{
	    for( z = 0; z < N_Z+2; z++)
	    {
	        
                if ( x == 0 ) { xx=N_X; flag = true; }
                if ( y == 0 ) { yy=N_Y; flag = true; }
                if ( z == 0 ) { zz=N_Z; flag = true; }

                if ( x == N_X+1 ) { xx=1; flag = true; }
                if ( y == N_Y+1 ) { yy=1; flag = true; }
                if ( z == N_Z+1 ) { zz=1; flag = true; }
	        
	        printf("\n   coord = ( %d, %d, %d ),  p = %e", x, y, z, p[0][x][y][z] );
	    }
	}
    }
    
    
    printf("\n______________________________________________________" );
    printf("\n>> after setting bc " );
        

        //#########################################################################################################################
        //##################################################### start #############################################################
        // explicit peridic boundary conditions
        // surfaces: xy
        for( x = 1; x < N_X+1; x++)
        {
  	  for( y = 1; y < N_Y+1; y++)
  	  {
  	    p[t_next][ x ][ y ][ 0 ] = p[t_next][ x ][ y ][N_Z];
  	    u[t_next][ x ][ y ][ 0 ][0] = u[t_next][ x ][ y ][N_Z][0];
  	    u[t_next][ x ][ y ][ 0 ][1] = u[t_next][ x ][ y ][N_Z][1];
  	    u[t_next][ x ][ y ][ 0 ][2] = u[t_next][ x ][ y ][N_Z][2];
  	    
  	    p[t_next][ x ][ y ][N_Z+1] = p[t_next][ x ][ y ][1];
  	    u[t_next][ x ][ y ][N_Z+1][0] = u[t_next][ x ][ y ][1][0];
  	    u[t_next][ x ][ y ][N_Z+1][1] = u[t_next][ x ][ y ][1][1];
  	    u[t_next][ x ][ y ][N_Z+1][2] = u[t_next][ x ][ y ][1][2];
  	    
  	  }
  	}
        
        // surfaces: xz			      
        for( x = 1; x < N_X+1; x++)
        {
  	  for( z = 1; z < N_Z+1; z++)
  	  {
  	    p[t_next][ x ][ 0 ][ z ] = p[t_next][ x ][N_Y][ z ];
  	    u[t_next][ x ][ 0 ][ z ][0] = u[t_next][ x ][N_Y][ z ][0];
  	    u[t_next][ x ][ 0 ][ z ][1] = u[t_next][ x ][N_Y][ z ][1];
  	    u[t_next][ x ][ 0 ][ z ][2] = u[t_next][ x ][N_Y][ z ][2];
  	    
  	    p[t_next][ x ][N_Y+1][ z ] = p[t_next][ x ][ 1 ][ z ];
  	    u[t_next][ x ][N_Y+1][ z ][0] = u[t_next][ x ][ 1 ][ z ][0];
  	    u[t_next][ x ][N_Y+1][ z ][1] = u[t_next][ x ][ 1 ][ z ][1];
  	    u[t_next][ x ][N_Y+1][ z ][2] = u[t_next][ x ][ 1 ][ z ][2];
  	    
  	  }
  	}
        
        // surfaces: yz
        for( y = 1; y < N_Y+1; y++)
  	{
  	  for( z = 1; z < N_Z+1; z++)
  	  {
  	    p[t_next][ 0 ][ y ][ z ] = p[t_next][N_X][ y ][ z ];
  	    u[t_next][ 0 ][ y ][ z ][0] = u[t_next][N_X][ y ][ z ][0];
  	    u[t_next][ 0 ][ y ][ z ][1] = u[t_next][N_X][ y ][ z ][1];
  	    u[t_next][ 0 ][ y ][ z ][2] = u[t_next][N_X][ y ][ z ][2];
  	    
  	    p[t_next][N_X+1][ y ][ z ] = p[t_next][ 1 ][ y ][ z ];
  	    u[t_next][N_X+1][ y ][ z ][0] = u[t_next][ 1 ][ y ][ z ][0];
  	    u[t_next][N_X+1][ y ][ z ][1] = u[t_next][ 1 ][ y ][ z ][1];
  	    u[t_next][N_X+1][ y ][ z ][2] = u[t_next][ 1 ][ y ][ z ][2];
  	    
  	  }
  	}
  	
  	
  	// lines x: 4
        for( x = 1; x < N_X+1; x++)
  	{
  	    // (x, 0, 0) 
  	    p[t_next][ x ][ 0 ][ 0 ] = p[t_next][ x ][N_Y][N_Z];
  	    u[t_next][ x ][ 0 ][ 0 ][0] = u[t_next][ x ][N_Y][N_Z][0];
  	    u[t_next][ x ][ 0 ][ 0 ][1] = u[t_next][ x ][N_Y][N_Z][1];
  	    u[t_next][ x ][ 0 ][ 0 ][2] = u[t_next][ x ][N_Y][N_Z][2];
  	    
  	    // (x, N_Y+1, N_Z+1)
  	    p[t_next][ x ][N_Y+1][N_Z+1] = p[t_next][ x ][ 1 ][ 1 ];
  	    u[t_next][ x ][N_Y+1][N_Z+1][0] = u[t_next][ x ][ 1 ][ 1 ][0];
  	    u[t_next][ x ][N_Y+1][N_Z+1][1] = u[t_next][ x ][ 1 ][ 1 ][1];
  	    u[t_next][ x ][N_Y+1][N_Z+1][2] = u[t_next][ x ][ 1 ][ 1 ][2];
  	    
  	    // (x, 0, N_Z+1)
  	    p[t_next][ x ][ 0 ][N_Z+1] = p[t_next][ x ][N_Y][ 1 ];
  	    u[t_next][ x ][ 0 ][N_Z+1][0] = u[t_next][ x ][N_Y][ 1 ][0];
  	    u[t_next][ x ][ 0 ][N_Z+1][1] = u[t_next][ x ][N_Y][ 1 ][1];
  	    u[t_next][ x ][ 0 ][N_Z+1][2] = u[t_next][ x ][N_Y][ 1 ][2];
  	    
  	    // (x, N_Y+1, 0)
  	    p[t_next][ x ][N_Y+1][ 0 ] = p[t_next][ x ][ 1 ][N_Z];
  	    u[t_next][ x ][N_Y+1][ 0 ][0] = u[t_next][ x ][ 1 ][N_Z][0];
  	    u[t_next][ x ][N_Y+1][ 0 ][1] = u[t_next][ x ][ 1 ][N_Z][1];
  	    u[t_next][ x ][N_Y+1][ 0 ][2] = u[t_next][ x ][ 1 ][N_Z][2];
  	}
        
  	// lines y: 4
        for( y = 1; y < N_Y+1; y++)
  	{
  	    // (0, y, 0) 
  	    p[t_next][ 0 ][ y ][ 0 ] = p[t_next][N_X][y][N_Z];
  	    u[t_next][ 0 ][ y ][ 0 ][0] = u[t_next][N_X][ y ][N_Z][0];
  	    u[t_next][ 0 ][ y ][ 0 ][1] = u[t_next][N_X][ y ][N_Z][1];
  	    u[t_next][ 0 ][ y ][ 0 ][2] = u[t_next][N_X][ y ][N_Z][2];
  	    
  	    // (N_X+1, y, N_Z+1)
  	    p[t_next][N_X+1][ y ][N_Z+1] = p[t_next][ 1 ][ y ][  1  ];
  	    u[t_next][N_X+1][ y ][N_Z+1][0] = u[t_next][ 1 ][ y ][  1  ][0];
  	    u[t_next][N_X+1][ y ][N_Z+1][1] = u[t_next][ 1 ][ y ][  1  ][1];
  	    u[t_next][N_X+1][ y ][N_Z+1][2] = u[t_next][ 1 ][ y ][  1  ][2];
  	    
  	    // (0, y, N_Z+1)
  	    p[t_next][ 0 ][ y ][N_Z+1] = p[t_next][N_X][ y ][ 1 ];
  	    u[t_next][ 0 ][ y ][N_Z+1][0] = u[t_next][N_X][ y ][ 1 ][0];
  	    u[t_next][ 0 ][ y ][N_Z+1][1] = u[t_next][N_X][ y ][ 1 ][1];
  	    u[t_next][ 0 ][ y ][N_Z+1][2] = u[t_next][N_X][ y ][ 1 ][2];
  	    
  	    // (N_X+1, y, 0)
  	    p[t_next][N_X+1][ y ][0] = p[t_next][ 1 ][ y ][N_Z];
  	    u[t_next][N_X+1][ y ][0][0] = u[t_next][ 1 ][ y ][N_Z][0];
  	    u[t_next][N_X+1][ y ][0][1] = u[t_next][ 1 ][ y ][N_Z][1];
  	    u[t_next][N_X+1][ y ][0][2] = u[t_next][ 1 ][ y ][N_Z][2];
  	    
  	}
        
  	// lines z: 4
        for( z = 1; z < N_Z+1; z++)
  	{
  	    // (0, 0, z) 
  	    p[t_next][ 0 ][ 0 ][ z ] = p[t_next][N_X][N_Y][ z ];
  	    u[t_next][ 0 ][ 0 ][ z ][0] = u[t_next][N_X][N_Y][ z ][0];
  	    u[t_next][ 0 ][ 0 ][ z ][1] = u[t_next][N_X][N_Y][ z ][1];
  	    u[t_next][ 0 ][ 0 ][ z ][2] = u[t_next][N_X][N_Y][ z ][2];
  	    
  	    // (N_X+1, N_Y+1, z)
  	    p[t_next][N_X+1][N_Y+1][ z ] = p[t_next][ 1 ][  1  ][ z ];
  	    u[t_next][N_X+1][N_Y+1][ z ][0] = u[t_next][ 1 ][  1  ][ z ][0];
  	    u[t_next][N_X+1][N_Y+1][ z ][1] = u[t_next][ 1 ][  1  ][ z ][1];
  	    u[t_next][N_X+1][N_Y+1][ z ][2] = u[t_next][ 1 ][  1  ][ z ][2];
  	    
  	    // (0, N_Y+1, z)
  	    p[t_next][ 0 ][N_Y+1][ z ] = p[t_next][N_X][ 1 ][ z ];
  	    u[t_next][ 0 ][N_Y+1][ z ][0] = u[t_next][N_X][ 1 ][ z ][0];
  	    u[t_next][ 0 ][N_Y+1][ z ][1] = u[t_next][N_X][ 1 ][ z ][1];
  	    u[t_next][ 0 ][N_Y+1][ z ][2] = u[t_next][N_X][ 1 ][ z ][2];
  	    
  	    // (N_X+1, 0, z)
  	    p[t_next][N_X+1][ 0 ][ z ] = p[t_next][ 1 ][N_Y][ z ];
  	    u[t_next][N_X+1][ 0 ][ z ][0] = u[t_next][ 1 ][N_Y][ z ][0];
  	    u[t_next][N_X+1][ 0 ][ z ][1] = u[t_next][ 1 ][N_Y][ z ][1];
  	    u[t_next][N_X+1][ 0 ][ z ][2] = u[t_next][ 1 ][N_Y][ z ][2];
  	    
  	}
  	
  	
  	
  	// points: 8
  	/*
  	// (0, 0, 0) <- (N_X, N_Y, N_Z)
        p[t_next][ 0 ][ 0 ][ 0 ] = p[t_next][N_X][N_Y][N_Z];
        u[t_next][ 0 ][ 0 ][ 0 ][0] = u[t_next][N_X][N_Y][N_Z][0];
        u[t_next][ 0 ][ 0 ][ 0 ][1] = u[t_next][N_X][N_Y][N_Z][1];
        u[t_next][ 0 ][ 0 ][ 0 ][2] = u[t_next][N_X][N_Y][N_Z][2];
  	
  	// (0, 0, N_Z+1) <- (N_X, N_Y, 1)
        p[t_next][ 0 ][ 0 ][N_Z+1] = p[t_next][N_X][N_Y][ 1 ];
        u[t_next][ 0 ][ 0 ][N_Z+1][0] = u[t_next][N_X][N_Y][ 1 ][0];
        u[t_next][ 0 ][ 0 ][N_Z+1][1] = u[t_next][N_X][N_Y][ 1 ][1];
        u[t_next][ 0 ][ 0 ][N_Z+1][2] = u[t_next][N_X][N_Y][ 1 ][2];
  	
  	// (0, N_Y+1, 0) <- (N_X, 1, N_Z)
        p[t_next][ 0 ][N_Y+1][ 0 ] = p[t_next][N_X][ 1 ][N_Z];
        u[t_next][ 0 ][N_Y+1][ 0 ][0] = u[t_next][N_X][ 1 ][N_Z][0];
        u[t_next][ 0 ][N_Y+1][ 0 ][1] = u[t_next][N_X][ 1 ][N_Z][1];
        u[t_next][ 0 ][N_Y+1][ 0 ][2] = u[t_next][N_X][ 1 ][N_Z][2];
        
  	// (0, N_Y+1, N_Z+1) <- (N_X, 1, 1)
        p[t_next][ 0 ][N_Y+1][N_Z+1] = p[t_next][N_X][ 1 ][ 1 ];
        u[t_next][ 0 ][N_Y+1][N_Z+1][0] = u[t_next][N_X][ 1 ][ 1 ][0];
        u[t_next][ 0 ][N_Y+1][N_Z+1][1] = u[t_next][N_X][ 1 ][ 1 ][1];
        u[t_next][ 0 ][N_Y+1][N_Z+1][2] = u[t_next][N_X][ 1 ][ 1 ][2];
        
  	// (N_X+1, 0, 0) <- (1, N_Y, N_Z)
        p[t_next][N_X+1][ 0 ][ 0 ] = p[t_next][ 1 ][N_Y][N_Z];
        u[t_next][N_X+1][ 0 ][ 0 ][0] = u[t_next][ 1 ][N_Y][N_Z][0];
        u[t_next][N_X+1][ 0 ][ 0 ][1] = u[t_next][ 1 ][N_Y][N_Z][1];
        u[t_next][N_X+1][ 0 ][ 0 ][2] = u[t_next][ 1 ][N_Y][N_Z][2];
        
  	// (N_X+1, 0, N_Z+1) <- (1, N_Y, 1)
        p[t_next][N_X+1][ 0 ][N_Z+1] = p[t_next][ 1 ][N_Y][ 1 ];
        u[t_next][N_X+1][ 0 ][N_Z+1][0] = u[t_next][ 1 ][N_Y][ 1 ][0];
        u[t_next][N_X+1][ 0 ][N_Z+1][1] = u[t_next][ 1 ][N_Y][ 1 ][1];
        u[t_next][N_X+1][ 0 ][N_Z+1][2] = u[t_next][ 1 ][N_Y][ 1 ][2];
        
  	// (N_X+1, N_Y+1, 0) <- (1, 1, N_Z)
        p[t_next][N_X+1][N_Y+1][ 0 ] = p[t_next][ 1 ][ 1 ][N_Z];
        u[t_next][N_X+1][N_Y+1][ 0 ][0] = u[t_next][ 1 ][ 1 ][N_Z][0];
        u[t_next][N_X+1][N_Y+1][ 0 ][1] = u[t_next][ 1 ][ 1 ][N_Z][1];
        u[t_next][N_X+1][N_Y+1][ 0 ][2] = u[t_next][ 1 ][ 1 ][N_Z][2];
        
  	// (N_X+1, N_Y+1, N_Z+1) <- (1, 1, 1)
        p[t_next][N_X+1][N_Y+1][N_Z+1] = p[t_next][ 1 ][ 1 ][ 1 ];
        u[t_next][N_X+1][N_Y+1][N_Z+1][0] = u[t_next][ 1 ][ 1 ][ 1 ][0];
        u[t_next][N_X+1][N_Y+1][N_Z+1][1] = u[t_next][ 1 ][ 1 ][ 1 ][1];
        u[t_next][N_X+1][N_Y+1][N_Z+1][2] = u[t_next][ 1 ][ 1 ][ 1 ][2];
        
        */
        
        //######################################################## end ############################################################
        //#########################################################################################################################
        
    
    
      
    // print p
    for( x = 0; x < N_X+2; x++)
    {
        for( y = 0; y < N_Y+2; y++)
	{
	    for( z = 0; z < N_Z+2; z++)
	    {
	        printf("\n   coord = ( %d, %d, %d ),  p = %e", x, y, z, p[0][x][y][z] );
	    }
	}
    }
        
        
        
    printf("\n\n");
    
    
    
    return 0;
}




