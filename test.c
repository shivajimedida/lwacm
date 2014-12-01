#include <stdio.h>

#define T_MAX     1000
#define N_X       10
#define N_Y       10
#define N_Z       10
#define ALPHA_MAX 18

int main()
{
    printf("\nprogram start\n");
    
    // variables
    int alpha = 0;
    int t = 0;
    
    
    int x = 0;
    int y = 0;
    int z = 0;
    
    
    // array to store value u
    //    | x | y | z | uxyz
    float u[10][10][10][3];
    
    
        // for all mesh point x
        for( x = 1; x <= N_X; x++)
        {
          for( y = 1; y <= N_Y; y++)
          {
            for( z = 1; z <= N_Z; z++)
            {
                
                u[x-1][y-1][z-1][0] = 0.1;
                u[x-1][y-1][z-1][1] = 0.1;
                u[x-1][y-1][z-1][2] = 0.1;
                
                
                printf("add = %f \n", u[x-1][y-1][z-1][0]+u[x-1][y-1][z-1][1]+u[x-1][y-1][z-1][2]);
                 
            }
          }
        }
    
    return 0;
}


