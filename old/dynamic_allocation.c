#include <stdio.h>
#include <stdlib.h>


int main()
{
    int i = 0;
    
    // define the 4 dimensional array  p[2][N_X][N_Y][N_Z];
    
    // domain boundary
    int N_X = 40;
    int N_Y = 40;
    int N_Z = 40;
    
    // iteration index
    int x = 0;
    int y = 0;
    int z = 0;
    
    long bytes;
    
    
    printf("\n>> start allocating mem for array of size %d\n", N_X);
    
    
    // 4 dimensional array pointer 
    double ****data;
    
    // allocate mem for array pointer
    data = malloc( 2 * sizeof(double ***) );
    
    if(data == NULL)
    {
        fprintf(stderr, "\n>> error allocating memory for array pointer(size = 2)\n\n");
        return 1;
    }
    
    // allocatw mem for t_now and t_next
    for(i = 0; i < 2; i++)
    {
        data[i] = malloc( (N_X+2) * sizeof(double**) );
        
        if(data[i] == NULL)
        {
            fprintf(stderr, "\n>> error allocating memory for t_now/t_next(size = %d)\n\n", N_X);
            return 1;
        }
        //else { printf("\n>> mem [%d] ok\n", i); }
        
        // allocatw mem for N_X
        for(x = 0; x < N_X+2; x++)
        {
            data[i][x] = malloc( (N_Y+2) * sizeof(double*) );
        
            if(data[i][x] == NULL)
            {
                fprintf(stderr, "\n>> error allocating memory for N_X(size = %d)\n\n", N_Y);
                return 1;
            }
            //else { printf("\n>> mem [%d][%d] ok\n", i, x); }
            
            // allocatw mem for N_X(actuall data)
            for(y = 0; y < N_Y+2; y++)
            {
                data[i][x][y] = malloc( (N_Z+2) * sizeof(double) );
        
                if(data[i][x][y] == NULL)
                {
                    fprintf(stderr, "\n>> error allocating memory for N_Y(size =%d)\n\n", N_Z);
                    return 1;
                }
                //else { printf("\n>> mem [%d][%d][%d] ok\n", i, x, y); }
            }
        }
    
    }
    
    
    
    bytes = 2 * (N_X+2) * (N_Y+2) * (N_Z+2) * sizeof(double);
    

    if( bytes > 1024*1024 )
    {
        printf("\n>> successfully allocated %ld MBytes memory for array data\n", bytes/1024/1024);
    }
    
    else if( bytes > 1024 )
    {
        printf("\n>> successfully allocated %ld KBytes memory for array data\n", bytes/1024);
    }
    
    else
    {
        printf("\n>> successfully allocated %ld Bytes memory for array data\n", bytes);
    }
    
    
    
    printf("\n>> initialize array data\n");
    
    // initialize p and u
    for( x = 0; x < N_X+2; x++)
    {
      for( y = 0; y < N_Y+2; y++)
      {
        for( z = 0; z < N_Z+2; z++)
        {
            data[0][x][y][z] = 0;
            data[1][x][y][z] = 1;
        }
      }
    }
    
    
    
    
    // free mem for z
    for(i = 0; i < 2; i++)
    {
        for(x = 0; x < N_X+2; x++)
        {
            for(y = 0; y < N_Y+2; y++)
            {
                free( data[i][x][y] );
            }  
        }
    }
    
    printf("\n>> z freed\n");
    
    // free mem for y
    for(i = 0; i < 2; i++)
    {
        for(x = 0; x < N_X+2; x++)
        {
            free( data[i][x] );
        }
    }
    
    printf("\n>> y freed\n");
    
    // free mem for x
    for(i = 0; i < 2; i++)
    {
        free( data[i] );
    }
    
    printf("\n>> x freed\n");
    
    // free data pointer
    free(data);
    
    
    printf("\n>> all memory freed\n");
    
    return 0;
}





