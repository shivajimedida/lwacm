#include <stdio.h>
#include <stdlib.h>

int main()
{
    
    long i = 0;
    
    int* array;
    
    for( i = 0; i < 1000000000; i+= 1000000)
    {
        array = malloc(sizeof(double) * i);
        
        if (array == NULL)
        {
            fprintf( stderr, "\nCould not allocate memory for double at size : %ld !!\n", i);
            return 1;
        }
        
        else
        {
            fprintf( stderr, "\nSuccessfully allocated memory for double at size : %ld\n", i);
        }
        
        free(array);
    }
    
    //printf("\nsize of double = %d\n\n", sizeof(double));
    
    
    // calloc() and malloc()
    
    return 0;
}



