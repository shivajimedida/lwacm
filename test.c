#include <time.h>
#include <stdio.h>

int main()
{
   int x = 0;
   int y = 0;
   int z = 0;
   
   int p[10][10];
   
   printf(">> init array...\n");
   
   for( x = 0; x < 10; x++)
   {
	for( y = 0; y < 10; y++)
	{
	    p[x][y] = 0;
        }
   }
   
   
   printf(">> test array...\n");
   for( x = 0; x < 10; x++)
	for( y = 0; y < 10; y++)
	  {
	    p[x][y] = 1;
	  }		
   
   
   printf(">> output array\n");
   for( x = 0; x < 10; x++)
   {
	printf("   ");
	
	for( y = 0; y < 10; y++)
	{
	    printf("p[%d][%d] = %d, ", x, y, p[x][y] );
        }
        
        printf("\n");
   }
   
   
   

   return(0);
}


