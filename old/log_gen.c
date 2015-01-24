#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int t_max = 0;
int n_x = 0;
int n_y = 0;
int n_z = 0;

int main()
{
    char *line;
    char line2[100];
    int i = 0;
    
    // read source code
    FILE *filein;
    filein = fopen("log_i7", "r");
    
    if (filein == NULL)
    {
        printf(">> cannot find source file");
        return 1;
    }
    
    
    
    // write to new file
    FILE *fileout;
    fileout = fopen("log_i7_post", "w+");
        
    if (fileout == NULL)
    {
        printf(">> cannot write output code");
        return 1;
    }    
    
    size_t size = 100;
    
    i = 5;
    
    while(-1 != getline( &line, &size, filein))
    {
        printf("%d : %s ", i, line);
        
        sprintf(line2, "%d    %s", i, line);
        fputs(line2, fileout);
        
        i++;
    }
    
   
    
    
    /*
    
    for(i = 5; i < 301; i++)
    {
        
        // write domain size
        sprintf(line, "%d    ", i);
        fputs(line, fileout);
        
    
    }
    
    */
    
    fclose(fileout);
    
    fclose(filein);
     
    printf(">> file created !\n");
    
    return 0;
}




