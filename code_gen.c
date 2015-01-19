#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int t_max = 2000;
int n_x = 10;
int n_y = 10;
int n_z = 10;


int main()
{
    char line[100];
    size_t len = 0;
    ssize_t read = 0;
    int i = 0;
    
    // read source code
    FILE *filein;
    filein = fopen("lwacm_raw", "r");
    
    if (filein == NULL)
    {
        printf(">> cannot find source file lwacm.c");
        return 1;
    }
    
    // read all file into buffer
    fseek(filein, 0, SEEK_END);
    
    long pos = ftell(filein);
    
    fseek(filein, 0, SEEK_SET);
    
    char *buffer = malloc(pos);
    fread(buffer, pos, 1, filein);
    
    fclose(filein);
    
    // define the run script to run all programs
    FILE *filerun;
    filerun = fopen("run.sh", "w+");    
    char runname[100];
    
    
    
    // output files
    printf(">> start generating source code\n");
    
    long domain = 0;
    char filename[100];
    
    
    for(i = 1; i < 100; i++)
    {
        t_max = 2000/i;
        
        n_x = i;
        n_y = i;
        n_z = i;
        
        // calculate domain
        domain = t_max * n_x * n_y * n_z;
        
        // make new filename
        sprintf(filename, "lwacm_run_%ld.c", domain);
        sprintf(runname, "./lwacm_run_%ld\n", domain);
        
        // write to new file
        FILE *fileout;
        fileout = fopen(filename, "w+");
        
        if (fileout == NULL)
        {
            printf(">> cannot write output code");
            return 1;
        }
        
        // write domain size
        sprintf(line, "\n#define  T_MAX  %d\n", t_max);
        fputs(line, fileout);
        
        sprintf(line, "#define  N_X  %d\n", n_x);
        fputs(line, fileout);
        
        sprintf(line, "#define  N_Y  %d\n", n_y);
        fputs(line, fileout);
        
        sprintf(line, "#define  N_Z  %d\n", n_z);
        fputs(line, fileout);
        
        // write remaining code
        fputs(buffer, fileout);
        
        fclose(fileout);
        
        fputs(runname, filerun);
    
    }
    
    
    fclose(filerun);
     
    printf(">> all files are created !\n");
    
    return 0;
}




