/* codes for reading and writing sac files
   By Liang Zhao, IGG, CAS
   June 28, 2006
*/

#include "sachead.h"
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/file.h>
#include <malloc.h>

void read_sac(sacfile,h,ss)
char sacfile[80];
SACHEAD *h;
float *ss;
{
    float *s=ss;
    int vector_size;
    FILE * fp;
    if ((fp = fopen( sacfile, "rb")) == NULL)
    {
        printf( "cannot open sac file: %s \n", sacfile);
    }
    vector_size = sizeof(*h);
    fread((char *) h, vector_size,1,fp);
    
    // read in data
    vector_size = h->npts*sizeof(float);

    fread((char *)s, vector_size,1,fp);
    fclose(fp);
}

void read_sac_h(sacfile,h,max)
char sacfile[80];
SACHEAD *h;
float *max;
{
    int vector_size;
    float *s;
    FILE * fp;
    if ((fp = fopen( sacfile, "rb")) == NULL)
    {
        printf( "cannot open sac file: %s \n", sacfile);
    }

    vector_size = sizeof(*h);
    fread((char *)h, vector_size,1,fp);
/*    printf("npts=%d size=%d", h->npts, vector_size);
    printf(" b=%-8.3f e=%-8.3f",h->b,h->e);
    printf(" min=%-10.3e max=%-10.3e\n",h->depmin,h->depmax);*/
    // read in data
    vector_size = h->npts*sizeof(float);
    s = (float *) malloc(vector_size);

    fread((char *)s, vector_size,1,fp);

    fclose(fp);

    int i,j;
    for(i=0;i<h->npts;i++){
         if(fabs(s[i]) > *max){
             *max = fabs(s[i]);
             j=i;
         }
    }
}

void h_v(sacfile,h)
char sacfile[80];
SACHEAD *h;
{
    int vector_size;
    float *s;
    FILE * fp;
    if ((fp = fopen( sacfile, "rb")) == NULL)
    {
        printf( "cannot open sac file: %s \n", sacfile);
    }

    vector_size = sizeof(*h);
    fread((char *)h, vector_size,1,fp);
/*    printf("npts=%d size=%d", h->npts, vector_size);
    printf(" b=%-8.3f e=%-8.3f",h->b,h->e);
    printf(" min=%-10.3e max=%-10.3e\n",h->depmin,h->depmax);*/
    // read in data
    vector_size = h->npts*sizeof(float);
    s = (float *) malloc(vector_size);

    fread((char *)s, vector_size,1,fp);

    fclose(fp);
}


void write_sac(sacfile,h,ss)
char sacfile[80];
SACHEAD h;
float *ss;
{
    float *s=ss;
    int vector_size;
    FILE * fp;

    if ((fp = fopen( sacfile, "wb")) == NULL)
    {
        printf( "cannot write sac file: %s \n", sacfile);
    }
    vector_size = sizeof(h);

    fwrite((char *) &h, vector_size,1,fp);
    
    vector_size = h.npts*sizeof(float);

//    printf("npts=%d vectorsize=%d \n", h.npts, vector_size);
    fwrite((char *)s, vector_size,1,fp);
    fclose(fp);
}
