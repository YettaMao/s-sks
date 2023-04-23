/*------------------------------------------------------------------
 *   a series of sac files ---> GMT-Plot coordinate system
 *
 *   copyright (c) 2006 by Liang Zhao, IGG,CAS
 *   June 28, 2006
 *------------------------------------------------------------------*/
/*
 * read sacfiles and plot the lines on a map
 *  ................-------------------------------------------
 *  |quake time     |kztime      |        |
 *  |<-   omarker ->|            |        |
 *  |<--     begin time       -->|        |
 *  |<--           end time            -->|
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

main(int argc, char ** argv)
{
    int i,num;
    int vector_size;
    float *s;

    char sacfile[50];

    SACHEAD *h = &sac_null ;
    if(argc != 2) 
    {
        printf("read_sac - Read an sac file\n\n");
        printf("Usage: read_sac file_name \n");
    }
    strcat(sacfile,argv[1]);
    printf("sacfile=%s \n", sacfile);

    h_v(sacfile,h);

    num = h->npts;
    printf("num=%d \n",num);

/* output to component files */
    float time;
    FILE * fopu;
    s = (float *) malloc(num*sizeof(float));
    read_sac(sacfile,h,s);

        fopu = fopen("sac.xy","wt");     
        fprintf(fopu,">> %d\n", num);
        for(i=0;i<num;i++){
               time = h->b - h->o + i * h->delta;
               fprintf(fopu,"%-9.3f %-12.6e\n", time, s[i]);               
        }
        fclose(fopu);
}
