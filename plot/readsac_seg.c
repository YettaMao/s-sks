/*------------------------------------------------------------------
 *   a series of sac files ---> GMT-Plot coordinate system
 *
 *   copyright (c) 2006 by Liang Zhao, IGG,CAS
 *   June 28, 2006
 *   20141118 read a segment of sac file
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
    float t_b, t_e;

    char sacfile[50];

    SACHEAD *h = &sac_null ;
    if(argc != 4) 
    {
        printf("read_sac - Read an sac file\n\n");
        printf("Usage: read_sac_seg file_name t_b t_e\n");
    }
    sscanf(&argv[1][0], "%s",sacfile);
    sscanf(&argv[2][0], "%f", &t_b);
    sscanf(&argv[3][0], "%f", &t_e);
//    printf("sacfile=%s time window= %f ~ %f \n", sacfile, t_b, t_e);

    h_v(sacfile,h);

    num = h->npts;
    printf("num=%d \n",num);

/* output to component files */
    float time,tempf;
    FILE * fopu, *fop;
    float maxv = -999.0;
    s = (float *) malloc(num*sizeof(float));
    read_sac(sacfile,h,s);

        fopu = fopen("sac.xy","wt");     
        fprintf(fopu,">> %d\n", num);
        for(i=0;i<num;i++){
               time = h->b - h->o + i * h->delta;
               if(time>=t_b && time <=t_e){
                  fprintf(fopu,"%-9.3f %-12.6e\n", time, s[i]);               
                  tempf = s[i];
                  if(tempf < 0) tempf *= -1;
                  if(tempf> maxv) maxv = tempf;
               }
        }
        fclose(fopu);
//        printf("max value= %e \n", maxv);
        fop = fopen("maxv.txt","wt");     
        fprintf(fop,"%-12.6e\n", maxv);
        fclose(fop);
        return maxv;
}
