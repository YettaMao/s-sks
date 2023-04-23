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

/* about transfer (x1,x2,x3) to (z,r,t) system
     x2 ^                     R ^  
        |                       |
        |                       |
        |                       |
        |-----------> x1        |-----------> T
   Since x1->R x2->t, t = x2*-1
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
    int index,i,j;
    int vector_size;
    int file_num;
    FILE *fop;
    float *s;
    float maxv;
    float amplify=5;

    char sacfile[100],t_index[20];
    char comp[3][20],plot[3][80];

    strcpy(comp[0],"u");
    strcpy(comp[1],"v");
    strcpy(comp[2],"w");
    strcpy(plot[0],"u.xy");
    strcpy(plot[1],"v.xy");
    strcpy(plot[2],"w.xy");

    SACHEAD *h = &sac_null ;
    if(argc != 2) 
    {
        printf("sac2plot - Plot seismograms on a GMT map \n\n");
        printf("Usage: sac2plot file_number \n");
    }
    sscanf(&argv[1][0], "%d", &file_num);

    maxv=-999;
    float g1,g2,unit;
    float a_max[1];
    a_max[0] = 0.;
    int num;
    for(index=0;index<file_num;index++){
        for(j=0;j<3;j++){
           fop = fopen("tt.txt","wt");
           fprintf(fop,"%d",index);
           fclose(fop);   
           fop = fopen("tt.txt","rt");
           fscanf(fop,"%s",t_index);
           fclose(fop); 

           strcpy(sacfile,t_index);
           if(j==0) strcat(sacfile,".ce");
           if(j==1) strcat(sacfile,".cn");
           if(j==2) strcat(sacfile,".cz");
           read_sac_h(sacfile,h,a_max);
//           if(fabs(h->depmin)>maxv) maxv= fabs(h->depmin);
//           if(fabs(h->depmax)>maxv) maxv= fabs(h->depmax);
           if(a_max[0] > maxv) {
               maxv = a_max[0];
               printf("index=%-3d j=%-3d a_max=%-10.3e sacfile=%s\n",index,j,a_max[0],sacfile);
           }
           if(index == 0) g1= h->gcarc;
           if(index == 1) g2= h->gcarc;
        }
    }
    unit = fabs(g2 - g1);
    num = h->npts;
    printf("maxv=%-10.3e unit=%-8.3f\n",maxv,unit);
    if(maxv>10) printf("Warning: a big value!\n");
    if(maxv <=0) printf("Warning: an error max value!\n");

/* output to component files */
    float scale,time;
    FILE * fopu;
    scale = unit/maxv*amplify;
    s = (float *) malloc(num*sizeof(float));
    for(j=0;j<3;j++){
        fopu = fopen(plot[j],"wt");     
        for(index=0;index<file_num;index++){
           fop = fopen("tt.txt","wt");
           fprintf(fop,"%d",index);
           fclose(fop);   
           fop = fopen("tt.txt","rt");
           fscanf(fop,"%s",t_index);
           fclose(fop); 

           strcpy(sacfile,t_index);
           if(j==0) strcat(sacfile,".ce");
           if(j==1) strcat(sacfile,".cn");
           if(j==2) strcat(sacfile,".cz");
           read_sac(sacfile,h,s);
           fprintf(fopu,">> line %f %d\n", h->gcarc,num);
           for(i=0;i<num;i++){
               time = h->b - h->o + i * h->delta;
               if(j==1) fprintf(fopu,"%-9.3f %-11.6e\n", time, h->gcarc - s[i]*scale);               
               else fprintf(fopu,"%-9.3f %-11.6e\n", time, h->gcarc + s[i]*scale);
           }
        }
        fclose(fopu);
    }
}
