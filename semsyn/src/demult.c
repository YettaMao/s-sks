#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


int main(int ac,char **av)
{
  int ntr,npts,i,j,k,check,size_f;
  float *data,*trace;
  int nex,nx, nl, nt, nz;
  int hpn;/*high polynomials number*/
  int nl_skip = 0, kdx = 1;
  FILE *fop, *open_file(), *fop1; char greenfile[50], greenfile1[50];
  FILE *fp; /* Liang 6/9/2006 */
  int COMP = 2;
  int diffile = 0;
  int finish;
  long int offset, rec_len;
  int read_record();

  size_f = sizeof(float);
  
  setpar(ac,av);
  mstpar("nex",       "d", &nex);
  mstpar("hpn",       "d", &hpn);
  getpar("nlskip",    "d", &nl_skip);
  getpar("kdx",       "d", &kdx);
  mstpar("nt",        "d", &nt);
  getpar("greenfile", "s", greenfile);
  getpar("diffile",   "d", &diffile);
  getpar("COMP",      "d", &COMP);
  if(diffile != 0)
      mstpar("greenfile1","s", greenfile1);
  endpar();

  nx = nex*hpn + 1;

  fop  = open_file(greenfile,"rb"); 
  fread(&nz,sizeof(int),1,fop);
  printf("nz=%d\n",nz);

/* output nz for plot, by Liang, June 9, 2006 */
   fp = fopen("dimension.txt","wt");
   fprintf(fp,"%d\n",nz);
   fclose(fp);   

   ntr  = COMP*(nx+nz);
   npts = nt;
   printf("ntr = %d, npts = %d\n",ntr,npts);    

  if( diffile != 1){
      if( (data = (float *)malloc(size_f * ntr * npts)) == NULL){
          fprintf(stderr,"cannot allocate memory for data\n");
          exit(-1);
      }
      check = fread(data,size_f,ntr*npts,fop);
      if (check != ntr*npts) {
        fprintf(stderr,"demult: too little data for specified nx and nt\n");
        exit(1);
      }
      fclose(fop);
  
      fop  = open_file(greenfile,"wb"); 
      trace = (float *)malloc(size_f * ntr); assert(trace != NULL);
/* ZHL: May 10, 2010 debug, modify the following line */
/* Result: It outputs the same waveforms as before */
    /* for (i=0;i<npts;++i) { */
      for (i=0;i<npts;i++) {
        for (j=i,k=0;j<ntr*npts;j+=npts,k++) {
          trace[k] = data[j];
        }
        check = fwrite(trace,size_f,ntr,fop);
        if (check != ntr) {
          fprintf(stderr,"demult: write failed for trace %d\n",i);
          exit(1);
        }
      }
      fclose(fop);
      free(data);
      finish = 0;
   } else {
       fop1 = open_file(greenfile1,"wb");
       fprintf(stderr,"ntr=%d\n",ntr);
       trace = (float *)malloc(size_f * ntr); assert(trace != NULL);
       rec_len = npts * size_f;
/*May 10, 2010 debug, modify the following line */
    /* for (i=0;i<npts;++i) { */
       for (i=0;i<npts;i++) {
           offset = i * size_f +sizeof(int);
           read_record(fop,offset,trace,ntr,size_f,rec_len);

           check = fwrite(trace,size_f,ntr,fop1);
           if (check != ntr) {
               fprintf(stderr,"demult: write failed for trace %d\n",i);
               exit(1);
           }
        }
        fclose(fop1);
        fclose(fop);
        finish = 1;
   }

  return finish;
}
