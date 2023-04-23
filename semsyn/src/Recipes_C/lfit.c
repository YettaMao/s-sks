/* adapted from recipes C */
static double sqrarg;
#define SQR(a) (sqrarg=(a),sqrarg*sqrarg)

void lfit(y,sig,ndata,a,ma,ia,covar,chisq,funcs)
int ndata,ma,ia[];
double y[],sig[],a[],**covar,*chisq, **funcs;
{
	int i, j, k, l, m, mfit=0;
	double ym,wt,sum,sig2i,**beta,*afunc;
	void gaussj(),covsrt(),nrerror(),free_dmatrix(),free_dvector();
	double **dmatrix(),*dvector();

	beta=dmatrix(1,ma,1,1);
	afunc=dvector(1,ma);
	for (j=1;j<=ma;j++) 
		if(ia[j]) mfit++;
	if (mfit == 0) nrerror("Bad LISTA permutation in LFIT-2");
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=mfit;k++) covar[j][k]=0.0;
		beta[j][1]=0.0;
	}
	for (i=1;i<=ndata;i++) {
	        for(j=1; j<=ma; j++)
		    afunc[j]=funcs[i][j];

		ym=y[i];
		if (mfit < ma){
			for(j=1; j<=ma; j++)
			     if(!ia[j]) ym -= a[j]*afunc[j];
		}
		sig2i=1.0/SQR(sig[i]);
		for (j=0,l=1;l<=mfit;l++) {
			if(ia[l]){
			    wt=afunc[l]*sig2i;
			    for (j++,k=0, m=1;m<=l;m++)
				if(ia[m])covar[j][++k] += wt*afunc[m];
			    beta[j][1] += ym*wt;
			}
		}
	}
	for (j=2;j<=mfit;j++)
		for (k=1;k<=j-1;k++)
			covar[k][j]=covar[j][k];
	gaussj(covar,mfit,beta,1);
	for (j=0,l=1;l<=ma;l++) a[l]=beta[++j][1];
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
	        for(j=1; j<=ma; j++)
		    afunc[j]=funcs[i][j];

		for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
		*chisq += SQR((y[i]-sum)/sig[i]);
	}
	covsrt(covar,ma,ia,mfit);
	free_dvector(afunc,1,ma);
	free_dmatrix(beta,1,ma,1,1);
}

#undef SQR
