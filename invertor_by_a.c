// General inversion function for matrix inversion using blockwise inversion methodology.  
// The program is meant to be adopted with user programs as per requirement.  
// A sample program `testinvertor.c' calls this `invertor_by_a.c' function, with random genrated matrix of different orders.

// To use: call the function "invermat" with the arguments: 
// 1. order the matrix, 2. input matrix (as 2 dimensional array) and 3. output matrix (as 2 dimensional array)
// The return value is 1 for successful calculation of inverse.

// Author: R. Thiru Senthil.
// The Institute of Mathematical Sciences, 
// IV Cross St, CIT Campus, Taramani, Chennai 600113, Tamil Nadu, India.
// Email: rtsenthil@imsc.res.in
// Presented at: ICHEP 2022
// Kindly cite as: 
// 1. Inspire Link: https://inspirehep.net/literature/2619671
// R.~Thiru Senthil, ``Invertor - Program to compute exact inversion of large matrices,'' PoS \textbf{ICHEP2022}, 1129 (2022)
// doi:10.22323/1.414.1129
// 2. Inspire Link: https://inspirehep.net/literature/2660850
// R. Thiru Senthil, ``Blockwise inversion and algorithms for inverting large partitioned matrices,'' [arXiv:2305.11103 [math.NA]].(Submitted)

// The invertor project details with downloads are available in the webpage: https://www.imsc.res.in/~rtsenthil/invertor.html
// and in github page: https://github.com/rthirusenthil/invertor

#include<stdlib.h>

int invertmatone(double** mata, double** inverta);
int invertmattwo(double** mata, double** inverta);
int invertmatthree(double** mata, double** inverta);
int invertmatfour(double** mata, double** inverta);

int invertblocks(int n, double** mata , double** inverta);

int invertor_by_a(int n, double** mata, double** inverta);

int byamatmul( double** mata, int ma, int na, double** matb, int mb, int nb, double** matres, int mc, int nc);
int byamatmulthree(double** mata, int ma, int na, double** matb, int mb, int nb, double** matc, int mc, int nc, double** matres, int mres, int nres);
int byascalarmul(double** mata, int ma, int na, double x);
int byamatsubtraction( double** mata, int ma, int na, double** matb, int mb, int nb, double** matres, int mc, int nc);

int invertmat(int n, double** mata, double** inverta)
{
        int invertstatus;
        int i,j, order;

        if(n<=0) return 0;
/*
        for(i=0;i<n;i++)
                for(j=0;j<n;j++)
                        inverta[i][j]=mata[i][j];
        order = n;
*/
        invertstatus = invertor_by_a(n, mata, inverta);
	if(invertstatus==0)
	{
		printf("\nUnable to invert matrix of order %d\n",n);
	}
        return invertstatus;
}


int invertor_by_a(int n, double** mata, double** inverta)
{
	int invertstatus;
	int i,j,k;
	
	if(n<=0) return 0;
	switch(n)
	{
		case 1:
			invertstatus=invertmatone(mata, inverta);
			break;
		case 2:
			invertstatus=invertmattwo(mata, inverta);
			break;
		case 3:
			invertstatus=invertmatthree(mata, inverta);
			break;
		case 4:
			invertstatus=invertmatfour(mata, inverta);
			break;
		/*case 5:
			invertstatus=invertmatfive(mata, inverta);
			break;
		case 6:
			invertstatus=invertmatsix(mata, inverta);
			break; */
		default:
			invertstatus=invertblocks(n, mata, inverta);
			break;
	}

	return invertstatus;
}


int invertblocks(int n, double** mata, double** inverta)
{
	int invertstatus;
	int i,j,k;
	

	double **mate, **matf, **matg, **math, **mats;
	int me, ne, mf, nf, mg, ng, mh, nh, ms, ns;
	
	double **matinve, **matinvs;
	int minve, ninve, minvs, ninvs;
	
	double **matsol1, **matsol2, **matsol3, **mattemp1, **mattemp2;
	int msol1, nsol1, msol2, nsol2, msol3, nsol3, mtemp1, ntemp1, mtemp2, ntemp2;

	
//	printf("\n We entered n by n inversion function\n");
		
	//Preparing the E, F, G, H Matrices
	me=n/2;
	ne=n/2;
	mf=me;
	nf=n-ne;
	mg=n-me;
	ng=ne;
	mh=n-me;
	nh=n-ne;
	
	mate=(double **) malloc(me * sizeof(*mate));
	for(i=0; i<me; i++) mate[i]=(double *)malloc(ne * sizeof(mate[i]));

	matf=(double **) malloc(mf * sizeof(*matf));
	for(i=0; i<mf; i++) matf[i]=(double *)malloc(nf * sizeof(matf[i]));

	matg=(double **) malloc(mg * sizeof(*matg));
	for(i=0; i<mg; i++) matg[i]=(double *)malloc(ng * sizeof(matg[i]));

	math=(double **) malloc(mh * sizeof(*math));
	for(i=0; i<mh; i++) math[i]=(double *)malloc(nh * sizeof(math[i]));

	for(i=0; i<me; i++)
	{
		for(j=0; j<ne; j++) mate[i][j]=mata[i][j];
		for(j=ne; j<n; j++) matf[i][j-ne]=mata[i][j];
	}
	for(i=me;i<n;i++)
	{
		for(j=0;j<ng;j++) matg[i-me][j]=mata[i][j];
		for(j=ng;j<n;j++) math[i-me][j-ng]=mata[i][j];
	}
	
/*	
	
	printf("\nWe are printint individual matrices here\n");
	printf("Matrix E\n");
	for(i=0;i<me;i++)
	{
		for(j=0;j<ne;j++) printf("\t%lf",mate[i][j]);
		printf("\n");
	}
	printf("Matrix F\n");
	for(i=0;i<mf;i++)
	{
		for(j=0;j<nf;j++) printf("\t%lf",matf[i][j]);
		printf("\n");
	}
	printf("Matrix G\n");
	for(i=0;i<mg;i++)
	{
		for(j=0;j<ng;j++) printf("\t%lf",matg[i][j]);
		printf("\n");
	}
	printf("Matrix H\n");
	for(i=0;i<mh;i++)
	{
		for(j=0;j<nh;j++) printf("\t%lf",math[i][j]);
		printf("\n");
	}


*/

	//Calculating E^-1 and freeing E	
	minve = me;
	ninve = ne;
	
	matinve=(double **) malloc(minve * sizeof(*matinve));
	for(i=0; i<minve; i++) matinve[i]=(double *)malloc(ninve * sizeof(matinve[i]));
	
	invertstatus=invertor_by_a(me, mate, matinve);
	if(invertstatus==0)
        {
                printf("\nUnable to invert matrix of order me = %d\n",me);
        }

//	printf("\n we finished mate inversion\n");	
	for(i=0;i<me;i++) free(mate[i]); free(mate);
	
	//Calculating S 
	ms=mh;
	ns=nh;

	mats=(double **) malloc(ms * sizeof(*mats));
	for(i=0; i<ms; i++) mats[i]=(double *)malloc(ns * sizeof(mats[i]));
		
	mtemp1=ms;
	ntemp1=ns;
	
	mattemp1=(double **) malloc(mtemp1 * sizeof(*mattemp1));
	for(i=0; i<mtemp1; i++) mattemp1[i]=(double *)malloc(ntemp1 * sizeof(mattemp1[i]));

	byamatmulthree(matg, mg, ng, matinve, minve, ninve, matf, mf, nf, mattemp1, mtemp1, ntemp1);

//	printf("\n we finished multiplication of g e^-1 f\n");	
	byamatsubtraction(math, mh, nh, mattemp1, mtemp1, ntemp1, mats, ms, ns);
	
//	printf("\n we finished mats calculation\n");
		
	for(i=0;i<mtemp1;i++) free(mattemp1[i]); free(mattemp1);	


	//Calculating S^-1 and freeing S
	minvs=ms;
	ninvs=ns;
	
	matinvs=(double **) malloc(minvs * sizeof(*matinvs));
	for(i=0; i<minvs; i++) matinvs[i]=(double *)malloc(ninvs * sizeof(matinvs[i]));

	invertstatus=invertor_by_a(ms,mats, matinvs);
	if(invertstatus==0)
        {
                printf("\nUnable to invert matrix of order ms= %d\n",ms);
        }

	//S^-1 calculated so we free S
	for(i=0;i<ms;i++) free(mats[i]); free(mats);
	
	//At this time we have E^-1, S^-1.  This S^-1 is Solution 4.
	//We have to prepare other three solutions.
	//Preparing Solution 3
	msol3=mg;
	nsol3=ng;
	
	matsol3=(double **) malloc(msol3 * sizeof(*matsol3));
	for(i=0; i<msol3; i++) matsol3[i]=(double *)malloc(nsol3 * sizeof(matsol3[i]));	

	byamatmulthree(matinvs, minvs, ninvs, matg, mg, ng, matinve, minve, ninve, matsol3, msol3, nsol3);
	
	byascalarmul(matsol3, msol3, nsol3, (double) -1.0);

	//At this time We have S^-1 which is solution 4 and Solution 3.
	//Preparing Solution 2
	msol2=mf;
	nsol2=nf;
	
	matsol2=(double **) malloc(msol2 * sizeof(*matsol2));
	for(i=0; i<msol2; i++) matsol2[i]=(double *)malloc(nsol2 * sizeof(matsol2[i]));	

	byamatmulthree(matinve, minve, ninve, matf, mf, nf, matinvs, minvs, ninvs, matsol2, msol2, nsol2);

	byascalarmul(matsol2, msol2, nsol2, (double) -1.0);

	//Preparing Solution 1
	msol1=me;
	nsol1=ne;
	
	matsol1=(double **) malloc(msol1 * sizeof(*matsol1));
	for(i=0; i<msol1; i++) matsol1[i]=(double *)malloc(nsol1 * sizeof(matsol1[i]));	
	
	mtemp1=msol1;
	ntemp1=nsol1;	
	
	mattemp1=(double **) malloc(mtemp1 * sizeof(*mattemp1));
	for(i=0;i<mtemp1;i++) mattemp1[i]=(double *)malloc(ntemp1 * sizeof(mattemp1[i]));

	byamatmulthree(matinve, minve, ninve, matf, mf, nf, matsol3, msol3, nsol3, mattemp1, mtemp1, ntemp1);
	
	byamatsubtraction(matinve, minve, ninve, mattemp1, mtemp1, ntemp1, matsol1, msol1, nsol1);
		
	for(i=0;i<mtemp1;i++) free(mattemp1[i]); free(mattemp1);


	//Preparing Full solution
	for(i=0;i<me;i++)
	{
		for(j=0;j<ne;j++) inverta[i][j]=matsol1[i][j];
		for(j=ne;j<n;j++) inverta[i][j]=matsol2[i][j-ne];
	}
	for(i=me;i<n;i++)
	{
		for(j=0;j<ne;j++) inverta[i][j]=matsol3[i-me][j];
		for(j=ne;j<n;j++) inverta[i][j]=matinvs[i-me][j-ne];
	}

	for(i=0;i<minvs;i++) free(matinvs[i]); free(matinvs);
	for(i=0;i<msol1;i++) free(matsol1[i]); free(matsol1);	
	for(i=0;i<msol2;i++) free(matsol2[i]); free(matsol2);	
	for(i=0;i<msol3;i++) free(matsol3[i]); free(matsol3);	
						 	
	return 1;
	
}

int invertmatone(double** mata, double** inverta)
{
	int n=1;
	double modmata;
	modmata=mata[0][0];
	
	if(modmata==0) 
	{
		printf("\nUnable to invert matrix of order 1\n");
		return 0;
	}

	inverta[0][0]=1/modmata;
	
	return 1;
}

int invertmattwo(double** mata, double** inverta)
{
	int n=2;

	double modmata;
	double a11=mata[0][0], a12=mata[0][1];
	double a21=mata[1][0], a22=mata[1][1];
	
	modmata=(-(a12*a21) + a11*a22);
	
	if(modmata==0) 
	{
		printf("\nUnable to invert matrix of order 2\n");
		return 0;
	}

	
	inverta[0][0]=a22/modmata;
	inverta[0][1]=-(a12/modmata);
	inverta[1][0]=-(a21/modmata);
	inverta[1][1]=a11/modmata;
	
	return 1;
}

int invertmatthree(double** mata, double** inverta)
{
	int n=3;

	double modmata;
	double a11=mata[0][0], a12=mata[0][1], a13=mata[0][2];
	double a21=mata[1][0], a22=mata[1][1], a23=mata[1][2];
	double a31=mata[2][0], a32=mata[2][1], a33=mata[2][2];
	
	modmata=(-(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33);
	
	if(modmata==0) 
	{
		printf("\nUnable to invert matrix of order 3\n");
		return 0;
	}

	
	inverta[0][0]=(-(a23*a32) + a22*a33)/modmata;
	
	inverta[0][1]=(a13*a32 - a12*a33)/modmata;
	
	inverta[0][2]=(-(a13*a22) + a12*a23)/modmata;
	
	inverta[1][0]=(a23*a31 - a21*a33)/modmata;
	
	inverta[1][1]=(-(a13*a31) + a11*a33)/modmata;
	
	inverta[1][2]=(a13*a21 - a11*a23)/modmata;
	
	inverta[2][0]=(-(a22*a31) + a21*a32)/modmata;
	
	inverta[2][1]=(a12*a31 - a11*a32)/modmata;
	
	inverta[2][2]=(-(a12*a21) + a11*a22)/modmata;

	return 1;
}

int invertmatfour(double** mata, double** inverta)
{
	int n=4;

	double modmata;
	double a11=mata[0][0], a12=mata[0][1], a13=mata[0][2], a14=mata[0][3];
	double a21=mata[1][0], a22=mata[1][1], a23=mata[1][2], a24=mata[1][3];
	double a31=mata[2][0], a32=mata[2][1], a33=mata[2][2], a34=mata[2][3];
	double a41=mata[3][0], a42=mata[3][1], a43=mata[3][2], a44=mata[3][3];
	
	
	modmata=(a14*a23*a32*a41 - a13*a24*a32*a41 - a14*a22*a33*a41 + a12*a24*a33*a41 + a13*a22*a34*a41 - a12*a23*a34*a41 - a14*a23*a31*a42 + a13*a24*a31*a42 + a14*a21*a33*a42 - a11*a24*a33*a42 - a13*a21*a34*a42 + a11*a23*a34*a42 + a14*a22*a31*a43 - a12*a24*a31*a43 - a14*a21*a32*a43 + a11*a24*a32*a43 + a12*a21*a34*a43 - a11*a22*a34*a43 - a13*a22*a31*a44 + a12*a23*a31*a44 + a13*a21*a32*a44 - a11*a23*a32*a44 - a12*a21*a33*a44 + a11*a22*a33*a44);
	
	if(modmata==0) 
	{
		printf("\nUnable to invert matrix of order 4\n");
		return 0;
	}


	inverta[0][0]=(-(a24*a33*a42) + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 + a22*a33*a44)/modmata;
	
	inverta[0][1]=(a14*a33*a42 - a13*a34*a42 - a14*a32*a43 + a12*a34*a43 + a13*a32*a44 - a12*a33*a44)/modmata;
	
	inverta[0][2]=(-(a14*a23*a42) + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 + a12*a23*a44)/modmata;
	
	inverta[0][3]=(a14*a23*a32 - a13*a24*a32 - a14*a22*a33 + a12*a24*a33 + a13*a22*a34 - a12*a23*a34)/modmata;
	
	inverta[1][0]= (a24*a33*a41 - a23*a34*a41 - a24*a31*a43 + a21*a34*a43 + a23*a31*a44 - a21*a33*a44)/modmata;
	
	inverta[1][1]=(-(a14*a33*a41) + a13*a34*a41 + a14*a31*a43 - a11*a34*a43 - a13*a31*a44 + a11*a33*a44)/modmata;
	
	inverta[1][2]=(a14*a23*a41 - a13*a24*a41 - a14*a21*a43 + a11*a24*a43 + a13*a21*a44 - a11*a23*a44)/modmata;
	
	inverta[1][3]=(-(a14*a23*a31) + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 + a11*a23*a34)/modmata;
	
	inverta[2][0]=(-(a24*a32*a41) + a22*a34*a41 + a24*a31*a42 - a21*a34*a42 - a22*a31*a44 + a21*a32*a44)/modmata;
	
	inverta[2][1]=(a14*a32*a41 - a12*a34*a41 - a14*a31*a42 + a11*a34*a42 + a12*a31*a44 - a11*a32*a44)/modmata;
	
	inverta[2][2]=(-(a14*a22*a41) + a12*a24*a41 + a14*a21*a42 - a11*a24*a42 - a12*a21*a44 + a11*a22*a44)/modmata;
	
	inverta[2][3]=(a14*a22*a31 - a12*a24*a31 - a14*a21*a32 + a11*a24*a32 + a12*a21*a34 - a11*a22*a34)/modmata;
	
	inverta[3][0]=(a23*a32*a41 - a22*a33*a41 - a23*a31*a42 + a21*a33*a42 + a22*a31*a43 - a21*a32*a43)/modmata;
	
	inverta[3][1]=(-(a13*a32*a41) + a12*a33*a41 + a13*a31*a42 - a11*a33*a42 - a12*a31*a43 + a11*a32*a43)/modmata;
	
	inverta[3][2]=(a13*a22*a41 - a12*a23*a41 - a13*a21*a42 + a11*a23*a42 + a12*a21*a43 - a11*a22*a43)/modmata;
	
	inverta[3][3]=(-(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33)/modmata;

	return 1;
}
int byamatmul( double** mata, int ma, int na, double** matb, int mb, int nb, double** matres, int mc, int nc)
{
	int i,j,k,l;
	if(na!=mb) return 0;

	for(i=0;i<ma;i++)
		for(j=0;j<nb;j++)
			for(matres[i][j]=k=0;k<na;k++)
				matres[i][j]+=mata[i][k]*matb[k][j];
	return 1;
}

int byamatmulthree(double** mata, int ma, int na, double** matb, int mb, int nb, double** matc, int mc, int nc, double** matres, int mres, int nres)
{
	int i;
	double **tempmat;

	tempmat = (double **) malloc( ma * sizeof(double *));
	for(i=0;i<ma;i++) tempmat[i]=(double *) malloc(nb * sizeof(double));

	i=byamatmul(mata, ma, na, matb, mb, nb, tempmat, ma, nb);
	if(i==0) return 0;

	i=byamatmul(tempmat, ma, nb, matc, mc, nc, matres, ma, nc);
	if(i==0) return 0;

	for(i=0;i<ma;i++) free(tempmat[i]);
	free(tempmat);
	return 1;
}


int byascalarmul(double** mata, int ma, int na, double x)
{
	int i,j;
	for(i=0;i<ma;i++)
		for(j=0;j<na;j++)
			mata[i][j]*=x;
	return 1;
}



int byamatsubtraction( double** mata, int ma, int na, double** matb, int mb, int nb, double** matres, int mc, int nc)
{
	int i,j,k,l;
	if(ma!=mb) return 0;

	for(i=0;i<ma;i++)
		for(j=0;j<na;j++)
				matres[i][j]=mata[i][j]-matb[i][j];
	return 1;
}



