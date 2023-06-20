// Sample program for matrix inversion using blockwise inversion methodology.  
// It is for performing inversion of random generated matrix. user programs as per requirement.  
// This sample program `testinvertor.c' can call the different invertor functions.  Keep only one as uncommented in the list of include functions.

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

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

//Keep only one of the following include function as uncommented for performing inversion by that method.
//#include "invertor_by_a.c"
#include "invertor_inplace_by_a.c"
//#include "invertor_by_ad.c"
//#include "invertor_by_prll.c"

int testfunc(double** , int , int , double** , int , int );

int matmul( double** , int , int , double** , int , int , double** , int , int );

int matmulthree(double** mata, int ma, int na, double** matb, int mb, int nb, double** matc, int mc, int nc, double** matres, int mres, int nres);

int matmulfive(double** mata, int ma, int na, double** matb, int mb, int nb, double** matc, int mc, int nc, double** matd, int md, int nd, double** mate, int me, int ne, double** matres, int mres, int nres);

int invertmatfive(double**, double** );
int invertmatsix(double**, double** );
//int invertmateleven(double**, double** );
int scalarmul(double **, int, int, double );

int mataddition( double** , int , int , double** , int , int , double** , int , int );
int matsubtraction( double** , int , int , double** , int , int , double** , int , int );


double randomrange(int ip, double min, double max)
{
        double val;
        int pw;
        pw = (int)log10((double)ip*1.0);
        val= (double)(ip/(pow(10.0,pw+1.0)));
        //Scaling to the range from 0 to 1 to min to max:
        val= min + val*(max-min);
        return val;
}
                
int main()
{
	int n=23;
	int i,j,k;

	double **matsmallres;

	double **p1;
	double **p2;
	
	//For printing time
	time_t timer,t;
	char buffer[26];
	struct tm* tm_info;
	clock_t start, end;
	double cpu_time_used;

	printf("\n%d by %d matrix inversion\n", n , n);
	p1=(double **)calloc(n,sizeof(*p1));
	for(i=0;i<n;i++) p1[i]=(double *) calloc(n,sizeof(p1[i]));

	//p2 here is for multiplication to verify inversion	
	p2=(double **)calloc(n,sizeof(*p2));
	for(i=0;i<n;i++) p2[i]=(double *) calloc(n,sizeof(p2[i]));
	
	//for(i=0;i<n;i++) for(j=0;j<n;j++) p1[i][j]=(double)rand();

	matsmallres=(double **)calloc(n,sizeof(*matsmallres));
	for(i=0;i<n;i++) matsmallres[i]=(double *) calloc(n,sizeof(matsmallres[i]));

	printf("\ninput matrix of order %d * %d\n",n,n);
	srand((unsigned)time(&t)); //Initializing random number generator
	//srand((unsigned)7); //Initializing random number generator
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			p1[i][j]= randomrange(rand(),-10.0,10.0); // rand()/(RAND_MAX/range); //range= (upper - lower );
				printf("%lf\t",p1[i][j]);
		}
		printf("\n");
	}

	//We print time here
	timer = time(NULL);
	tm_info = localtime(&timer);

	strftime (buffer, 100, "%Y-%m-%d %H:%M:%S.000", localtime (&timer));
	printf("\nBegining time\n");
	puts(buffer);
	
	start = clock();	
        invertmat(n,p1,matsmallres);
        end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("\n CPU time used: %f seconds\n",cpu_time_used);

	timer = time(NULL);
	tm_info = localtime(&timer);

	//strftime(buffer, 26, "%Y-%m-%d %H:%M:%S", tm_info);
	strftime (buffer, 100, "%Y-%m-%d %H:%M:%S.000", localtime (&timer));
	printf("\nEnd time\n");
	puts(buffer);
	printf("\n---------------------\n");

//Uncomment the following section to print the inverted matrix and to verify the inversion.

	printf("\ninverted matrix\n");
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++) printf("%lf\t",matsmallres[i][j]);
		printf("\n");
	}


	matmul(p1,n,n,matsmallres,n,n,p2,n,n);
	
	
	printf("\nverification of inversion\n");
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++) 
		//printf("%lf\t",p2[i][j]);
		printf("%lf\t",p2[i][j]);
		//{
			//if(p2[i][j]>=1.000100 || p2[i][j]<=-0.0001000) printf("\n%d, %d \t%lf",i,j,p2[i][j]);
			//if(p2[i][j]<=0.000900 && p2[i][j]>=0.0001000) printf("\n%d, %d \t%lf",i,j,p2[i][j]);
		//}
		
		printf("\n");
	}
	printf("\nverification of inversion completed");


	for(i=0;i<n;i++) free(p1[i]); free(p1);
	for(i=0;i<n;i++) free(matsmallres[i]); free(matsmallres);
	for(i=0;i<n;i++) free(p2[i]); free(p2);
	return 0;
}

int testfunc(double** mata, int ma, int na, double** matres, int mres, int nres)
{
	int i1, j1;	
	for(i1=0;i1<ma; i1++)
	{
		for(j1=0; j1<na; j1++)
			matres[i1][j1]=mata[i1][j1]*mata[i1][j1];
	}
	return 1;
}


int matmul( double** mata, int ma, int na, double** matb, int mb, int nb, double** matres, int mc, int nc)
{
	int i,j,k,l;
	if(na!=mb) return 0;

	for(i=0;i<ma;i++)
		for(j=0;j<nb;j++)
			for(matres[i][j]=k=0;k<na;k++)
				matres[i][j]+=mata[i][k]*matb[k][j];
	return 1;
}

int matmulthree(double** mata, int ma, int na, double** matb, int mb, int nb, double** matc, int mc, int nc, double** matres, int mres, int nres)
{
	int i;
	double **tempmat;

	tempmat = (double **) malloc( ma * sizeof(double *));
	for(i=0;i<ma;i++) tempmat[i]=(double *) malloc(nb * sizeof(double));

	i=matmul(mata, ma, na, matb, mb, nb, tempmat, ma, nb);
	if(i==0) return 0;

	i=matmul(tempmat, ma, nb, matc, mc, nc, matres, ma, nc);
	if(i==0) return 0;

	for(i=0;i<ma;i++) free(tempmat[i]);
	free(tempmat);
	return 1;
}

int matmulfive(double** mata, int ma, int na, double** matb, int mb, int nb, double** matc, int mc, int nc, double** matd, int md, int nd, double** mate, int me, int ne, double** matres, int mres, int nres)
{
	int i;
	double **tempmat;

	tempmat = (double **) malloc( ma * sizeof(double *));
	for(i=0;i<ma;i++) tempmat[i]=(double *) malloc(nc * sizeof(double));

	i=matmulthree(mata, ma, na, matb, mb, nb, matc, mc, nc, tempmat, ma, nc);
	if(i==0) return 0;

	i=matmulthree(tempmat, ma, nc, matd, md, nd, mate, me, ne,  matres, ma, ne);
	if(i==0) return 0;

	for(i=0;i<ma;i++) free(tempmat[i]);
	free(tempmat);
	return 1;
}

int scalarmul(double** mata, int ma, int na, double x)
{
	int i,j;
	for(i=0;i<ma;i++)
		for(j=0;j<na;j++)
			mata[i][j]*=x;
	return 1;
}



int mataddition( double** mata, int ma, int na, double** matb, int mb, int nb, double** matres, int mc, int nc)
{
	int i,j,k,l;
	if(ma!=mb) return 0;

	for(i=0;i<ma;i++)
		for(j=0;j<na;j++)
				matres[i][j]=mata[i][j]+matb[i][j];
	return 1;
}

int matsubtraction( double** mata, int ma, int na, double** matb, int mb, int nb, double** matres, int mc, int nc)
{
	int i,j,k,l;
	if(ma!=mb) return 0;

	for(i=0;i<ma;i++)
		for(j=0;j<na;j++)
				matres[i][j]=mata[i][j]-matb[i][j];
	return 1;
}



