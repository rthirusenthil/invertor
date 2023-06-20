// Parallel processing inversion program for matrix inversion using blockwise inversion methodology.  
// The program is meant to be adopted with user programs as per requirement.  
// A sample program `testinvertor.c' can call this `invertor_by_prll.c' function, with random genrated matrix of different orders.

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

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <omp.h>
#include<unistd.h>

int invertmat(int n, double** mata, double** inverta);
int invertcases(int n, int apos, double** a, int invapos, double** inverta);

int invertmatone(int pos, double** a, int ipos, double** inverta);
int invertmattwo(int pos, double** a, int ipos, double** inverta);
int invertmatthree(int pos, double** a, int ipos, double** inverta);
int invertmatfour(int pos, double** a, int ipos, double** inverta);


//int invertfunction(int loopidloc, int* loopid, int blocksize, int* blocks, int* blockspos, double** mata, double*** mirrormat, double** inverta);
//int schurfunction(int loopidloc, int* loopid, int blocksize, int* blocks, int* blockspos, double** mata, double*** mirrormat, double** inverta);

int invertblocks(int order, double** mata, double** inverta);

struct block
{
	int blkm;
	int blkn;
	double **blkelement;
};

struct blockmat
{
	int brows;
	int bcols;
	struct block **blk;
};

struct mirrorstruct
{
	int morder; // order of the mirror mirror0 is 1, mirror1 is 2, mirror2 is 4, ... //order tells number of blocks
	int mblocksize;//total number of blocks in the mirror
	//int *mblockorder; //order of the blocks in the mirror
	struct blockmat *sqrblocks; //number of blocks are given by mblocksize
	struct blockmat *lpartnerblocks; // for each sqr block, one left partner block and
	struct blockmat *rpartnerblocks; // one right partner block
};

//typedef double **block;
//typedef block **blockmat;
//We are setting the precision here
//double PRECP=pow(10.0,PREC);
///double PRECP=1e100;
/*
int printer(int order, double** mata, double** inverta, int sizeofmirror, struct mirrorstruct *mptr)
{
	int i,j,k,l,m,n;
	printf("\n--------------------------------------\n");
	printf("printing the matrix A\n");
	for(i=0;i<order;i++)
	{
		for(j=0;j<order;j++) printf("%lf\t",mata[i][j]);	
		printf("\n");
	}
	printf("printing the inverted matrix invertA\n");
	for(i=0;i<order;i++)
	{
		for(j=0;j<order;j++) printf("%lf\t",inverta[i][j]);	
		printf("\n");
	}
	for(m=0;m<sizeofmirror;m++)
	{
		printf("printing the mirror matrix %d\n",m);
		printf("sqr blocks\n");
		for(i=0;i<mptr[m].mblocksize;i++)
		{
			for(j=0;j<mptr[m].sqrblocks[i].brows;j++)
			for(k=0;k<mptr[m].sqrblocks[i].bcols;k++)
			{
				printf("sqr blk of %d:\t%d %d : \n",i,j,k);
				for(l=0;l<mptr[m].sqrblocks[i].blk[j][k].blkm;l++)
				{
					for(n=0;n<mptr[m].sqrblocks[i].blk[j][k].blkn;n++)
						printf("%lf\t",mptr[m].sqrblocks[i].blk[j][k].blkelement[l][n]);	
					printf("\n");
				}
			}
		}
		printf("left and right partner blocks\n");
		for(i=0;(i<mptr[m].mblocksize/2);i++)
		{
			for(j=0;j<(mptr[m].rpartnerblocks[i].brows);j++)
			for(k=0;k<(mptr[m].rpartnerblocks[i].bcols);k++)
			{
				printf("right partner blk of %d:\t %d %d : \n",i,j,k);
				for(l=0;l<mptr[m].rpartnerblocks[i].blk[j][k].blkm;l++)
				{
					for(n=0;n<mptr[m].rpartnerblocks[i].blk[j][k].blkn;n++)
						printf("%lf\t",mptr[m].rpartnerblocks[i].blk[j][k].blkelement[l][n]);	
					printf("\n");
				}
			}
			for(j=0;j<(mptr[m].lpartnerblocks[i].brows);j++)
			for(k=0;k<(mptr[m].lpartnerblocks[i].bcols);k++)
			{
				printf("left partner blk of %d:\t %d %d : \n",i,j,k);
				for(l=0;l<mptr[m].lpartnerblocks[i].blk[j][k].blkm;l++)
				{
					for(n=0;n<mptr[m].lpartnerblocks[i].blk[j][k].blkn;n++)
						printf("%lf\t",mptr[m].lpartnerblocks[i].blk[j][k].blkelement[l][n]);	
					printf("\n");
				}
			}
		}
	}
	return 1;
}
*/
int invertmat(int n, double** mata, double** inverta)
{
	int invertstatus=0;
	int i,j;
	//for(n=5;n<=125;n=n+5)
	//{
	if(n<=0) return 0;
	if(n<=4) 
	{	
		//for(i=0;i<n;i++) for(j=0;j<n;j++) inverta[i][j]=mata[i][j];
		invertstatus = invertcases(n, 0, mata, 0, inverta);
	}
	else invertstatus = invertblocks(n, mata, inverta);
	//}
	return invertstatus;
}

int invertcases(int n, int apos, double** a, int invapos, double** inverta)
{
	int invertstatus;
	int i,j,k;
	
	if(n<=0) return 0;
	switch(n)
	{
		case 1:
			invertstatus=invertmatone(apos,a,invapos,inverta);
			break;
		case 2:
			invertstatus=invertmattwo(apos,a,invapos,inverta);
			break;
		case 3:
			invertstatus=invertmatthree(apos,a,invapos,inverta);
			break;
		case 4:
			invertstatus=invertmatfour(apos,a,invapos,inverta);
			break;
		//default:
			//invertstatus=invertblocks(n, mata, inverta);
			//break;
	}

	return invertstatus;
}

int invertmatone(int pos, double** a, int ipos, double** inverta)
{
	//int n=1;
	double modmata;

	//double PRECP=pow(10.0,PREC);

	modmata=a[pos+0][pos+0];
	
	if(modmata==0) return 0;
	
	inverta[ipos+0][ipos+0]=(1/modmata);
	
	return 1;
}

int invertmattwo(int pos, double** a, int ipos, double** inverta)
{
	//int n=2;

	double modmata;
	//double a11=a[pos+0][pos+0], a12=a[pos+0][pos+1];
	//double a21=a[pos+1][pos+0], a22=a[pos+1][pos+1];
	
	//double PRECP=pow(10.0,PREC);

	//modmata=floor(PRECP*(-(a12*a21) + a11*a22))/PRECP;
	modmata=(-(a[pos+0][pos+1]*a[pos+1][pos+0]) + a[pos+0][pos+0]*a[pos+1][pos+1]);
	
	if(modmata==0) 
	{ 
		printf("\n********** Unable to invert in blockpos=%d with order=%d ****************\n", pos, 2); 
		return 0;
	}
	
	//inverta[pos+0][pos+0]=floor(PRECP*(a22/modmata))/PRECP;
	//inverta[pos+0][pos+1]=floor(PRECP*(-(a12/modmata)))/PRECP;
	//inverta[pos+1][pos+0]=floor(PRECP*(-(a21/modmata)))/PRECP;
	//inverta[pos+1][pos+1]=floor(PRECP*(a11/modmata))/PRECP;
	
	inverta[ipos+0][ipos+0]=(a[pos+1][pos+1]/modmata);
	inverta[ipos+0][ipos+1]=(-(a[pos+0][pos+1]/modmata));
	inverta[ipos+1][ipos+0]=(-(a[pos+1][pos+0]/modmata));
	inverta[ipos+1][ipos+1]=(a[pos+0][pos+0]/modmata);
	
	return 1;
}

int invertmatthree(int pos, double** a, int ipos, double** inverta)
{
	//int n=3;

	double modmata;
	//double a11=a[pos+0][pos+0], a12=a[pos+0][pos+1], a13=a[pos+0][pos+2];
	//double a21=a[pos+1][pos+0], a22=a[pos+1][pos+1], a23=a[pos+1][pos+2];
	//double a31=a[pos+2][pos+0], a32=a[pos+2][pos+1], a33=a[pos+2][pos+2];

	//double PRECP=pow(10.0,PREC);

//	printf("\ninverting matthree is called position = %d \n", pos);	
	//modmata=floor(PRECP*(-(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33))/PRECP;
	modmata=(-(a[pos+0][pos+2]*a[pos+1][pos+1]*a[pos+2][pos+0]) + a[pos+0][pos+1]*a[pos+1][pos+2]*a[pos+2][pos+0] + a[pos+0][pos+2]*a[pos+1][pos+0]*a[pos+2][pos+1] - a[pos+0][pos+0]*a[pos+1][pos+2]*a[pos+2][pos+1] - a[pos+0][pos+1]*a[pos+1][pos+0]*a[pos+2][pos+2] + a[pos+0][pos+0]*a[pos+1][pos+1]*a[pos+2][pos+2]);
	
	if(modmata==0) 
	{ 
		printf("\n********** Unable to invert in blockpos=%d with order=%d ****************\n", pos, 3); 
		return 0;
	}
	
	//inverta[pos+0][pos+0]=floor(PRECP*((-(a23*a32) + a22*a33)/modmata))/PRECP;
	
	//inverta[pos+0][pos+1]=floor(PRECP*((a13*a32 - a12*a33)/modmata))/PRECP;
	
	//inverta[pos+0][pos+2]=floor(PRECP*((-(a13*a22) + a12*a23)/modmata))/PRECP;
	
	//inverta[pos+1][pos+0]=floor(PRECP*((a23*a31 - a21*a33)/modmata))/PRECP;
	
	//inverta[pos+1][pos+1]=floor(PRECP*((-(a13*a31) + a11*a33)/modmata))/PRECP;
	
	//inverta[pos+1][pos+2]=floor(PRECP*((a13*a21 - a11*a23)/modmata))/PRECP;
	
	//inverta[pos+2][pos+0]=floor(PRECP*((-(a22*a31) + a21*a32)/modmata))/PRECP;
	
	//inverta[pos+2][pos+1]=floor(PRECP*((a12*a31 - a11*a32)/modmata))/PRECP;
	
	//inverta[pos+2][pos+2]=floor(PRECP*((-(a12*a21) + a11*a22)/modmata))/PRECP;

	inverta[ipos+0][ipos+0]=((-(a[pos+1][pos+2]*a[pos+2][pos+1]) + a[pos+1][pos+1]*a[pos+2][pos+2])/modmata);
	
	inverta[ipos+0][ipos+1]=((a[pos+0][pos+2]*a[pos+2][pos+1] - a[pos+0][pos+1]*a[pos+2][pos+2])/modmata);
	
	inverta[ipos+0][ipos+2]=((-(a[pos+0][pos+2]*a[pos+1][pos+1]) + a[pos+0][pos+1]*a[pos+1][pos+2])/modmata);
	
	inverta[ipos+1][ipos+0]=((a[pos+1][pos+2]*a[pos+2][pos+0] - a[pos+1][pos+0]*a[pos+2][pos+2])/modmata);
	
	inverta[ipos+1][ipos+1]=((-(a[pos+0][pos+2]*a[pos+2][pos+0]) + a[pos+0][pos+0]*a[pos+2][pos+2])/modmata);
	
	inverta[ipos+1][ipos+2]=((a[pos+0][pos+2]*a[pos+1][pos+0] - a[pos+0][pos+0]*a[pos+1][pos+2])/modmata);
	
	inverta[ipos+2][ipos+0]=((-(a[pos+1][pos+1]*a[pos+2][pos+0]) + a[pos+1][pos+0]*a[pos+2][pos+1])/modmata);
	
	inverta[ipos+2][ipos+1]=((a[pos+0][pos+1]*a[pos+2][pos+0] - a[pos+0][pos+0]*a[pos+2][pos+1])/modmata);
	
	inverta[ipos+2][ipos+2]=((-(a[pos+0][pos+1]*a[pos+1][pos+0]) + a[pos+0][pos+0]*a[pos+1][pos+1])/modmata);

//	printf("\ninverted 3 by 3\n");
//	for(i=0;i<3;i++)
//	{
//		for(j=0;j<3;j++) printf("\t%lf",inverta[pos+i][pos+j]);
//		printf("\n");
//	}

	return 1;
}

int invertmatfour(int pos, double** a, int ipos, double** inverta)
{
	//int n=4;

	double modmata;
	//double a11=a[pos+0][pos+0], a12=a[pos+0][pos+1], a13=a[pos+0][pos+2], a14=a[pos+0][pos+3];
	//double a21=a[pos+1][pos+0], a22=a[pos+1][pos+1], a23=a[pos+1][pos+2], a24=a[pos+1][pos+3];
	//double a31=a[pos+2][pos+0], a32=a[pos+2][pos+1], a33=a[pos+2][pos+2], a34=a[pos+2][pos+3];
	//double a41=a[pos+3][pos+0], a42=a[pos+3][pos+1], a43=a[pos+3][pos+2], a44=a[pos+3][pos+3];
	
	//double PRECP=pow(10.0,PREC);
	
	//modmata=floor(PRECP*(a14*a23*a32*a41 - a13*a24*a32*a41 - a14*a22*a33*a41 + a12*a24*a33*a41 + a13*a22*a34*a41 - a12*a23*a34*a41 - a14*a23*a31*a42 + a13*a24*a31*a42 + a14*a21*a33*a42 - a11*a24*a33*a42 - a13*a21*a34*a42 + a11*a23*a34*a42 + a14*a22*a31*a43 - a12*a24*a31*a43 - a14*a21*a32*a43 + a11*a24*a32*a43 + a12*a21*a34*a43 - a11*a22*a34*a43 - a13*a22*a31*a44 + a12*a23*a31*a44 + a13*a21*a32*a44 - a11*a23*a32*a44 - a12*a21*a33*a44 + a11*a22*a33*a44))/PRECP;
	modmata=(a[pos+0][pos+3]*a[pos+1][pos+2]*a[pos+2][pos+1]*a[pos+3][pos+0] - a[pos+0][pos+2]*a[pos+1][pos+3]*a[pos+2][pos+1]*a[pos+3][pos+0] - a[pos+0][pos+3]*a[pos+1][pos+1]*a[pos+2][pos+2]*a[pos+3][pos+0] + a[pos+0][pos+1]*a[pos+1][pos+3]*a[pos+2][pos+2]*a[pos+3][pos+0] + a[pos+0][pos+2]*a[pos+1][pos+1]*a[pos+2][pos+3]*a[pos+3][pos+0] - a[pos+0][pos+1]*a[pos+1][pos+2]*a[pos+2][pos+3]*a[pos+3][pos+0] - a[pos+0][pos+3]*a[pos+1][pos+2]*a[pos+2][pos+0]*a[pos+3][pos+1] + a[pos+0][pos+2]*a[pos+1][pos+3]*a[pos+2][pos+0]*a[pos+3][pos+1] + a[pos+0][pos+3]*a[pos+1][pos+0]*a[pos+2][pos+2]*a[pos+3][pos+1] - a[pos+0][pos+0]*a[pos+1][pos+3]*a[pos+2][pos+2]*a[pos+3][pos+1] - a[pos+0][pos+2]*a[pos+1][pos+0]*a[pos+2][pos+3]*a[pos+3][pos+1] + a[pos+0][pos+0]*a[pos+1][pos+2]*a[pos+2][pos+3]*a[pos+3][pos+1] + a[pos+0][pos+3]*a[pos+1][pos+1]*a[pos+2][pos+0]*a[pos+3][pos+2] - a[pos+0][pos+1]*a[pos+1][pos+3]*a[pos+2][pos+0]*a[pos+3][pos+2] - a[pos+0][pos+3]*a[pos+1][pos+0]*a[pos+2][pos+1]*a[pos+3][pos+2] + a[pos+0][pos+0]*a[pos+1][pos+3]*a[pos+2][pos+1]*a[pos+3][pos+2] + a[pos+0][pos+1]*a[pos+1][pos+0]*a[pos+2][pos+3]*a[pos+3][pos+2] - a[pos+0][pos+0]*a[pos+1][pos+1]*a[pos+2][pos+3]*a[pos+3][pos+2] - a[pos+0][pos+2]*a[pos+1][pos+1]*a[pos+2][pos+0]*a[pos+3][pos+3] + a[pos+0][pos+1]*a[pos+1][pos+2]*a[pos+2][pos+0]*a[pos+3][pos+3] + a[pos+0][pos+2]*a[pos+1][pos+0]*a[pos+2][pos+1]*a[pos+3][pos+3] - a[pos+0][pos+0]*a[pos+1][pos+2]*a[pos+2][pos+1]*a[pos+3][pos+3] - a[pos+0][pos+1]*a[pos+1][pos+0]*a[pos+2][pos+2]*a[pos+3][pos+3] + a[pos+0][pos+0]*a[pos+1][pos+1]*a[pos+2][pos+2]*a[pos+3][pos+3]);
	
	if(modmata==0) 
	{ 
		printf("\n********** Unable to invert in blockpos=%d with order=%d ****************\n", pos, 4); 
		return 0;
	}
	

	//inverta[pos+0][pos+0]=floor(PRECP*((-(a24*a33*a42) + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 + a22*a33*a44)/modmata))/PRECP;
	inverta[ipos+0][ipos+0]=((-(a[pos+1][pos+3]*a[pos+2][pos+2]*a[pos+3][pos+1]) + a[pos+1][pos+2]*a[pos+2][pos+3]*a[pos+3][pos+1] + a[pos+1][pos+3]*a[pos+2][pos+1]*a[pos+3][pos+2] - a[pos+1][pos+1]*a[pos+2][pos+3]*a[pos+3][pos+2] - a[pos+1][pos+2]*a[pos+2][pos+1]*a[pos+3][pos+3] + a[pos+1][pos+1]*a[pos+2][pos+2]*a[pos+3][pos+3])/modmata);
	
	//inverta[pos+0][pos+1]=floor(PRECP*((a14*a33*a42 - a13*a34*a42 - a14*a32*a43 + a12*a34*a43 + a13*a32*a44 - a12*a33*a44)/modmata))/PRECP;
	inverta[ipos+0][ipos+1]=((a[pos+0][pos+3]*a[pos+2][pos+2]*a[pos+3][pos+1] - a[pos+0][pos+2]*a[pos+2][pos+3]*a[pos+3][pos+1] - a[pos+0][pos+3]*a[pos+2][pos+1]*a[pos+3][pos+2] + a[pos+0][pos+1]*a[pos+2][pos+3]*a[pos+3][pos+2] + a[pos+0][pos+2]*a[pos+2][pos+1]*a[pos+3][pos+3] - a[pos+0][pos+1]*a[pos+2][pos+2]*a[pos+3][pos+3])/modmata);
	
	//inverta[pos+0][pos+2]=floor(PRECP*((-(a14*a23*a42) + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 + a12*a23*a44)/modmata))/PRECP;
	inverta[ipos+0][ipos+2]=((-(a[pos+0][pos+3]*a[pos+1][pos+2]*a[pos+3][pos+1]) + a[pos+0][pos+2]*a[pos+1][pos+3]*a[pos+3][pos+1] + a[pos+0][pos+3]*a[pos+1][pos+1]*a[pos+3][pos+2] - a[pos+0][pos+1]*a[pos+1][pos+3]*a[pos+3][pos+2] - a[pos+0][pos+2]*a[pos+1][pos+1]*a[pos+3][pos+3] + a[pos+0][pos+1]*a[pos+1][pos+2]*a[pos+3][pos+3])/modmata);
	
	//inverta[pos+0][pos+3]=floor(PRECP*((a14*a23*a32 - a13*a24*a32 - a14*a22*a33 + a12*a24*a33 + a13*a22*a34 - a12*a23*a34)/modmata))/PRECP;
	inverta[ipos+0][ipos+3]=((a[pos+0][pos+3]*a[pos+1][pos+2]*a[pos+2][pos+1] - a[pos+0][pos+2]*a[pos+1][pos+3]*a[pos+2][pos+1] - a[pos+0][pos+3]*a[pos+1][pos+1]*a[pos+2][pos+2] + a[pos+0][pos+1]*a[pos+1][pos+3]*a[pos+2][pos+2] + a[pos+0][pos+2]*a[pos+1][pos+1]*a[pos+2][pos+3] - a[pos+0][pos+1]*a[pos+1][pos+2]*a[pos+2][pos+3])/modmata);
	
	//inverta[pos+1][pos+0]=floor(PRECP*( (a24*a33*a41 - a23*a34*a41 - a24*a31*a43 + a21*a34*a43 + a23*a31*a44 - a21*a33*a44)/modmata))/PRECP;
	inverta[ipos+1][ipos+0]=( (a[pos+1][pos+3]*a[pos+2][pos+2]*a[pos+3][pos+0] - a[pos+1][pos+2]*a[pos+2][pos+3]*a[pos+3][pos+0] - a[pos+1][pos+3]*a[pos+2][pos+0]*a[pos+3][pos+2] + a[pos+1][pos+0]*a[pos+2][pos+3]*a[pos+3][pos+2] + a[pos+1][pos+2]*a[pos+2][pos+0]*a[pos+3][pos+3] - a[pos+1][pos+0]*a[pos+2][pos+2]*a[pos+3][pos+3])/modmata);
	
	//inverta[pos+1][pos+1]=floor(PRECP*((-(a14*a33*a41) + a13*a34*a41 + a14*a31*a43 - a11*a34*a43 - a13*a31*a44 + a11*a33*a44)/modmata))/PRECP;
	inverta[ipos+1][ipos+1]=((-(a[pos+0][pos+3]*a[pos+2][pos+2]*a[pos+3][pos+0]) + a[pos+0][pos+2]*a[pos+2][pos+3]*a[pos+3][pos+0] + a[pos+0][pos+3]*a[pos+2][pos+0]*a[pos+3][pos+2] - a[pos+0][pos+0]*a[pos+2][pos+3]*a[pos+3][pos+2] - a[pos+0][pos+2]*a[pos+2][pos+0]*a[pos+3][pos+3] + a[pos+0][pos+0]*a[pos+2][pos+2]*a[pos+3][pos+3])/modmata);
	
	//inverta[pos+1][pos+2]=floor(PRECP*((a14*a23*a41 - a13*a24*a41 - a14*a21*a43 + a11*a24*a43 + a13*a21*a44 - a11*a23*a44)/modmata))/PRECP;
	inverta[ipos+1][ipos+2]=((a[pos+0][pos+3]*a[pos+1][pos+2]*a[pos+3][pos+0] - a[pos+0][pos+2]*a[pos+1][pos+3]*a[pos+3][pos+0] - a[pos+0][pos+3]*a[pos+1][pos+0]*a[pos+3][pos+2] + a[pos+0][pos+0]*a[pos+1][pos+3]*a[pos+3][pos+2] + a[pos+0][pos+2]*a[pos+1][pos+0]*a[pos+3][pos+3] - a[pos+0][pos+0]*a[pos+1][pos+2]*a[pos+3][pos+3])/modmata);
	
	//inverta[pos+1][pos+3]=floor(PRECP*((-(a14*a23*a31) + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 + a11*a23*a34)/modmata))/PRECP;
	inverta[ipos+1][ipos+3]=((-(a[pos+0][pos+3]*a[pos+1][pos+2]*a[pos+2][pos+0]) + a[pos+0][pos+2]*a[pos+1][pos+3]*a[pos+2][pos+0] + a[pos+0][pos+3]*a[pos+1][pos+0]*a[pos+2][pos+2] - a[pos+0][pos+0]*a[pos+1][pos+3]*a[pos+2][pos+2] - a[pos+0][pos+2]*a[pos+1][pos+0]*a[pos+2][pos+3] + a[pos+0][pos+0]*a[pos+1][pos+2]*a[pos+2][pos+3])/modmata);
	
	//inverta[pos+2][pos+0]=floor(PRECP*((-(a24*a32*a41) + a22*a34*a41 + a24*a31*a42 - a21*a34*a42 - a22*a31*a44 + a21*a32*a44)/modmata))/PRECP;
	inverta[ipos+2][ipos+0]=((-(a[pos+1][pos+3]*a[pos+2][pos+1]*a[pos+3][pos+0]) + a[pos+1][pos+1]*a[pos+2][pos+3]*a[pos+3][pos+0] + a[pos+1][pos+3]*a[pos+2][pos+0]*a[pos+3][pos+1] - a[pos+1][pos+0]*a[pos+2][pos+3]*a[pos+3][pos+1] - a[pos+1][pos+1]*a[pos+2][pos+0]*a[pos+3][pos+3] + a[pos+1][pos+0]*a[pos+2][pos+1]*a[pos+3][pos+3])/modmata);
	
	//inverta[pos+2][pos+1]=floor(PRECP*((a14*a32*a41 - a12*a34*a41 - a14*a31*a42 + a11*a34*a42 + a12*a31*a44 - a11*a32*a44)/modmata))/PRECP;
	inverta[ipos+2][ipos+1]=((a[pos+0][pos+3]*a[pos+2][pos+1]*a[pos+3][pos+0] - a[pos+0][pos+1]*a[pos+2][pos+3]*a[pos+3][pos+0] - a[pos+0][pos+3]*a[pos+2][pos+0]*a[pos+3][pos+1] + a[pos+0][pos+0]*a[pos+2][pos+3]*a[pos+3][pos+1] + a[pos+0][pos+1]*a[pos+2][pos+0]*a[pos+3][pos+3] - a[pos+0][pos+0]*a[pos+2][pos+1]*a[pos+3][pos+3])/modmata);
	
	//inverta[pos+2][pos+2]=floor(PRECP*((-(a14*a22*a41) + a12*a24*a41 + a14*a21*a42 - a11*a24*a42 - a12*a21*a44 + a11*a22*a44)/modmata))/PRECP;
	inverta[ipos+2][ipos+2]=((-(a[pos+0][pos+3]*a[pos+1][pos+1]*a[pos+3][pos+0]) + a[pos+0][pos+1]*a[pos+1][pos+3]*a[pos+3][pos+0] + a[pos+0][pos+3]*a[pos+1][pos+0]*a[pos+3][pos+1] - a[pos+0][pos+0]*a[pos+1][pos+3]*a[pos+3][pos+1] - a[pos+0][pos+1]*a[pos+1][pos+0]*a[pos+3][pos+3] + a[pos+0][pos+0]*a[pos+1][pos+1]*a[pos+3][pos+3])/modmata);
	
	//inverta[pos+2][pos+3]=floor(PRECP*((a14*a22*a31 - a12*a24*a31 - a14*a21*a32 + a11*a24*a32 + a12*a21*a34 - a11*a22*a34)/modmata))/PRECP;
	inverta[ipos+2][ipos+3]=((a[pos+0][pos+3]*a[pos+1][pos+1]*a[pos+2][pos+0] - a[pos+0][pos+1]*a[pos+1][pos+3]*a[pos+2][pos+0] - a[pos+0][pos+3]*a[pos+1][pos+0]*a[pos+2][pos+1] + a[pos+0][pos+0]*a[pos+1][pos+3]*a[pos+2][pos+1] + a[pos+0][pos+1]*a[pos+1][pos+0]*a[pos+2][pos+3] - a[pos+0][pos+0]*a[pos+1][pos+1]*a[pos+2][pos+3])/modmata);
	
	//inverta[pos+3][pos+0]=floor(PRECP*((a23*a32*a41 - a22*a33*a41 - a23*a31*a42 + a21*a33*a42 + a22*a31*a43 - a21*a32*a43)/modmata))/PRECP;
	inverta[ipos+3][ipos+0]=((a[pos+1][pos+2]*a[pos+2][pos+1]*a[pos+3][pos+0] - a[pos+1][pos+1]*a[pos+2][pos+2]*a[pos+3][pos+0] - a[pos+1][pos+2]*a[pos+2][pos+0]*a[pos+3][pos+1] + a[pos+1][pos+0]*a[pos+2][pos+2]*a[pos+3][pos+1] + a[pos+1][pos+1]*a[pos+2][pos+0]*a[pos+3][pos+2] - a[pos+1][pos+0]*a[pos+2][pos+1]*a[pos+3][pos+2])/modmata);
	
	//inverta[pos+3][pos+1]=floor(PRECP*((-(a13*a32*a41) + a12*a33*a41 + a13*a31*a42 - a11*a33*a42 - a12*a31*a43 + a11*a32*a43)/modmata))/PRECP;
	inverta[ipos+3][ipos+1]=((-(a[pos+0][pos+2]*a[pos+2][pos+1]*a[pos+3][pos+0]) + a[pos+0][pos+1]*a[pos+2][pos+2]*a[pos+3][pos+0] + a[pos+0][pos+2]*a[pos+2][pos+0]*a[pos+3][pos+1] - a[pos+0][pos+0]*a[pos+2][pos+2]*a[pos+3][pos+1] - a[pos+0][pos+1]*a[pos+2][pos+0]*a[pos+3][pos+2] + a[pos+0][pos+0]*a[pos+2][pos+1]*a[pos+3][pos+2])/modmata);
	
	//inverta[pos+3][pos+2]=floor(PRECP*((a13*a22*a41 - a12*a23*a41 - a13*a21*a42 + a11*a23*a42 + a12*a21*a43 - a11*a22*a43)/modmata))/PRECP;
	inverta[ipos+3][ipos+2]=((a[pos+0][pos+2]*a[pos+1][pos+1]*a[pos+3][pos+0] - a[pos+0][pos+1]*a[pos+1][pos+2]*a[pos+3][pos+0] - a[pos+0][pos+2]*a[pos+1][pos+0]*a[pos+3][pos+1] + a[pos+0][pos+0]*a[pos+1][pos+2]*a[pos+3][pos+1] + a[pos+0][pos+1]*a[pos+1][pos+0]*a[pos+3][pos+2] - a[pos+0][pos+0]*a[pos+1][pos+1]*a[pos+3][pos+2])/modmata);
	
	//inverta[pos+3][pos+3]=floor(PRECP*((-(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33)/modmata))/PRECP;
	inverta[ipos+3][ipos+3]=((-(a[pos+0][pos+2]*a[pos+1][pos+1]*a[pos+2][pos+0]) + a[pos+0][pos+1]*a[pos+1][pos+2]*a[pos+2][pos+0] + a[pos+0][pos+2]*a[pos+1][pos+0]*a[pos+2][pos+1] - a[pos+0][pos+0]*a[pos+1][pos+2]*a[pos+2][pos+1] - a[pos+0][pos+1]*a[pos+1][pos+0]*a[pos+2][pos+2] + a[pos+0][pos+0]*a[pos+1][pos+1]*a[pos+2][pos+2])/modmata);


	return 1;
}



int invertblocks(int order, double** mata, double** inverta)
{

	int i,j,k,l,m,n;
	int mcid, msid, ud, udmc, msiditr;
	int nloop, threadid;
	int **loopid, loopidsize, noofloops;
	int *blocks, *blockspos, blocksize;

	int ii, jj, kk, iblk, itr, blkchoice, itrn;

	//struct mirrorstruct mirror[((int)log2(pow(2,(int)log2(order)-1)))];
	struct mirrorstruct *mirror;
	int mirrorsize;

	int *morder, *mblocksize;

	double temp, temp1, temp2, temp3, temp4;
	double wtime;

	blocksize = (int) log2(order);
	blocksize = pow(2,(int)log2(order)-1);
	blocks=(int *) malloc(blocksize*sizeof(int));
	blockspos=(int *) malloc(blocksize*sizeof(int));
	
	////printf("\n order = %d, blocksize = %d", order, blocksize);

        for(i=0;i<blocksize;i++)
        {       
                n=(int)log2(order)-1;
                if(order<((int)pow(2,n)*3)) blocks[i]=(i<((int)pow(2,n)- order%((int)pow(2,n))))?2:3;
                else blocks[i]=(i<((int)pow(2,n)- order%((int)pow(2,n))))?3:4;
        }

        for(blockspos[0]=0,i=1;i<blocksize; i++) blockspos[i]=blockspos[i-1]+blocks[i-1];

        ////printf("\nPrinting the blocks array and blockspos array of size %d\n",blocksize);
        /////for(i=0;i<blocksize;i++)
        ////printf("%d\t%d\n",blocks[i],blockspos[i]);


	loopidsize = (int)log2(blocksize)+1;
	noofloops=blocksize*2;
	//loopid = (int *) calloc(loopidsize,sizeof(int));
	//loopid as two dim array
	loopid=(int **)calloc(noofloops,sizeof(int *));
	for(i=0;i<noofloops;i++) loopid[i]=(int *)calloc(loopidsize,sizeof(int));

	////printf("\n number of loops %d\n",noofloops);
	for(i=1;i<noofloops;i++)
	{
		//for(k=0;k<8;k++) printf("\tk is %d, %d, %d :",k,((int)pow(2*1.0,k*1.0)), (int)(i/(pow(2*1.0,k*1.0))));
		for(k=0;((int)(i/pow(2*1.0,k*1.0)))>0;k++);
		loopid[i][0]=k;
		
		//Calculation for loopid[i]:
		for(j=1;j<loopid[i][0];j++)
		{
			l=0;
			if(loopid[i][j-1]>0) l=(int)(i%((int)(pow(2*1.0,loopid[i][j-1]-1.0))));	
			for(k=0;((int)(l/pow(2*1.0,k*1.0)))>0;k++);
			loopid[i][j]=k;
		}
	}

	//Decleration of mirrors
	
	mirror=(struct mirrorstruct *)malloc(((int)log2(pow(2,(int)log2(order)-1))) * sizeof(struct mirrorstruct));
	//printf("Mirrors count = %d  \n",sizeof(mirror)/sizeof(mirror[0]));
	//printf("Mirrors count = %d  \n",((int)log2(pow(2,(int)log2(order)-1))));
	mirrorsize=(int)log2(pow(2,(int)log2(order)-1));
	//for(i=0;i<(sizeof(mirror)/sizeof(mirror[0]));i++)
        morder=(int *)calloc(mirrorsize,sizeof(int));
        mblocksize=(int *)calloc(mirrorsize,sizeof(int));

	for(i=0;i<mirrorsize;i++)
	//for(i=0;i<(sizeof(mirror)/sizeof(mirror[0]));i++)
	{
                morder[i]=(int)pow(2,i);
                mblocksize[i]=(int)(blocksize/morder[i]);//blocksize interms of number of blocks comparing blocks array

		mirror[i].morder=(int)pow(2,i);
		mirror[i].mblocksize=(int)(blocksize/mirror[i].morder);//blocksize interms of number of blocks comparing blocks array
		//mirror[i].mblocks=(int *)malloc(((int)blocksize/mirror[i].morder)*sizeof(int));
		//mirror[i].mblockspos=(int *)malloc(((int)blocksize/mirror[i].morder)*sizeof(int));
		
		mirror[i].sqrblocks=(struct blockmat *)calloc(mirror[i].mblocksize , sizeof(struct block));
		//mirror[i].sqrblocks=(struct block *)malloc(mirror[i].mblocksize * sizeof(struct block));
		mirror[i].lpartnerblocks=(struct blockmat *)calloc((int)(mirror[i].mblocksize/2) , sizeof(struct block));
		mirror[i].rpartnerblocks=(struct blockmat *)calloc((int)(mirror[i].mblocksize/2) , sizeof(struct block));
		for(j=0; j<mirror[i].mblocksize; j++)
		{
			mirror[i].sqrblocks[j].brows=mirror[i].morder;
			mirror[i].sqrblocks[j].bcols=mirror[i].morder;
			mirror[i].sqrblocks[j].blk=(struct block **)calloc(mirror[i].sqrblocks[j].brows,sizeof(struct block *));
			for(k=0; k< mirror[i].sqrblocks[j].brows;k++)
			{
				mirror[i].sqrblocks[j].blk[k]=(struct block *)calloc(mirror[i].sqrblocks[j].bcols,sizeof(struct block));
				for(l=0;l<mirror[i].sqrblocks[j].bcols;l++)
				{
						mirror[i].sqrblocks[j].blk[k][l].blkm=blocks[k+(j*mirror[i].morder)];
						mirror[i].sqrblocks[j].blk[k][l].blkn=blocks[l+(j*mirror[i].morder)];
						mirror[i].sqrblocks[j].blk[k][l].blkelement=(double **)calloc(mirror[i].sqrblocks[j].blk[k][l].blkm,sizeof(double *));
						for(m=0;m<(mirror[i].sqrblocks[j].blk[k][l].blkm);m++)
							mirror[i].sqrblocks[j].blk[k][l].blkelement[m]=(double *)calloc(mirror[i].sqrblocks[j].blk[k][l].blkn,sizeof(double));
				}
			}
			if((j%2)==0)
			{
				mirror[i].rpartnerblocks[(int)j/2].brows=mirror[i].morder;
				mirror[i].rpartnerblocks[(int)j/2].bcols=mirror[i].morder;
				mirror[i].rpartnerblocks[(int)j/2].blk=(struct block **)calloc(mirror[i].rpartnerblocks[(int)j/2].brows,sizeof(struct block *));
				for(k=0; k< mirror[i].rpartnerblocks[(int)j/2].brows;k++)
				{
					mirror[i].rpartnerblocks[(int)j/2].blk[k]=(struct block *)calloc(mirror[i].rpartnerblocks[(int)j/2].bcols,sizeof(struct block));
					for(l=0;l<mirror[i].rpartnerblocks[(int)j/2].bcols;l++)
					{
						mirror[i].rpartnerblocks[(int)j/2].blk[k][l].blkm=blocks[k+(j*mirror[i].morder)];
						mirror[i].rpartnerblocks[(int)j/2].blk[k][l].blkn=blocks[l+((j+1)*mirror[i].morder)];
						mirror[i].rpartnerblocks[(int)j/2].blk[k][l].blkelement=(double **)calloc(mirror[i].rpartnerblocks[(int)j/2].blk[k][l].blkm,sizeof(double *));
						for(m=0;m<(mirror[i].rpartnerblocks[(int)j/2].blk[k][l].blkm);m++)
							mirror[i].rpartnerblocks[(int)j/2].blk[k][l].blkelement[m]=(double *)calloc(mirror[i].rpartnerblocks[(int)j/2].blk[k][l].blkn,sizeof(double));
					}
				}
			}
			else
			{
				mirror[i].lpartnerblocks[(int)j/2].brows=mirror[i].morder;
				mirror[i].lpartnerblocks[(int)j/2].bcols=mirror[i].morder;
				mirror[i].lpartnerblocks[(int)j/2].blk=(struct block **)calloc(mirror[i].lpartnerblocks[(int)j/2].brows,sizeof(struct block *));
				for(k=0; k< mirror[i].lpartnerblocks[(int)j/2].brows;k++)
				{
					mirror[i].lpartnerblocks[(int)j/2].blk[k]=(struct block *)calloc(mirror[i].lpartnerblocks[(int)j/2].bcols,sizeof(struct block));
					for(l=0;l<mirror[i].lpartnerblocks[(int)j/2].bcols;l++)
					{
						mirror[i].lpartnerblocks[(int)j/2].blk[k][l].blkm=blocks[k+(j*mirror[i].morder)];
						mirror[i].lpartnerblocks[(int)j/2].blk[k][l].blkn=blocks[l+((j-1)*mirror[i].morder)];
						mirror[i].lpartnerblocks[(int)j/2].blk[k][l].blkelement=(double **)calloc(mirror[i].lpartnerblocks[(int)j/2].blk[k][l].blkm,sizeof(double *));
						for(m=0;m<(mirror[i].lpartnerblocks[(int)j/2].blk[k][l].blkm);m++)
							mirror[i].lpartnerblocks[(int)j/2].blk[k][l].blkelement[m]=(double *)calloc(mirror[i].lpartnerblocks[(int)j/2].blk[k][l].blkn,sizeof(double));
					}
				}
			}
		}
		////printf("\n Mirror %d creation : morder = %d , mblocksize = %d",i, mirror[i].morder, mirror[i].mblocksize);
	}

	//int printer(int order, double** mata, double** inverta, int sizeofmirror, struct mirrorstruct *mptr)
	//printer(order, mata, inverta, (int)(sizeof(mirror)/sizeof(mirror[0])), mirror);
	//printf("\ninversion with omp parallel, PREC = %d\n",PREC);	

	#pragma omp parallel private(wtime, i,j,nloop,blkchoice, msid, mcid, udmc) shared(order, noofloops, PRECP, blocks, blockspos, morder, mblocksize, mirror, mata, inverta)
	//for(i=1;i<noofloops;i++)
	{	i=1; //nloop=0;
		//wtime = omp_get_wtime();
		do
		{
			//printf("\n----------------------------loop %d---------------------------\n",i);
			//for(j=0;j<loopidsize;j++) printf("\t%d",loopid[i][j]); printf("\n--------------\n");
			//Here we perform the operation for given loopid:
			//Finding the loopid location for which the operation has to be done.
			for(j=loopidsize-1;(j>=0)&&(loopid[i][j]==0);j--);
			#pragma omp barrier	
			if(j==0)
			{
				if(loopid[i][0]==1)
				{
					//Inversion using A
					#pragma omp for private(k) schedule(dynamic)
					for(k=0;k<blocksize;k++)
					{
						invertcases(blocks[k], blockspos[k], mata, blockspos[k], inverta);
					}
				}
				else
				{
					//mirrorchoiceid is mata //mirrorstoreid is zero.
					msid = loopid[i][0]-2;
					//1. Calculation of -A^-1B and -D^-1C using A (for B and D) and storing at Mirror (id=loopid[0]-2) 
					#pragma omp for collapse(3) private(k,ii,jj,l,m) schedule(static,8)
					for(k=0;k<(mblocksize[msid]/2);k++)
					//for(k=0;k<(mirror[msid].mblocksize/2);k++)
					{
						for(ii=0; ii<morder[msid]; ii++)
						//for(ii=0; ii<mirror[msid].morder; ii++)
						{
							for(jj=0; jj<morder[msid]; jj++)
							//for(jj=0; jj<mirror[msid].morder; jj++)
							{
								//for(l=0; l<rpartblkm[msid][k][ii][jj]; l++)
								for(l=0; l<mirror[msid].rpartnerblocks[k].blk[ii][jj].blkm; l++)
									//for(m=0; m<rpartblkn[msid][k][ii][jj]; m++)
									for(m=0; m<mirror[msid].rpartnerblocks[k].blk[ii][jj].blkn; m++)
										//mrpartner[msid][k][ii][jj][l][m]=0.0;
										#pragma omp atomic write
										mirror[msid].rpartnerblocks[k].blk[ii][jj].blkelement[l][m]=0.0;
							}
						}
					}
					#pragma omp for collapse(3) private(k,ii,jj,l,m) schedule(static,8)
					for(k=0;k<(mblocksize[msid]/2);k++)
					//for(k=0;k<(mirror[msid].mblocksize/2);k++)
					{
						for(ii=0; ii<morder[msid]; ii++)
						//for(ii=0; ii<mirror[msid].morder; ii++)
						{
							for(jj=0; jj<morder[msid]; jj++)
							//for(jj=0; jj<mirror[msid].morder; jj++)
							{
								//for(l=0; l<lpartblkm[msid][k][ii][jj]; l++)
								for(l=0; l<mirror[msid].lpartnerblocks[k].blk[ii][jj].blkm; l++)
									//for(m=0; m<lpartblkn[msid][k][ii][jj]; m++)
									for(m=0; m<mirror[msid].lpartnerblocks[k].blk[ii][jj].blkn; m++)
										#pragma omp atomic write
										//mlpartner[msid][k][ii][jj][l][m]=0.0;
										mirror[msid].lpartnerblocks[k].blk[ii][jj].blkelement[l][m]=0.0;
							}
						}
					}
					#pragma omp for collapse(4) private(k,ii,jj,kk,l,m,n) reduction(+:temp) schedule(static)
					for(k=0;k<(mblocksize[msid]/2);k++)
					//for(k=0;k<(mirror[msid].mblocksize/2);k++)
					{
						for(ii=0; ii<morder[msid]; ii++)
						//for(ii=0; ii<mirror[msid].morder; ii++)
						{
							for(jj=0; jj<morder[msid]; jj++)
							//for(jj=0; jj<mirror[msid].morder; jj++)
							{
								for(kk=0; kk<morder[msid]; kk++)
								{
									//printf("\nmsid=%d, k=%d, ii=%d, jj=%d kk=%d\n",msid,k,ii,jj,kk);
									//for(l=0; l<rpartblkm[msid][k][jj][kk]; l++)
									for(l=0; l<mirror[msid].rpartnerblocks[k].blk[jj][kk].blkm; l++)
										//for(m=0; m<rpartblkn[msid][k][jj][kk]; m++)
										//	for(n=0;n<blocks[(k*2*morder[msid])+(jj+ii)%morder[msid]];n++)
										for(m=0; m<mirror[msid].rpartnerblocks[k].blk[jj][kk].blkn; m++)
									//for(l=0; l<mirror[msid].rpartnerblocks[k].blk[jj][kk].blkm; l++)
											{
											//for(n=0; n<mirror[msid].rpartnerblocks[k].blk[jj][kk].blkm; n++)
											//for(n=0;n<blocks[(k*2*mirror[msid].morder)+(jj+ii)%mirror[msid].morder];n++)
											//for(n=0, temp=mrpartner[msid][k][jj][kk][l][m];n<blocks[(k*2*morder[msid])+(jj+ii)%morder[msid]];n++)
											//temp=mrpartner[msid][k][jj][kk][l][m];
											temp=0.0;//mrpartner[msid][k][jj][kk][l][m];
											for(n=0;n<blocks[(k*2*morder[msid])+(jj+ii)%morder[msid]];n++)
											//mrpartner[msid][k][jj][kk][l][m]+=floor(PRECP*(inverta[blockspos[(k*2*morder[msid])+jj]+l][blockspos[(k*2*morder[msid])+(jj+ii)%morder[msid]]+n]*mata[blockspos[(k*2*morder[msid])+(jj+ii)%morder[msid]]+n][blockspos[((k*2+1)* morder[msid])+kk]+m]*(-1.0)))/PRECP;
											temp +=(inverta[blockspos[(k*2*morder[msid])+jj]+l][blockspos[(k*2*morder[msid])+(jj+ii)%morder[msid]]+n]*mata[blockspos[(k*2*morder[msid])+(jj+ii)%morder[msid]]+n][blockspos[((k*2+1)* morder[msid])+kk]+m]*(-1.0));
											//mirror[msid].rpartnerblocks[k].blk[jj][kk].blkelement[l][m]+=floor(PRECP*(inverta[blockspos[(k*2*morder[msid])+jj]+l][blockspos[(k*2*morder[msid])+(jj+ii)%morder[msid]]+n]*mata[blockspos[(k*2*morder[msid])+(jj+ii)%morder[msid]]+n][blockspos[((k*2+1)* morder[msid])+kk]+m]*(-1.0)))/PRECP;
										        #pragma omp atomic update
											mirror[msid].rpartnerblocks[k].blk[jj][kk].blkelement[l][m]+=temp;
											//printf("\n\t Test k %d ii %d\tjj %d\tkk %d\t l %d m %d temp %lf\n",k,ii,jj,kk,l,m,temp);
											}
											
									
								//mirror[msid].mirrormat[mirror[msid].mblockspos[k]+l][mirror[msid].mblockspos[k+1]+m]+=floor(PRECP*(-1.0*inverta[mirror[msid].mblockspos[k]+l][mirror[msid].mblockspos[k]+n]*mata[mirror[msid].mblockspos[k]+n][mirror[msid].mblockspos[k+1]+m]))/PRECP;
								//}
								//for(kk=0; kk<morder[msid]; kk++)
								//{
									//printf("\nmsid=%d, k=%d, ii=%d, jj=%d kk=%d\n",msid,k,ii,jj,kk);
								}
							}
						}
					}
					#pragma omp for collapse(4) private(k,ii,jj,kk,l,m,n) reduction(+:temp)
					for(k=0;k<(mblocksize[msid]/2);k++)
					//for(k=0;k<(mirror[msid].mblocksize/2);k++)
					{
						for(ii=0; ii<morder[msid]; ii++)
						//for(ii=0; ii<mirror[msid].morder; ii++)
						{
							for(jj=0; jj<morder[msid]; jj++)
							//for(jj=0; jj<mirror[msid].morder; jj++)
							{
								for(kk=0; kk<morder[msid]; kk++)
								{
									//printf("\nmsid=%d, k=%d, ii=%d, jj=%d kk=%d\n",msid,k,ii,jj,kk);
									//for(l=0; l<lpartblkm[msid][k][jj][kk]; l++)
									for(l=0; l<mirror[msid].lpartnerblocks[k].blk[jj][kk].blkm; l++)
										//for(m=0; m<lpartblkn[msid][k][jj][kk]; m++)
										for(m=0; m<mirror[msid].lpartnerblocks[k].blk[jj][kk].blkn; m++)
											{
											//temp=mlpartner[msid][k][jj][kk][l][m];
											temp=0.0;//=mlpartner[msid][k][jj][kk][l][m];
											//for(n=0; n<mirror[msid].lpartnerblocks[k].blk[jj][kk].blkm; n++)
											//for(n=0, temp=mlpartner[msid][k][jj][kk][l][m];n<blocks[((k*2+1)*morder[msid])+(jj+ii)%morder[msid]];n++)
											for(n=0;n<blocks[((k*2+1)*morder[msid])+(jj+ii)%morder[msid]];n++)
											//mlpartner[msid][k][jj][kk][l][m]+=floor(PRECP*(inverta[blockspos[((k*2+1)*morder[msid])+jj]+l][blockspos[((k*2+1)*morder[msid])+(jj+ii)%morder[msid]]+n]*mata[blockspos[((k*2+1)*morder[msid])+(jj+ii)%morder[msid]]+n][blockspos[((2*k)*morder[msid])+kk]+m]*(-1.0)))/PRECP;
											temp+=(inverta[blockspos[((k*2+1)*morder[msid])+jj]+l][blockspos[((k*2+1)*morder[msid])+(jj+ii)%morder[msid]]+n]*mata[blockspos[((k*2+1)*morder[msid])+(jj+ii)%morder[msid]]+n][blockspos[((2*k)*morder[msid])+kk]+m]*(-1.0));
										        #pragma omp atomic update
											mirror[msid].lpartnerblocks[k].blk[jj][kk].blkelement[l][m]+=temp;
											//mlpartner[msid][k][jj][kk][l][m]+=temp;
											//printf("\n\t Test k %d ii %d\tjj %d\tkk %d\t l %d m %d temp %lf\n",k,ii,jj,kk,l,m,temp);
											}
								//mirror[msid].mirrormat[mirror[msid].mblockspos[k]+l][mirror[msid].mblockspos[k+1]+m]+=floor(PRECP*(-1.0*inverta[mirror[msid].mblockspos[k]+l][mirror[msid].mblockspos[k]+n]*mata[mirror[msid].mblockspos[k]+n][mirror[msid].mblockspos[k+1]+m]))/PRECP;
								}
							}
						}
					}
					//#pragma omp barrier
					//2. Calculation of S_A and S_D using A at the location Mirror (id=loopid[0]-2) which is msid
					#pragma omp for collapse(3) private(k,ii,jj,l,m)
					for(k=0;k<(mblocksize[msid]);k++)
					//for(k=0;k<(mirror[msid].mblocksize);k++)
					{
						for(ii=0; ii<morder[msid]; ii++)
						//for(ii=0; ii<mirror[msid].morder; ii++)
						{
							for(jj=0; jj<morder[msid]; jj++)
							//for(jj=0; jj<mirror[msid].morder; jj++)
							{
								//for(l=0; l<sqrblkm[msid][k][ii][jj]; l++)
								for(l=0; l<mirror[msid].sqrblocks[k].blk[ii][jj].blkm; l++)
								//for(l=0; l<mirror[msid].sqrblocks[k].blk[ii][jj].blkm; l++)
									//for(m=0; m<sqrblkn[msid][k][ii][jj]; m++)
									for(m=0; m<mirror[msid].sqrblocks[k].blk[ii][jj].blkn; m++)
									//for(m=0; m<mirror[msid].sqrblocks[k].blk[ii][jj].blkn; m++)
										#pragma omp atomic write
										mirror[msid].sqrblocks[k].blk[ii][jj].blkelement[l][m]=mata[blockspos[k*morder[msid]+ii]+l][blockspos[k*morder[msid]+jj]+m];
										//msquare[msid][k][ii][jj][l][m]=mata[blockspos[k*morder[msid]+ii]+l][blockspos[k*morder[msid]+jj]+m];
										//mirror[msid].sqrblocks[k].blk[ii][jj].blkelement[l][m]=mata[blockspos[k*mirror[msid].morder+ii]+l][blockspos[k*mirror[msid].morder+jj]+m];
							}
						}
					}
					#pragma omp for collapse(4) private(k,ii,jj,kk,l,m,n) reduction(+:temp)
					for(k=0;k<(mblocksize[msid]/2);k++)
					{
						//for(ii=0; ii<mirror[msid].rpartnerblocks[k].brows; ii++)
						for(ii=0; ii<morder[msid]; ii++)
						{
							//for(jj=0; jj<mirror[msid].rpartnerblocks[k].bcols; jj++)
							for(jj=0; jj<morder[msid]; jj++)
							{
								for(kk=0; kk<morder[msid]; kk++)
								{
									//printf("\nmsid=%d, k=%d, ii=%d, jj=%d kk=%d\n",msid,k,ii,jj,kk);
									//for(l=0; l<sqrblkm[msid][2*k][jj][kk]; l++)
									for(l=0; l<mirror[msid].sqrblocks[2*k].blk[jj][kk].blkm; l++)
										//for(m=0; m<sqrblkn[msid][2*k][jj][kk]; m++)
										for(m=0; m<mirror[msid].sqrblocks[2*k].blk[jj][kk].blkn; m++)
											{
											//temp=msquare[msid][2*k][jj][kk][l][m];
											temp=0.0;//msquare[msid][2*k][jj][kk][l][m];
											//for(n=0; n<lpartblkm[msid][k][(jj+ii)% morder[msid]][kk]; n++)
											for(n=0; n<mirror[msid].lpartnerblocks[k].blk[(jj+ii)% morder[msid]][kk].blkm; n++)
											//for(n=0; n<mirror[msid].lpartnerblocks[k].blk[(jj+ii)%mirror[msid].morder][kk].blkm; n++)
											//temp+=floor(PRECP*(mata[blockspos[(k*2*morder[msid])+jj]+l][blockspos[((k*2+1)*morder[msid])+(jj+ii)% morder[msid]]+n] * mlpartner[msid][k][(jj+ii) % morder[msid]][kk][n][m]))/PRECP;
											temp+=(mata[blockspos[(k*2*morder[msid])+jj]+l][blockspos[((k*2+1)*morder[msid])+(jj+ii)% morder[msid]]+n] * mirror[msid].lpartnerblocks[k].blk[(jj+ii) % morder[msid]][kk].blkelement[n][m]);
										        #pragma omp atomic update
											mirror[msid].sqrblocks[2*k].blk[jj][kk].blkelement[l][m]+=temp;
											//msquare[msid][2*k][jj][kk][l][m]+=temp;
											}
											//mirror[msid].sqrblocks[2*k].blk[jj][kk].blkelement[l][m]+=floor(PRECP*(mata[blockspos[(k*2*mirror[msid].morder)+jj]+l][blockspos[((k*2+1)*mirror[msid].morder)+(jj+ii)%mirror[msid].morder]+n] * mirror[msid].lpartnerblocks[k].blk[(jj+ii)%mirror[msid].morder][kk].blkelement[n][m]))/PRECP;
								//mirror[msid].mirrormat[mirror[msid].mblockspos[k]+l][mirror[msid].mblockspos[k+1]+m]+=floor(PRECP*(-1.0*inverta[mirror[msid].mblockspos[k]+l][mirror[msid].mblockspos[k]+n]*mata[mirror[msid].mblockspos[k]+n][mirror[msid].mblockspos[k+1]+m]))/PRECP;
								}
							}
						}
					}
					#pragma omp for collapse(4) private(k,ii,jj,kk,l,m,n) reduction(+:temp)
					for(k=0;k<(mblocksize[msid]/2);k++)
					{
						//for(ii=0; ii<mirror[msid].rpartnerblocks[k].brows; ii++)
						for(ii=0; ii<morder[msid]; ii++)
						{
							//for(jj=0; jj<mirror[msid].rpartnerblocks[k].bcols; jj++)
							for(jj=0; jj<morder[msid]; jj++)
							{
								for(kk=0; kk<morder[msid]; kk++)
								{
									//for(l=0; l<sqrblkm[msid][2*k+1][jj][kk]; l++)
									for(l=0; l<mirror[msid].sqrblocks[2*k+1].blk[jj][kk].blkm; l++)
										//for(m=0; m < sqrblkn[msid][2*k+1][jj][kk]; m++)
										for(m=0; m < mirror[msid].sqrblocks[2*k+1].blk[jj][kk].blkn; m++)
											{
											//temp=msquare[msid][2*k+1][jj][kk][l][m];
											temp=0.0;//msquare[msid][2*k+1][jj][kk][l][m];
											//for(n=0; n < rpartblkm[msid][k][(jj+ii)%morder[msid]][kk]; n++)
											for(n=0; n < mirror[msid].rpartnerblocks[k].blk[(jj+ii)%morder[msid]][kk].blkm; n++)
											//for(n=0; n<mirror[msid].rpartnerblocks[k].blk[(jj+ii)%mirror[msid].morder][kk].blkm; n++)
											//temp+=floor(PRECP*(mata[blockspos[((k*2+1)* morder[msid])+jj]+l][blockspos[((k*2)* morder[msid])+(jj+ii)% morder[msid]]+n] * mrpartner[msid][k][(jj+ii)% morder[msid]][kk][n][m]))/PRECP;
											temp+=(mata[blockspos[((k*2+1)* morder[msid])+jj]+l][blockspos[((k*2)* morder[msid])+(jj+ii)% morder[msid]]+n] * mirror[msid].rpartnerblocks[k].blk[(jj+ii)% morder[msid]][kk].blkelement[n][m]);
										        #pragma omp atomic update
											mirror[msid].sqrblocks[2*k+1].blk[jj][kk].blkelement[l][m]+=temp;
											//msquare[msid][2*k+1][jj][kk][l][m]+=temp;
											}
											//mirror[msid].sqrblocks[2*k+1].blk[jj][kk].blkelement[l][m]+=floor(PRECP*(mata[blockspos[((k*2+1)*mirror[msid].morder)+jj]+l][blockspos[((k*2)*mirror[msid].morder)+(jj+ii)%mirror[msid].morder]+n] * mirror[msid].rpartnerblocks[k].blk[(jj+ii)%mirror[msid].morder][kk].blkelement[n][m]))/PRECP;
								}
							}
						}
					}
				}

			} //end of j==0 condition  
			
			else
			{ 
				if(loopid[i][j]==1)
				{
					//1. Performing the Inverse for S_A and S_D stored in the mirror (id=loopid[j-1]-2)
					//#pragma atomic
					//#pragma critical
					//{
					//#pragma omp barrier
					mcid=loopid[i][j-1]-2;
					//#pragma omp barrier
					//}
					//#pragma omp for collapse(2) private(k,ii)
					//#pragma omp for private(k,ii)
					#pragma omp for collapse(2) private(k,ii)
					for(k=0;k<(mblocksize[mcid]);k++)
					{
						for(ii=0; ii<morder[mcid]; ii++)
						{
							//printf("invertmirror: \tmcid=%d\tk=%d\tii=%d\tmirror[mcid].morder=%d\tk*mirror[mcid].morder+ii=%d\tblocks=%d\tblockspos=%d\n",mcid,k,ii,mirror[mcid].morder, k*mirror[mcid].morder+ii,blocks[k*mirror[mcid].morder+ii],blockspos[k*mirror[mcid].morder+ii]);
							
							//invertcases(blocks[k*morder[mcid]+ii], 0, msquare[mcid][k][ii][ii], blockspos[k * morder[mcid]+ii], inverta);
							invertcases(blocks[k*morder[mcid]+ii], 0, mirror[mcid].sqrblocks[k].blk[ii][ii].blkelement, blockspos[k * morder[mcid]+ii], inverta);
						}
					} 
					//#pragma omp for private(k)
					///for(k=0;k<blocksize;k++)
					///invertcases(blocks[k], blockspos[k], mirror[mcid].mirrormat, inverta);

					//2. Up Down arrow calculation 1->2, 2->3, 3->4, ...
					//#pragma atomic 
					//#pragma critical
					//#pragma omp barrier
 
					for(udmc=j,mcid=-1;udmc>0;udmc--)
					{
						if(loopid[i][udmc]==(loopid[i][udmc-1]-1)) mcid=j-udmc;
						else udmc=-1;
					}
					
					//printf("\n Up Down arrow decision with ud = %d: mcid=%d",ud, mcid);
					//#pragma omp barrier
					//threads for all the columns of the inverted matrix
					//#pragma barrier	
					
					if(mcid>=0)
					{
					#pragma omp for private(k,ii,jj,kk,msiditr,l,m,n,ud, itr,itrn) reduction(+:temp1,temp2,temp3,temp4)
					for(k=0;k<blocksize;k++)
					//for(k=(mcid>=0)?0:blocksize;k<blocksize;k++)
					//for(k=(mcid>=0)?0:order;k<blocksize;k++)
					{
						//printf("\n We are performing up down arrow upto mcid=%d with k=%d", mcid,k);
						//loop for number of iterations: 2^mcid.  i.e. if mcid=0,1,2,3,... then no of iterations=1,2,4,8,...
						for(itr=0;itr<=mcid;itr++)
						//for(itr=0;itr<((int)pow(2,mcid));itr++)
						{
							//printf("column k=%d\titeration count itr = %d\tmcid=%d\n",k,itr,mcid);
							//printf("mcid=%d\tk=%d\tii=%d\tjj=%d\n",mcid,k,ii,jj);
							if(itr==0) 
								jj=k;
							else 
							{ //itrn=itr-1;
                                                                if(((int)(k/morder[itr-1]))%2==0)
                                                                        jj=((int)(k/morder[itr-1]))*morder[itr-1]+morder[itr-1];
                                                                else
                                                                        jj=((int)(k/morder[itr-1]))*morder[itr-1]-morder[itr-1];
                                                       
							//	printf("\n in itr=%d testing morder of itr =%d",itr, morder[itr-1]);
							}


						//	printf("column k=%d\titeration count itr = %d\tmcid=%d\tmorder[itr-1]=%d\tjj=%d\n",k,itr,mcid,morder[itr-1],jj);
									//jj=(((int)((int)(k/((int)pow(2,itr-1)))*((int)pow(2,itr-1)))%2)==0)?((int)(k/((int)pow(2,itr-1)))*((int)pow(2,itr-1)))+((int)pow(2,itr-1))+(k%((int)pow(2,itr-1))):((int)(k/((int)pow(2,itr-1)))*((int)pow(2,itr-1)))-((int)pow(2,itr-1))+(k%((int)pow(2,itr-1)));

							//for(msid=itr;msid<=mcid;msid++)


							//selection I block which has to be multiplied.
							//we are keeping the information at kk.  So, (kk,k) is the location of I block
							if(itr==0) 
							{	//jj=k;
									for(msiditr=itr;msiditr<=mcid;msiditr++) //only number of blocks in each column
									{
										//(ii,k)is the place in inverta to store.
										//ii=((k/((int)pow(2,msid)))%2==0)?k+((int)pow(2,msid))-(k%((int)pow(2,msid))):k-((int)pow(2,msid))-(k%((int)pow(2,msid)));
										//ii=((k/((int)pow(2,msid)))%2==0)?k+((int)pow(2,msid))-(k%((int)pow(2,msid))):k-((int)pow(2,msid))-(k%((int)pow(2,msid)));
										if(((int)(k/morder[msiditr]))%2==0)
											ii=((int)(k/morder[msiditr]))*morder[msiditr]+morder[msiditr];
										else
											ii=((int)(k/morder[msiditr]))*morder[msiditr]-morder[msiditr];
										//ii=k+(int)pow(2,msid);
										//ii is block selection at k th column.  ii itself is the row information based on ii and k.
										//so, Storing at the block (ii,k).  ii has to be looped to order of msid.
										for(kk=0;kk<(morder[msiditr]);kk++)
												for(l=0;l<blocks[ii+kk];l++)
												for(m=0;m<blocks[k];m++)
													#pragma omp atomic update
													inverta[blockspos[ii+kk]+l][blockspos[k]+m]*=0.0;
										for(kk=0;kk<morder[msiditr];kk++)
										{
											//mirror column is I block location%msid order
											//mirror row is looped to msid order
											//storing: ii+kk, k.  kk is blk row.
											//if((k/((int)pow(2,msid)))%2==1)
											if(ii<jj)
											{
												//printf("\nBlock choice is rpartner blk=%d for k=%d, msid=%d\n",(int)k/((int)pow(2,msid+1)),k,msid);
												//blkchoice=(int)k/((int)pow(2,msid+1));
												//blkchoice=(int)k/(morder[msid]*2);
												//printf("\nBlock choice is rpartner blk=%d for k=%d, msid=%d\n",blkchoice,k,msid);
												//for(l=0;l< rpartblkm[msid][(int)k/((int)pow(2,msid+1))][kk][jj%((int)pow(2,msid))];l++)
												for(l=0;l<blocks[ii+kk];l++)
												//for(l=0;l<mirror[msid].rpartnerblocks[k/((int)pow(2,msid+1))].blk[0][0].blkm;l++)//.blk[kk][jj%((int)pow(2,msid))].blkm;l++)
												//for(l=0;l<mirror[msid].rpartnerblocks[blkchoice].blk[kk][jj%((int)pow(2,msid))].blkm;l++)
												for(m=0;m<blocks[k];m++)
												{
													temp1=0.0;
													//inverta[blockspos[ii+kk]+l][blockspos[k]+m]=0.0;
												for(n=0;n<blocks[jj];n++)
												//for(n=0;n< rpartblkn[msid][(int)k/((int)pow(2,msid+1))][kk][jj%((int)pow(2,msid))];n++)
												{
													///inverta[blockspos[ii+kk]+l][blockspos[k]+m]+=floor(PRECP*(mrpartner[msiditr][(int)k/((int)pow(2,msiditr+1))][kk][jj%((int)pow(2,msiditr))][l][n] * inverta[blockspos[jj]+n][blockspos[k]+m]))/PRECP;
													//temp1+=floor(PRECP*(mrpartner[msiditr][(int)k/((int)pow(2,msiditr+1))][kk][jj%((int)pow(2,msiditr))][l][n] * inverta[blockspos[jj]+n][blockspos[k]+m]))/PRECP;
													temp1+=(mirror[msiditr].rpartnerblocks[(int)k/((int)pow(2,msiditr+1))].blk[kk][jj%((int)pow(2,msiditr))].blkelement[l][n] * inverta[blockspos[jj]+n][blockspos[k]+m]);
													//printf("\n Check mult: mirror msid=%d rpartner=%d blk(%d,%d) blkelement=%lf, \tinverta location=(%d,%d)\tvalue=%lf\tstoring:inva(%d,%d)", msid,(int)k/((int)pow(2,msid+1)),kk,jj%((int)pow(2,msid)),mirror[msid].rpartnerblocks[(int)k/((int)pow(2,msid+1))].blk[kk][jj%((int)pow(2,msid))].blkelement[l][n],blockspos[jj]+n,blockspos[k]+m,inverta[blockspos[jj]+n][blockspos[k]+m],blockspos[ii+kk]+l,blockspos[k]+m);
												}
													#pragma omp atomic update
													inverta[blockspos[ii+kk]+l][blockspos[k]+m]+=temp1;
												}
											}
											else
											{
												//for(l=0;l < lpartblkm[msid][(int)k/((int)pow(2,msid+1))][kk][jj%((int)pow(2,msid))];l++)
												for(l=0;l<blocks[ii+kk];l++)
												//for(l=0;l<mirror[msid].rpartnerblocks[k/((int)pow(2,msid+1))].blk[0][0].blkm;l++)//.blk[kk][jj%((int)pow(2,msid))].blkm;l++)
												//for(l=0;l<mirror[msid].rpartnerblocks[blkchoice].blk[kk][jj%((int)pow(2,msid))].blkm;l++)
												for(m=0;m < blocks[k];m++)
												{
													temp2=0.0;
													//inverta[blockspos[ii+kk]+l][blockspos[k]+m]=0.0;
												//for(n=0;n<blocks[jj];n++)
												for(n=0;n<blocks[jj];n++)
												//for(n=0;n < lpartblkn[msid][(int)k/((int)pow(2,msid+1))][kk][jj%((int)pow(2,msid))];n++)
												{
													///inverta[blockspos[ii+kk]+l][blockspos[k]+m]+=floor(PRECP*( mlpartner[msiditr][(int)k/((int)pow(2,msiditr+1))][kk][jj%((int)pow(2,msiditr))][l][n] * inverta[blockspos[jj]+n][blockspos[k]+m]))/PRECP;
													//temp2+=floor(PRECP*( mlpartner[msiditr][(int)k/((int)pow(2,msiditr+1))][kk][jj%((int)pow(2,msiditr))][l][n] * inverta[blockspos[jj]+n][blockspos[k]+m]))/PRECP;
													temp2+=( mirror[msiditr].lpartnerblocks[(int)k/((int)pow(2,msiditr+1))].blk[kk][jj%((int)pow(2,msiditr))].blkelement[l][n] * inverta[blockspos[jj]+n][blockspos[k]+m]);
													//printf("\n Check mult: mirror msid=%d lpartner=%d blk(%d,%d) blkelement=%lf, \tinverta location=(%d,%d)\tvalue=%lf\tstoring:inva(%d,%d)", msid,(int)k/((int)pow(2,msid+1)),kk,jj%((int)pow(2,msid)),mirror[msid].lpartnerblocks[(int)k/((int)pow(2,msid+1))].blk[kk][jj%((int)pow(2,msid))].blkelement[l][n],blockspos[jj]+n,blockspos[k]+m,inverta[blockspos[jj]+n][blockspos[k]+m],blockspos[ii+kk]+l,blockspos[k]+m);
												//for(l=0;l<mirror[msid].lpartnerblocks[(int)k/((int)pow(2,msid+1))].blk[kk][jj%((int)pow(2,msid))].blkm;l++)
												//for(m=0;m<blocks[k];m++)
												//{
												//	inverta[blockspos[ii+kk]+l][blockspos[k]+m]=0.0;
												//for(n=0;n<blocks[jj];n++)
												//{
												//	inverta[blockspos[ii+kk]+l][blockspos[k]+m]+=mirror[msid].lpartnerblocks[(int)k/((int)pow(2,msid+1))].blk[kk][jj%((int)pow(2,msid))].blkelement[l][n]*inverta[blockspos[jj]+n][blockspos[k]+m];
												}
													#pragma omp atomic update
													inverta[blockspos[ii+kk]+l][blockspos[k]+m]+=temp2;
												}
												//for(l=0;mirror[msid].lpartnerblock[].blk[][].blkm;l++)
											}
										}	
									//printf("\t check in itr=0 k=%d loop mcid=%d\t msid =%d\t jj=%d\tii=%d\t\n",k, mcid,msid,jj,ii);
									//selection of I block to get multiplied. Based on k and itr
	
									//column for mirrors is row of I in the virtual block.
									}
							}
							else 
							{
								//jj=((int)(k/((int)pow(2,itr-1)))%2==0)?k+((int)pow(2,itr-1))-(k%((int)pow(2,itr-1))):k-((int)pow(2,itr-1))-(k%((int)pow(2,itr-1)));
								//jj=(((int)((int)(k/((int)pow(2,itr-1)))*((int)pow(2,itr-1)))%2)==0)?((int)(k/((int)pow(2,itr-1)))*((int)pow(2,itr-1)))+((int)pow(2,itr-1))+(k%((int)pow(2,itr-1))):((int)(k/((int)pow(2,itr-1)))*((int)pow(2,itr-1)))-((int)pow(2,itr-1))+(k%((int)pow(2,itr-1)));
								for(ud=0; ud<((int)pow(2,itr-1));ud++)
								{
									//kk= ((k%((int)pow(2,itr)))%2==0)?k+((int)pow(2,itr))-(((int)pow(2,itr))):k-((int)pow(2,itr))-(((int)pow(2,itr)));
									//jj+=ud;   ///rts we cannot change jj inthis place. it gives segmentation.
									//for(ii=0;ii<((int)(pow(2,mcid)*2));ii++) //only number of blocks in each column
									for(msiditr=itr;msiditr<=mcid;msiditr++) //only number of blocks in each column
									{
										//(ii,k)is the place in inverta to store.
										////ii=((k/((int)pow(2,msid)))%2==0)?k+((int)pow(2,msid))-(k%((int)pow(2,msid))):k-((int)pow(2,msid))-(k%((int)pow(2,msid)));
										if(((int)(k/morder[msiditr]))%2==0)
											ii=((int)(k/morder[msiditr]))*morder[msiditr]+morder[msiditr];
										else
											ii=((int)(k/morder[msiditr]))*morder[msiditr]-morder[msiditr];
										//ii is block selection at k th column.  ii itself is the row information based on ii and k.
										//so, Storing at the block (ii,k).  ii has to be looped to order of msid.
										for(kk=0;kk<((int)pow(2,msiditr));kk++)
										{
											//mirror column is I block location%msid order
											//mirror row is looped to msid order
											if(ii<jj)
											{
												//for(l=0;l < rpartblkm[msid][(int)k/((int)pow(2,msid+1))][kk][jj%((int)pow(2,msid))];l++)
												for(l=0;l < blocks[ii+kk];l++)
												for(m=0;m<blocks[k];m++)
												{
													temp3=0.0;
												//for(n=0;n<blocks[jj];n++)
												//for(n=0;n < rpartblkn[msiditr][(int)k/((int)pow(2,msiditr+1))][kk][(jj+ud)%((int)pow(2,msiditr))];n++)
												for(n=0;n < mirror[msiditr].rpartnerblocks[(int)k/((int)pow(2,msiditr+1))].blk[kk][(jj+ud)%((int)pow(2,msiditr))].blkn;n++)
												{
													//inverta[blockspos[ii+kk]+l][blockspos[k]+m]+=floor(PRECP*( mrpartner[msid][(int)k/((int)pow(2,msid+1))][kk][jj%((int)pow(2,msid))][l][n] * inverta[blockspos[jj]+n][blockspos[k]+m]))/PRECP;
													///inverta[blockspos[ii+kk]+l][blockspos[k]+m]+=floor(PRECP*( mrpartner[msiditr][(int)(k/((int)pow(2,msiditr+1)))][kk][(int)((jj+ud)%((int)pow(2,msiditr)))][l][n] * inverta[blockspos[jj+ud]+n][blockspos[k]+m]))/PRECP;
													//temp3+=floor(PRECP*( mrpartner[msiditr][(int)(k/((int)pow(2,msiditr+1)))][kk][(int)((jj+ud)%((int)pow(2,msiditr)))][l][n] * inverta[blockspos[jj+ud]+n][blockspos[k]+m]))/PRECP;
													temp3+=( mirror[msiditr].rpartnerblocks[(int)(k/((int)pow(2,msiditr+1)))].blk[kk][(int)((jj+ud)%((int)pow(2,msiditr)))].blkelement[l][n] * inverta[blockspos[jj+ud]+n][blockspos[k]+m]);
													//printf("\n Testing: %d %d %d %d %d %d", blockspos[ii+kk], blockspos[k], (int)(k/((int)pow(2,msid+1))),(int)(jj%((int)pow(2,msid))),blockspos[jj],blockspos[k]);
													//printf("\n Testing: %lf", mrpartner[msid][(int)(k/((int)pow(2,msid+1)))][kk][(int)(jj%((int)pow(2,msid)))][l][n]);// * inverta[blockspos[jj]+n][blockspos[k]+m]))/PRECP;
													//printf("\n Testing: %lf \t %d\t%d", inverta[blockspos[jj]+n][blockspos[k]+m], blockspos[jj]+n,blockspos[k]+m);
													//printf("\n Testing: jj+ud=%d blockspos[jj+ud]=%d \t %d\t%d", jj+ud,blockspos[jj+ud], blockspos[jj+ud]+n,blockspos[k]+m);
												}
													#pragma omp atomic update
													inverta[blockspos[ii+kk]+l][blockspos[k]+m]+=temp3;
												}
											}
											else
											{
												//for(l=0;l< lpartblkm[msid][(int)k/((int)pow(2,msid+1))][kk][jj%((int)pow(2,msid))];l++)
												for(l=0;l< blocks[ii+kk];l++)
												for(m=0;m<blocks[k];m++)
												{
													temp4=0.0;
												//for(n=0;n<blocks[jj];n++)
												//for(n=0;n< lpartblkn[msiditr][(int)k/((int)pow(2,msiditr+1))][kk][(jj+ud)%((int)pow(2,msiditr))];n++)
												for(n=0;n< mirror[msiditr].lpartnerblocks[(int)k/((int)pow(2,msiditr+1))].blk[kk][(jj+ud)%((int)pow(2,msiditr))].blkn;n++)
												{
													//inverta[blockspos[ii+kk]+l][blockspos[k]+m]+=floor(PRECP*(mlpartner[msid][(int)k/((int)pow(2,msid+1))][kk][jj%((int)pow(2,msid))][l][n] * inverta[blockspos[jj]+n][blockspos[k]+m]))/PRECP;
													///inverta[blockspos[ii+kk]+l][blockspos[k]+m]+=floor(PRECP*(mlpartner[msiditr][(int)(k/((int)pow(2,msiditr+1)))][kk][(int)((jj+ud)%((int)pow(2,msiditr)))][l][n] * inverta[blockspos[jj+ud]+n][blockspos[k]+m]))/PRECP;
													//temp4+=floor(PRECP*(mlpartner[msiditr][(int)(k/((int)pow(2,msiditr+1)))][kk][(int)((jj+ud)%((int)pow(2,msiditr)))][l][n] * inverta[blockspos[jj+ud]+n][blockspos[k]+m]))/PRECP;
													temp4+=(mirror[msiditr].lpartnerblocks[(int)(k/((int)pow(2,msiditr+1)))].blk[kk][(int)((jj+ud)%((int)pow(2,msiditr)))].blkelement[l][n] * inverta[blockspos[jj+ud]+n][blockspos[k]+m]);
												}
													//inverta[blockspos[ii+kk]+l][blockspos[k]+m]+=floor(PRECP*(mirror[msid].lpartnerblocks[(int)k/((int)pow(2,msid+1))].blk[kk][jj%((int)pow(2,msid))].blkelement[l][n]*inverta[blockspos[jj]+n][blockspos[k]+m]))/PRECP;
												//for(l=0;mirror[msid].lpartnerblock[].blk[][].blkm;l++)
													#pragma omp atomic update
													inverta[blockspos[ii+kk]+l][blockspos[k]+m]+=temp4;
													
												}
											}
										}	
									//printf("\t in ii loop mcid=%d\tk=%d\tii=%d\tjj=%d\n",mcid,k,ii,jj);
									//selection of I block to get multiplied. Based on k and itr
	
									//column for mirrors is row of I in the virtual block.
									}
								}
							}
						}
					 }  } 
				} 
				else
				{
					//1. Calculation of -A^-1B and -D^-1C using mirror (id=loopid[j-1]-2) (for B and D) and storing at Mirror (id=loopid[j]-2), 
					//#pragma atomic
					//#pragma critical
					//{
					//mirrorchoiceid
					//#pragma omp barrier
					mcid = loopid[i][j-1]-2;
					//mirrorstoreid
					msid = loopid[i][j]-2;
					//}
					//#pragma omp barrier
					#pragma omp for collapse(3) private(k,ii,jj,l,m)
					for(k=0;k<(mblocksize[msid]/2);k++)
					{
						for(ii=0; ii<morder[msid]; ii++)
						{
							for(jj=0; jj<morder[msid]; jj++)
							{
								//for(l=0; l< rpartblkm[msid][k][ii][jj]; l++)
								for(l=0; l< mirror[msid].rpartnerblocks[k].blk[ii][jj].blkm; l++)
									//for(m=0; m< rpartblkn[msid][k][ii][jj]; m++)
									for(m=0; m< mirror[msid].rpartnerblocks[k].blk[ii][jj].blkn; m++)
										//mrpartner[msid][k][ii][jj][l][m]=0.0;
										mirror[msid].rpartnerblocks[k].blk[ii][jj].blkelement[l][m]=0.0;
							}
						}
					}
					#pragma omp for collapse(3) private(k,ii,jj,l,m)
					for(k=0;k<(mblocksize[msid]/2);k++)
					{
						for(ii=0; ii<morder[msid]; ii++)
						{
							for(jj=0; jj<morder[msid]; jj++)
							{
								//for(l=0; l< lpartblkm[msid][k][ii][jj]; l++)
								for(l=0; l< mirror[msid].lpartnerblocks[k].blk[ii][jj].blkm; l++)
									//for(m=0; m< lpartblkn[msid][k][ii][jj]; m++)
									for(m=0; m< mirror[msid].lpartnerblocks[k].blk[ii][jj].blkn; m++)
										//mlpartner[msid][k][ii][jj][l][m]=0.0;
										mirror[msid].lpartnerblocks[k].blk[ii][jj].blkelement[l][m]=0.0;
							}
						}
					}
					#pragma omp for collapse(4) private(k,ii,jj,kk,l,m,n) reduction(+:temp)
					for(k=0;k<(mblocksize[msid]/2);k++)
					{
						//for(ii=0; ii<mirror[msid].rpartnerblocks[k].brows; ii++)
						for(ii=0; ii<morder[msid]; ii++)
						{
							//for(jj=0; jj<mirror[msid].rpartnerblocks[k].bcols; jj++)
							for(jj=0; jj<morder[msid]; jj++)
							{
								for(kk=0; kk< morder[msid]; kk++)
								{
									//printf("\nmsid=%d, k=%d, ii=%d, jj=%d kk=%d\n",msid,k,ii,jj,kk);
									//for(l=0; l< rpartblkm[msid][k][jj][kk]; l++)
									for(l=0; l< mirror[msid].rpartnerblocks[k].blk[jj][kk].blkm; l++)
										//for(m=0; m<mirror[msid].rpartnerblocks[k].blk[jj][kk].blkn; m++)
											//for(n=0; n<mirror[msid].rpartnerblocks[k].blk[jj][kk].blkm; n++)
											//for(n=0; n<blocks[((k*2)*mirror[msid].morder)+(jj+ii)%mirror[msid].morder];n++)
											//for(m=0;m< sqrblkn[mcid][(int)((2*k*morder[msid])/morder[mcid])][((2*k*morder[msid])%morder[mcid])+(jj+ii)%morder[msid]][(((2*k)* morder[msid])% morder[mcid])+ morder[msid]+kk];m++)
											for(m=0;m< mirror[mcid].sqrblocks[(int)((2*k*morder[msid])/morder[mcid])].blk[((2*k*morder[msid])%morder[mcid])+(jj+ii)%morder[msid]][(((2*k)* morder[msid])% morder[mcid])+ morder[msid]+kk].blkn;m++)
											{
											temp=0.0;
											//for(m=0;m<mirror[mcid].sqrblocks[(int)((2*k*mirror[msid].morder)/mirror[mcid].morder)].blk[((2*k*mirror[msid].morder)%mirror[mcid].morder)+(jj+ii)%mirror[msid].morder][(((2*k)*mirror[msid].morder)%mirror[mcid].morder)+mirror[msid].morder+kk].blkn;m++)
											//for(n=0;n< sqrblkm[mcid][(int)((2*k*morder[msid])/morder[mcid])][((2*k*morder[msid])%morder[mcid])+(jj+ii)%morder[msid]][(((2*k)*morder[msid])% morder[mcid])+ morder[msid]+kk];n++)
											for(n=0;n< mirror[mcid].sqrblocks[(int)((2*k*morder[msid])/morder[mcid])].blk[((2*k*morder[msid])%morder[mcid])+(jj+ii)%morder[msid]][(((2*k)*morder[msid])% morder[mcid])+ morder[msid]+kk].blkm;n++)
											//for(n=0;n<mirror[mcid].sqrblocks[(int)((2*k*mirror[msid].morder)/mirror[mcid].morder)].blk[((2*k*mirror[msid].morder)%mirror[mcid].morder)+(jj+ii)%mirror[msid].morder][(((2*k)*mirror[msid].morder)%mirror[mcid].morder)+mirror[msid].morder+kk].blkm;n++)
											///mrpartner[msid][k][jj][kk][l][m]+=floor(PRECP*(inverta[blockspos[(k*2* morder[msid])+jj]+l][blockspos[(k*2*morder[msid])+(jj+ii)% morder[msid]]+n]* msquare[mcid][(int)((2*k*morder[msid])/morder[mcid])][((2*k*morder[msid])%morder[mcid])+(jj+ii)%morder[msid]][((2*k*morder[msid])%morder[mcid])+morder[msid]+kk][n][m]*(-1.0)))/PRECP; ///rrrts change k to k+1.
											//temp+=floor(PRECP*(inverta[blockspos[(k*2* morder[msid])+jj]+l][blockspos[(k*2*morder[msid])+(jj+ii)% morder[msid]]+n]* msquare[mcid][(int)((2*k*morder[msid])/morder[mcid])][((2*k*morder[msid])%morder[mcid])+(jj+ii)%morder[msid]][((2*k*morder[msid])%morder[mcid])+morder[msid]+kk][n][m]*(-1.0)))/PRECP; ///rrrts change k to k+1.
											temp+=(inverta[blockspos[(k*2* morder[msid])+jj]+l][blockspos[(k*2*morder[msid])+(jj+ii)% morder[msid]]+n]* mirror[mcid].sqrblocks[(int)((2*k*morder[msid])/morder[mcid])].blk[((2*k*morder[msid])%morder[mcid])+(jj+ii)%morder[msid]][((2*k*morder[msid])%morder[mcid])+morder[msid]+kk].blkelement[n][m]*(-1.0)); ///rrrts change k to k+1.
											//mirror[msid].rpartnerblocks[k].blk[jj][kk].blkelement[l][m]+=floor(PRECP*(inverta[blockspos[(k*2*mirror[msid].morder)+jj]+l][blockspos[(k*2*mirror[msid].morder)+(jj+ii)%mirror[msid].morder]+n]*mirror[mcid].sqrblocks[(int)((2*k*mirror[msid].morder)/mirror[mcid].morder)].blk[((2*k*mirror[msid].morder)%mirror[mcid].morder)+(jj+ii)%mirror[msid].morder][(((2*k+1)*mirror[msid].morder)%mirror[mcid].morder)+kk].blkelement[n][m]*(-1.0)))/PRECP; ///rrrts change k to k+1.
//.blk[(k%mirror[mcid].morder)*mirror[msid].morder][(k%mirror[mcid].morder+1)*mirror[msid].morder].blkelement[n][m]*(-1.0);
								// *mata[blockspos[(k*2*mirror[msid].morder)+(jj+ii)%mirror[msid].morder]+n][blockspos[((k*2+1)*mirror[msid].morder)+kk]+m]*(-1.0);
											#pragma omp atomic update
											mirror[msid].rpartnerblocks[k].blk[jj][kk].blkelement[l][m]+=temp;
											//mrpartner[msid][k][jj][kk][l][m]+=temp;
								 			}
								}
							}
						}
					}
					#pragma omp for collapse(4) private(k,ii,jj,kk,l,m,n) reduction(+:temp)
					for(k=0;k<(mblocksize[msid]/2);k++)
					{
						//for(ii=0; ii<mirror[msid].rpartnerblocks[k].brows; ii++)
						for(ii=0; ii<morder[msid]; ii++)
						{
							//for(jj=0; jj<mirror[msid].rpartnerblocks[k].bcols; jj++)
							for(jj=0; jj<morder[msid]; jj++)
							{
								for(kk=0; kk< morder[msid]; kk++)
								{
									//printf("\nmsid=%d, k=%d, ii=%d, jj=%d kk=%d\n",msid,k,ii,jj,kk);
									//for(l=0; l< lpartblkm[msid][k][jj][kk]; l++)
									for(l=0; l< mirror[msid].lpartnerblocks[k].blk[jj][kk].blkm; l++)
										//for(m=0; m<mirror[msid].lpartnerblocks[k].blk[jj][kk].blkn; m++)
											//for(n=0; n<mirror[msid].lpartnerblocks[k].blk[jj][kk].blkm; n++)
											//for(n=0; n<blocks[((k*2+1)*mirror[msid].morder)+(jj+ii)%mirror[msid].morder];n++)
											//for(m=0;m< sqrblkn[mcid][(int)(((2*k)*morder[msid])/morder[mcid])][(((2*k+1)*morder[msid])%morder[mcid])+(jj+ii)%morder[msid]][(((2*k+1)*morder[msid])%morder[mcid])-morder[msid]+kk];m++)
											for(m=0;m< mirror[mcid].sqrblocks[(int)(((2*k)*morder[msid])/morder[mcid])].blk[(((2*k+1)*morder[msid])%morder[mcid])+(jj+ii)%morder[msid]][(((2*k+1)*morder[msid])%morder[mcid])-morder[msid]+kk].blkn;m++)
											//for(m=0;m<mirror[mcid].sqrblocks[(int)(((2*k)*mirror[msid].morder)/mirror[mcid].morder)].blk[(((2*k+1)*mirror[msid].morder)%mirror[mcid].morder)+(jj+ii)%mirror[msid].morder][(((2*k+1)*mirror[msid].morder)%mirror[mcid].morder)-mirror[msid].morder+kk].blkn;m++)
											{
											temp=0.0;
											//for(n=0;n< sqrblkm[mcid][(int)((2*k*morder[msid])/morder[mcid])][(((2*k+1)*morder[msid])%morder[mcid])+(jj+ii)%morder[msid]][(((2*k+1)* morder[msid])%morder[mcid])-morder[msid]+kk];n++)
											for(n=0;n< mirror[mcid].sqrblocks[(int)((2*k*morder[msid])/morder[mcid])].blk[(((2*k+1)*morder[msid])%morder[mcid])+(jj+ii)%morder[msid]][(((2*k+1)* morder[msid])%morder[mcid])-morder[msid]+kk].blkm;n++)
											//for(n=0;n<mirror[mcid].sqrblocks[(int)((2*k*mirror[msid].morder)/mirror[mcid].morder)].blk[(((2*k+1)*mirror[msid].morder)%mirror[mcid].morder)+(jj+ii)%mirror[msid].morder][(((2*k+1)*mirror[msid].morder)%mirror[mcid].morder)-mirror[msid].morder+kk].blkm;n++)
											///mlpartner[msid][k][jj][kk][l][m]+=floor(PRECP*(inverta[blockspos[((k*2+1)*morder[msid])+jj]+l][blockspos[((k*2+1)*morder[msid])+(jj+ii)% morder[msid]]+n]* msquare[mcid][(int)(((2*k)*morder[msid])/morder[mcid])][(((2*k+1)*morder[msid])%morder[mcid])+(jj+ii)%morder[msid]][(((2*k+1)*morder[msid])%morder[mcid])-morder[msid]+kk][n][m]*(-1.0)))/PRECP;
											//temp+=floor(PRECP*(inverta[blockspos[((k*2+1)*morder[msid])+jj]+l][blockspos[((k*2+1)*morder[msid])+(jj+ii)% morder[msid]]+n]* msquare[mcid][(int)(((2*k)*morder[msid])/morder[mcid])][(((2*k+1)*morder[msid])%morder[mcid])+(jj+ii)%morder[msid]][(((2*k+1)*morder[msid])%morder[mcid])-morder[msid]+kk][n][m]*(-1.0)))/PRECP;
											temp+=(inverta[blockspos[((k*2+1)*morder[msid])+jj]+l][blockspos[((k*2+1)*morder[msid])+(jj+ii)% morder[msid]]+n]* mirror[mcid].sqrblocks[(int)(((2*k)*morder[msid])/morder[mcid])].blk[(((2*k+1)*morder[msid])%morder[mcid])+(jj+ii)%morder[msid]][(((2*k+1)*morder[msid])%morder[mcid])-morder[msid]+kk].blkelement[n][m]*(-1.0));
											//mirror[msid].lpartnerblocks[k].blk[jj][kk].blkelement[l][m]+=floor(PRECP*(inverta[blockspos[((k*2+1)*mirror[msid].morder)+jj]+l][blockspos[((k*2+1)*mirror[msid].morder)+(jj+ii)%mirror[msid].morder]+n]*mirror[mcid].sqrblocks[(int)((2*k*mirror[msid].morder)/mirror[mcid].morder)].blk[(((2*k+1)*mirror[msid].morder)%mirror[mcid].morder)+(jj+ii)%mirror[msid].morder][(((2*k)*mirror[msid].morder)%mirror[mcid].morder)+kk].blkelement[n][m]*(-1.0)))/PRECP;///rrrts change 2k+1 to 2k.
//.blk[(k%mirror[mcid].morder)*mirror[msid].morder][(k%mirror[mcid].morder+1)*mirror[msid].morder].blkelement[n][m]*(-1.0);
								// *mata[blockspos[(k*2*mirror[msid].morder)+(jj+ii)%mirror[msid].morder]+n][blockspos[((k*2+1)*mirror[msid].morder)+kk]+m]*(-1.0);
											//mirror[msid].lpartnerblocks[k].blk[jj][kk].blkelement[l][m]+=inverta[blockspos[((k*2+1)*mirror[msid].morder)+jj]+l][blockspos[((k*2+1)*mirror[msid].morder)+(jj+ii)%mirror[msid].morder]+n]*mata[blockspos[((k*2+1)*mirror[msid].morder)+(jj+ii)%mirror[msid].morder]+n][blockspos[((k*2)*mirror[msid].morder)+kk]+m]*(-1.0);
											#pragma omp atomic update
											mirror[msid].lpartnerblocks[k].blk[jj][kk].blkelement[l][m]+=temp;
											//mlpartner[msid][k][jj][kk][l][m]+=temp;
											}
								}
							}
						}
					}
					
					//#pragma omp barrier
					//2. Calculation of S_A and S_D for the mirror (id=loopid[j-1]-2)  at the location Mirror (id=loopid[j]-2)
					#pragma omp for collapse(3) private(k,ii,jj,l,m)
					for(k=0;k<(mblocksize[msid]);k++)
					{
						for(ii=0; ii<morder[msid]; ii++)
						{
							for(jj=0; jj<morder[msid]; jj++)
							{
								//for(l=0; l< sqrblkm[msid][k][ii][jj]; l++)
								for(l=0; l< mirror[msid].sqrblocks[k].blk[ii][jj].blkm; l++)
									//for(m=0; m< sqrblkn[msid][k][ii][jj]; m++)
									for(m=0; m< mirror[msid].sqrblocks[k].blk[ii][jj].blkn; m++)
										//msquare[msid][k][ii][jj][l][m]= msquare[mcid][(int)((k*morder[msid])/morder[mcid])][(((k)*morder[msid])%morder[mcid])+ii][(int)((k)*morder[msid])%morder[mcid]+jj][l][m];
										mirror[msid].sqrblocks[k].blk[ii][jj].blkelement[l][m]= mirror[mcid].sqrblocks[(int)((k*morder[msid])/morder[mcid])].blk[(((k)*morder[msid])%morder[mcid])+ii][(int)((k)*morder[msid])%morder[mcid]+jj].blkelement[l][m];
							}
						}
					}
					#pragma omp for collapse(4) private(k,ii,jj,kk,l,m,n) reduction(+:temp)
					for(k=0;k<(mblocksize[msid]/2);k++)
					{
						//for(ii=0; ii<mirror[msid].rpartnerblocks[k].brows; ii++)
						for(ii=0; ii< morder[msid]; ii++)
						{
							//for(jj=0; jj<mirror[msid].rpartnerblocks[k].bcols; jj++)
							for(jj=0; jj< morder[msid]; jj++)
							{
								for(kk=0; kk<morder[msid]; kk++)
								{
									//printf("\nmsid=%d, k=%d, ii=%d, jj=%d kk=%d\n",msid,k,ii,jj,kk);
									//for(l=0; l< sqrblkm[msid][2*k][jj][kk]; l++)
									for(l=0; l<mirror[msid].sqrblocks[2*k].blk[jj][kk].blkm; l++)
										//for(m=0; m<mirror[msid].sqrblocks[2*k].blk[jj][kk].blkn; m++)
										//for(m=0; m< lpartblkn[msid][k][(jj+ii)% morder[msid]][kk]; m++)
										for(m=0; m<mirror[msid].lpartnerblocks[k].blk[(jj+ii)%morder[msid]][kk].blkn; m++)
										{
											temp=0.0;
											//for(n=0; n< lpartblkm[msid][k][(jj+ii)%morder[msid]][kk]; n++)
											for(n=0; n<mirror[msid].lpartnerblocks[k].blk[(jj+ii)%morder[msid]][kk].blkm; n++)
											//for(n=0;n<mirror[mcid].sqrblocks[(int)(((2*k)*mirror[msid].morder)/mirror[mcid].morder)].blk[(((2*k)*mirror[msid].morder)%mirror[mcid].morder)+(jj)][((2*k)*mirror[msid].morder)%mirror[mcid].morder+mirror[msid].morder+(jj+ii)%mirror[msid].morder].blkn;n++)
											//mirror[msid].sqrblocks[2*k].blk[jj][kk].blkelement[l][m]+= floor(PRECP*(mirror[mcid].sqrblocks[(int)((2*k*mirror[msid].morder)/mirror[mcid].morder)].blk[(((2*k)*mirror[msid].morder)%mirror[mcid].morder)+(jj)][(((2*k)*mirror[msid].morder)%mirror[mcid].morder)+mirror[msid].morder+(jj+ii)%mirror[msid].morder].blkelement[l][n] * mirror[msid].lpartnerblocks[k].blk[(jj+ii)%mirror[msid].morder][kk].blkelement[n][m]))/PRECP;
											///msquare[msid][2*k][jj][kk][l][m]+= floor(PRECP*( msquare[mcid][(int)(((2*k)*morder[msid])/morder[mcid])][(((2*k)*morder[msid])%morder[mcid])+(jj)][((2*k)*morder[msid])%morder[mcid]+morder[msid]+(jj+ii)%morder[msid]][l][n] * mlpartner[msid][k][(jj+ii)%morder[msid]][kk][n][m]))/PRECP;
											//temp+= floor(PRECP*( msquare[mcid][(int)(((2*k)*morder[msid])/morder[mcid])][(((2*k)*morder[msid])%morder[mcid])+(jj)][((2*k)*morder[msid])%morder[mcid]+morder[msid]+(jj+ii)%morder[msid]][l][n] * mlpartner[msid][k][(jj+ii)%morder[msid]][kk][n][m]))/PRECP;
											temp+= ( mirror[mcid].sqrblocks[(int)(((2*k)*morder[msid])/morder[mcid])].blk[(((2*k)*morder[msid])%morder[mcid])+(jj)][((2*k)*morder[msid])%morder[mcid]+morder[msid]+(jj+ii)%morder[msid]].blkelement[l][n] * mirror[msid].lpartnerblocks[k].blk[(jj+ii)%morder[msid]][kk].blkelement[n][m]);
								//mirror[msid].mirrormat[mirror[msid].mblockspos[k]+l][mirror[msid].mblockspos[k+1]+m]+=floor(PRECP*(-1.0*inverta[mirror[msid].mblockspos[k]+l][mirror[msid].mblockspos[k]+n]*mata[mirror[msid].mblockspos[k]+n][mirror[msid].mblockspos[k+1]+m]))/PRECP;
											#pragma omp atomic update
											mirror[msid].sqrblocks[2*k].blk[jj][kk].blkelement[l][m]+=temp;
											//msquare[msid][2*k][jj][kk][l][m]+=temp;
										}
								}
							}
						}
					}
					#pragma omp for collapse(4) private(k,ii,jj,kk,l,m,n) reduction(+:temp)
					for(k=0;k<(mblocksize[msid]/2);k++)
					{
						//for(ii=0; ii<mirror[msid].rpartnerblocks[k].brows; ii++)
						for(ii=0; ii< morder[msid]; ii++)
						{
							//for(jj=0; jj<mirror[msid].rpartnerblocks[k].bcols; jj++)
							for(jj=0; jj< morder[msid]; jj++)
							{
								for(kk=0; kk<morder[msid]; kk++)
								{
									//printf("\nmsid=%d, k=%d, ii=%d, jj=%d kk=%d\n",msid,k,ii,jj,kk);
									//for(l=0; l< sqrblkm[msid][2*k+1][jj][kk]; l++) //rrrtschange added +1 with 2*k
									for(l=0; l< mirror[msid].sqrblocks[2*k+1].blk[jj][kk].blkm; l++) //rrrtschange added +1 with 2*k
									//for(l=0; l<mirror[msid].sqrblocks[2*k+1].blk[jj][kk].blkm; l++) //rrrtschange added +1 with 2*k
										//for(m=0; m<mirror[msid].sqrblocks[2*k+1].blk[jj][kk].blkn; m++) //rrrtschange
										//for(m=0; m< rpartblkn[msid][k][(jj+ii)%morder[msid]][kk]; m++)
										for(m=0; m< mirror[msid].rpartnerblocks[k].blk[(jj+ii)%morder[msid]][kk].blkn; m++)
										//for(m=0; m<mirror[msid].rpartnerblocks[k].blk[(jj+ii)%mirror[msid].morder][kk].blkn; m++)
										{
											temp=0.0;
											//for(n=0; n<mirror[msid].rpartnerblocks[k].blk[(jj+ii)%mirror[msid].morder][kk].blkm; n++)
											//for(n=0;n<(sqrblkn[mcid][(int)(((2*k+1)*morder[msid])/morder[mcid])][(((2*k+1)*morder[msid])%morder[mcid])+(jj)][((2*k+1)*morder[msid])%morder[mcid]-morder[msid]+(jj+ii)%morder[msid]]);n++)
											for(n=0;n<(mirror[mcid].sqrblocks[(int)(((2*k+1)*morder[msid])/morder[mcid])].blk[(((2*k+1)*morder[msid])%morder[mcid])+(jj)][((2*k+1)*morder[msid])%morder[mcid]-morder[msid]+(jj+ii)%morder[msid]].blkn);n++)
											//for(n=0;n<(mirror[mcid].sqrblocks[(int)(((2*k+1)*mirror[msid].morder)/mirror[mcid].morder)].blk[(((2*k+1)*mirror[msid].morder)%mirror[mcid].morder)+(jj)][((2*k+1)*mirror[msid].morder)%mirror[mcid].morder-mirror[msid].morder+(jj+ii)%mirror[msid].morder].blkn);n++)
											///msquare[msid][2*k+1][jj][kk][l][m] += floor(PRECP*( msquare[mcid][(int)(((2*k+1)*morder[msid])/morder[mcid])][(((2*k+1)*morder[msid])%morder[mcid])+(jj)][((2*k+1)*morder[msid])%morder[mcid]-morder[msid]+(jj+ii)%morder[msid]][l][n] * mrpartner[msid][k][(jj+ii)% morder[msid]][kk][n][m]))/PRECP;
											//temp+= floor(PRECP*( msquare[mcid][(int)(((2*k+1)*morder[msid])/morder[mcid])][(((2*k+1)*morder[msid])%morder[mcid])+(jj)][((2*k+1)*morder[msid])%morder[mcid]-morder[msid]+(jj+ii)%morder[msid]][l][n] * mrpartner[msid][k][(jj+ii)% morder[msid]][kk][n][m]))/PRECP;
											temp+= ( mirror[mcid].sqrblocks[(int)(((2*k+1)*morder[msid])/morder[mcid])].blk[(((2*k+1)*morder[msid])%morder[mcid])+(jj)][((2*k+1)*morder[msid])%morder[mcid]-morder[msid]+(jj+ii)%morder[msid]].blkelement[l][n] * mirror[msid].rpartnerblocks[k].blk[(jj+ii)% morder[msid]][kk].blkelement[n][m]);
											//mirror[msid].sqrblocks[2*k+1].blk[jj][kk].blkelement[l][m]+= floor(PRECP*(mirror[mcid].sqrblocks[(int)(((2*k+1)*mirror[msid].morder)/mirror[mcid].morder)].blk[(((2*k+1)*mirror[msid].morder)%mirror[mcid].morder)+(jj)][((2*k+1)*mirror[msid].morder)%mirror[mcid].morder-mirror[msid].morder+(jj+ii)%mirror[msid].morder].blkelement[l][n] * mirror[msid].rpartnerblocks[k].blk[(jj+ii)%mirror[msid].morder][kk].blkelement[n][m]))/PRECP;
											#pragma omp atomic update
											mirror[msid].sqrblocks[2*k+1].blk[jj][kk].blkelement[l][m]+=temp;
											//msquare[msid][2*k+1][jj][kk][l][m]+=temp;
										}
								}
							}
						}
					}
				} 
			} //end of j!=0 condition 
			
			i++;
		} while(i<noofloops);
		
		//wtime = omp_get_wtime() - wtime;
		//printf( "Time taken by thread %d is %f\n", omp_get_thread_num(), wtime );
		//End of operation		
	}
	//printer(order, mata, inverta, (int)(sizeof(mirror)/sizeof(mirror[0])), mirror);

	for(i=0;i<mirrorsize;i++)
	{
		for(j=0; j<mirror[i].mblocksize; j++)
		{
			for(k=0; k< mirror[i].sqrblocks[j].brows;k++)
			{
				for(l=0;l<mirror[i].sqrblocks[j].bcols;l++)
				{
					for(m=0;m<(mirror[i].sqrblocks[j].blk[k][l].blkm);m++) free(mirror[i].sqrblocks[j].blk[k][l].blkelement[m]);
					free(mirror[i].sqrblocks[j].blk[k][l].blkelement);
				}	
				free(mirror[i].sqrblocks[j].blk[k]);
			}
			free(mirror[i].sqrblocks[j].blk);
		}
		free(mirror[i].sqrblocks);
		for(j=0; j<(mirror[i].mblocksize/2); j++)
		{
			for(k=0; k< mirror[i].rpartnerblocks[j].brows;k++)
			{
				for(l=0;l<mirror[i].rpartnerblocks[j].bcols;l++)
				{
					for(m=0;m<(mirror[i].rpartnerblocks[j].blk[k][l].blkm);m++) free(mirror[i].rpartnerblocks[j].blk[k][l].blkelement[m]);
					free(mirror[i].rpartnerblocks[j].blk[k][l].blkelement);
				}	
				free(mirror[i].rpartnerblocks[j].blk[k]);
			}
			free(mirror[i].rpartnerblocks[j].blk);
		}
		free(mirror[i].rpartnerblocks);
		for(j=0; j<(mirror[i].mblocksize/2); j++)
		{
			for(k=0; k< mirror[i].lpartnerblocks[j].brows;k++)
			{
				for(l=0;l<mirror[i].lpartnerblocks[j].bcols;l++)
				{
					for(m=0;m<(mirror[i].lpartnerblocks[j].blk[k][l].blkm);m++) free(mirror[i].lpartnerblocks[j].blk[k][l].blkelement[m]);
					free(mirror[i].lpartnerblocks[j].blk[k][l].blkelement);
				}	
				free(mirror[i].lpartnerblocks[j].blk[k]);
			}
			free(mirror[i].lpartnerblocks[j].blk);
		}
		free(mirror[i].lpartnerblocks);
	}
	free(mirror);

	for(i=0;i<noofloops;i++) free(loopid[i]);

	free(loopid);
	free(blocks);
	free(blockspos);
	return 1;
}



