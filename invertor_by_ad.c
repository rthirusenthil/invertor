// General inversion function for matrix inversion using blockwise inversion methodology.  
// The program is meant to be adopted with user programs as per requirement.  
// A sample program `testinvertor.c' can call this `invertor_by_ad.c' function, with random genrated matrix of different orders.

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

//int invertblocks(int n, double** mata , double** inverta);

int invertinplace(int order, double** mat, int pos);
int inplaceblocksbya(int order, double** mat, int pos);
int inplaceblocksbyd(int order, double** mat, int pos);
int inplaceleftmatmul(double** mat, int ordera, int aposmn, int nb, int bposm, int bposn);
int inplacerightmatmul(double** mat, int orderb, int bposmn, int ma, int aposm, int aposn);
int schurcomplement(double** mat, int order, int matpos, int xposm, int xposn, int xn, int yposm, int yposn, int ym);

int invertbyaandd(int order, double** mat, double** invertmat, int pos);
int invertblockaandd(int order, double** mat, double** invertmat, int pos);
int schurcompforad(double** mat, int order, int matpos, int xposm, int xposn, int xn, int yposm, int yposn, int ym);
int schurad(double** mat, double** invertmat, int order, int matpos, int xposm, int xposn, int xn, int yposm, int yposn, int ym);
int invertmat(int n, double** mata, double** inverta)
{
	int invertstatus;
	int i,j, order;
	
	if(n<=0) return 0;
	
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			inverta[i][j]=mata[i][j];
	order = n;
	
	//order=2;
//	invertstatus = invertinplace(order, inverta, 0);
	invertstatus = invertbyaandd(order, mata, inverta, 0);
	return invertstatus;
}
int invertbyaandd(int order, double** mat, double** invertmat, int pos)
{
	int invertstatus;
	double modmat, a11, a12, a13, a21, a22, a23, a31, a32, a33;
	
	switch(order)
	{
		case 1:
			if (mat[pos][pos]==0) return 0;
			else 
				{
					invertmat[pos][pos]/=(invertmat[pos][pos]*invertmat[pos][pos]);
					invertstatus=1;
				}
			break;
		case 2:
		
			modmat=(-(invertmat[pos+0][pos+1]*invertmat[pos+1][pos+0]) + invertmat[pos+0][pos+0]*invertmat[pos+1][pos+1]);
			if(modmat==0) return 0;
	
			invertmat[pos+1][pos+1]*=invertmat[pos+0][pos+0]/modmat;
			
			invertmat[pos+0][pos+0]=(invertmat[pos+1][pos+1]/invertmat[pos+0][pos+0]);
			invertmat[pos+0][pos+1]/=(-1*modmat);
			invertmat[pos+1][pos+0]/=(-1*modmat);
			invertmat[pos+1][pos+1]/=(invertmat[pos+0][pos+0]*modmat);
			invertstatus=1;
			break;
		case 3:
			a11=invertmat[pos+0][pos+0];
			a12=invertmat[pos+0][pos+1];
			a13=invertmat[pos+0][pos+2];
			a21=invertmat[pos+1][pos+0];
			a22=invertmat[pos+1][pos+1];
			a23=invertmat[pos+1][pos+2];
			a31=invertmat[pos+2][pos+0];
			a32=invertmat[pos+2][pos+1];
			a33=invertmat[pos+2][pos+2];
			
			modmat=(-(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33);
			if(modmat==0) return 0;
	
			invertmat[pos+0][pos+0]=(-(a23*a32) + a22*a33)/modmat;
			invertmat[pos+0][pos+1]=(a13*a32 - a12*a33)/modmat;			
			invertmat[pos+0][pos+2]=(-(a13*a22) + a12*a23)/modmat;
			invertmat[pos+1][pos+0]=(a23*a31 - a21*a33)/modmat;
			invertmat[pos+1][pos+1]=(-(a13*a31) + a11*a33)/modmat;
			invertmat[pos+1][pos+2]=(a13*a21 - a11*a23)/modmat;
			invertmat[pos+2][pos+0]=(-(a22*a31) + a21*a32)/modmat;
			invertmat[pos+2][pos+1]=(a12*a31 - a11*a32)/modmat;
			invertmat[pos+2][pos+2]=(-(a12*a21) + a11*a22)/modmat; 

			break;
		default:
			//printf("\ncalling block function");
			invertstatus=invertblockaandd(order, mat, invertmat, pos);
			break;
	}
	return invertstatus;
}

int invertblockaandd(int order, double** mat, double** invertmat, int pos)
{
	//By simultaneous inverse of A and D:
	
	int invertstatus;

	int i,j;
	int ordera, orderd, mb, nb, mc, nc;
	int bposm, bposn, cposm, cposn;
	
	//step-1: Preparing the blocks A, B, C, D
	ordera=order/2;
	orderd=order-ordera;
	
	mb=ordera;
	nb=orderd;
	
	mc=orderd;
	nc=ordera;
	
	bposm=pos;
	bposn=pos+ordera;
	cposm=pos+ordera;
	cposn=pos;
	
	//step-2: Calculating A^-1 and D^-1
	invertstatus=invertbyaandd(ordera, mat, invertmat, pos);
	invertstatus=invertbyaandd(orderd, mat, invertmat, pos+ordera);

	
	//step-3: Calculating -1*A^-1*B and -1*D^-1*C
//int inplaceleftmatmul(double** mat, int ordera, int aposmn, int nb, int bposm, int bposn)
	invertstatus=inplaceleftmatmul(invertmat, ordera, pos, nb, bposm, bposn);
	invertstatus=inplaceleftmatmul(invertmat, orderd, pos+ordera, nc, cposm, cposn);

	//step-4: Calculating A^-1B * D^-1C for A and D^-1C * A^-1B for D
	
	
	//step-4: Replacing A and D in place of A^-1 and D^-1
//	for(i=0;i<ordera;i++) for(j=0;j<ordera;j++) invertmat[pos+i][pos+j]=mat[pos+i][pos+j];
//	for(i=0;i<orderd;i++) for(j=0;j<orderd;j++) invertmat[pos+ordera+i][pos+ordera+j]=mat[pos+ordera+i][pos+ordera+j];


	//step-5 Calculating Schur Complements
//schurcompforad(double** mat, int order, int matpos, int xposm, int xposn, int xn, int yposm, int yposn, int ym)
	//schurcompforad(invertmat, ordera, pos, bposm, bposn, nb, cposm, cposn, mc);
	//schurcompforad(invertmat, orderd, pos+ordera, cposm, cposn, nc, bposm, bposn, mb);
	schurad(mat,invertmat, ordera, pos, bposm, bposn, nb, cposm, cposn, mc);
	schurad(mat,invertmat, orderd, pos+ordera, cposm, cposn, nc, bposm, bposn, mb);
	
	
	//step-6 Calculating Inverse of Schur Complements
	///invertstatus=invertbyaandd(ordera, mat, invertmat, pos);
	///invertstatus=invertbyaandd(orderd, mat, invertmat, pos+ordera);
	
	invertstatus=invertinplace(ordera,  invertmat, pos);
	invertstatus=invertinplace(orderd,  invertmat, pos+ordera);

	//step-7 Multiplying location of B and C with inverted schur complements at locations D and A.
	invertstatus=inplacerightmatmul(invertmat, orderd, pos+ordera, mb, bposm, bposn);
	invertstatus=inplacerightmatmul(invertmat, ordera, pos, mc, cposm, cposn); 	
	
	//printf("\ninverse computed =======\n");	
	return invertstatus;
}

int schurad(double** mat, double** invertmat, int order, int matpos, int xposm, int xposn, int xn, int yposm, int yposn, int ym)
{
	int i,j,k;
	double *temp;
	
	temp = (double *) malloc(order*sizeof(double));
	
	for(i=0;i<order;i++)
		for(j=0;j<order;j++)
		{
			for(invertmat[matpos+i][matpos+j]=k=0;k<xn;k++)  //xn = ym
			invertmat[matpos+i][matpos+j]+=invertmat[xposm+i][xposn+k]*invertmat[yposm+k][yposn+j];
		}
		
	for(i=0;i<order;i++)
	{
		for(j=0;j<order;j++) temp[j]=invertmat[matpos+j][matpos+i];
		
		
		for(j=0;j<order;j++)
		{
			for(invertmat[matpos+j][matpos+i]=k=0;k<order;k++)
			{
				invertmat[matpos+j][matpos+i]+=mat[matpos+j][matpos+k] * temp[k];//invertmat[matpos+k][matpos+j];
			}
		}
	}

	for(i=0;i<order;i++)
		for(j=0;j<order;j++)
			invertmat[matpos+i][matpos+j]=(mat[matpos+i][matpos+j]-invertmat[matpos+i][matpos+j]);

	free(temp);
	return 1;
}

int schurcompforad(double** mat, int order, int matpos, int xposm, int xposn, int xn, int yposm, int yposn, int ym)
{
	//At the place A: The term will be A - A * (A^-1B) * (D^-1C)
	//At the place D: The term will be D - D * (D^-1C) * (A^-1B)
	
	int i,j,k, im, in;

	double *temp, *atemp, mulres;
	temp =(double *) malloc(order*sizeof(double));
	atemp=(double *)malloc(order*sizeof(double));
	
	//for(j=0;j<order; j++) temp[j]=0;
	
	for(im=0;im<order;im++)
	{
		for(in=0;in<order;in++) atemp[in]=mat[matpos+in][matpos+im];
		
		for(in=0;in<order;in++)
		{
			//for(i=0;i<order;i++) //in is the column id to take the ith column of X * Y.
			//{
				for(j=0;j<order;j++)
				{
					temp[j]=0;
					for(k=0; k<xn ; k++)  //xn = ym
					{
						temp[j]+=mat[xposm+j][xposn+k] * mat[yposm+k][yposn+in] ;
					}
					// Now the temp array contains the i'th column in the multiplication x * y;
					//if(in==j) temp[j]=1-temp[j];
					//else 
					//temp[j]*=-1;
				}
				//Now the temp array corresponds to ith column of (I - A^-1B * D^-1C)   
			//}
			mulres=0;
		
			for(k=0;k<order;k++) 
			{	
				mulres+=atemp[k]*temp[k];
			}
			mat[matpos+im][matpos+in]-=mulres;
		}
	}

	free(temp);
	free(atemp);
	
	return 1;
}

int invertinplace(int order, double** mat, int pos)
{
	int invertstatus;
	double modmat, a11, a12, a13, a21, a22, a23, a31, a32, a33;
	
	switch(order)
	{
		case 1:
			if (mat[pos][pos]==0) return 0;
			else 
				{
					mat[pos][pos]/=(mat[pos][pos]*mat[pos][pos]);
					invertstatus=1;
				}
			break;
		case 2:
		
			modmat=(-(mat[pos+0][pos+1]*mat[pos+1][pos+0]) + mat[pos+0][pos+0]*mat[pos+1][pos+1]);
			if(modmat==0) return 0;
	
			mat[pos+1][pos+1]*=mat[pos+0][pos+0]/modmat;
			
			mat[pos+0][pos+0]=(mat[pos+1][pos+1]/mat[pos+0][pos+0]);
			mat[pos+0][pos+1]/=(-1*modmat);
			mat[pos+1][pos+0]/=(-1*modmat);
			mat[pos+1][pos+1]/=(mat[pos+0][pos+0]*modmat);
			invertstatus=1;
			break;
		case 3:
			a11=mat[pos+0][pos+0];
			a12=mat[pos+0][pos+1];
			a13=mat[pos+0][pos+2];
			a21=mat[pos+1][pos+0];
			a22=mat[pos+1][pos+1];
			a23=mat[pos+1][pos+2];
			a31=mat[pos+2][pos+0];
			a32=mat[pos+2][pos+1];
			a33=mat[pos+2][pos+2];
			
			modmat=(-(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33);
			if(modmat==0) return 0;
	
			mat[pos+0][pos+0]=(-(a23*a32) + a22*a33)/modmat;
			mat[pos+0][pos+1]=(a13*a32 - a12*a33)/modmat;			
			mat[pos+0][pos+2]=(-(a13*a22) + a12*a23)/modmat;
			mat[pos+1][pos+0]=(a23*a31 - a21*a33)/modmat;
			mat[pos+1][pos+1]=(-(a13*a31) + a11*a33)/modmat;
			mat[pos+1][pos+2]=(a13*a21 - a11*a23)/modmat;
			mat[pos+2][pos+0]=(-(a22*a31) + a21*a32)/modmat;
			mat[pos+2][pos+1]=(a12*a31 - a11*a32)/modmat;
			mat[pos+2][pos+2]=(-(a12*a21) + a11*a22)/modmat; 

			break;
		default:
			//printf("\ncalling block function");
			invertstatus=inplaceblocksbyd(order, mat, pos);
			break;
	}
	return invertstatus;
}
int inplaceblocksbya(int order, double** mat, int pos)
{
	int invertstatus;

	int i,j;
	int ordera, orderd, mb, nb, mc, nc;
	int bposm, bposn, cposm, cposn;
	
	//step-1: Preparing the blocks A, B, C, D
	ordera=order/2;
	orderd=order-ordera;
	
	mb=ordera;
	nb=orderd;
	
	mc=orderd;
	nc=ordera;
	
	bposm=pos;
	bposn=pos+ordera;
	cposm=pos+ordera;
	cposn=pos;
	
	//step-2: Calculating A^-1
	invertstatus=invertinplace(ordera, mat, pos);
	
	//step-3: Calculating -1*A^-1*B
//int inplaceleftmatmul(double** mat, int ordera, int aposmn, int nb, int bposm, int bposn)
	invertstatus=inplaceleftmatmul(mat, ordera, pos, nb, bposm, bposn);
	
	//step-4: Calculating Schur complement S = D - C A^-1B
//int schurcomplement(double** mat, int order, int matpos, int xposm, int xposn, int xn, int yposm, int yposn, int ym)
	invertstatus=schurcomplement(mat, orderd, pos+ordera, cposm, cposn, nc, bposm, bposn, mb);
	
	//step-5: Calculating C * A^-1
//int inplacerightmatmul(double** mat, int orderb, int bposmn, int ma, int aposm, int aposn)
	invertstatus=inplacerightmatmul(mat, ordera, pos, mc, cposm, cposn);
	
	//step-6: Calculating S^-1 
	invertstatus=invertinplace(orderd, mat, pos+ordera);
	
	//step-7: Calculatin S^-1 * CA^-1
	invertstatus=inplaceleftmatmul(mat, orderd, pos+ordera, nc, cposm, cposn);
	
	//step-8: Calculating Schur completed at A: A^-1 + A^-1B * S^-1CA^-1
	invertstatus=schurcomplement(mat, ordera, pos, bposm, bposn, nb, cposm, cposn, mc);
	
	//step-9: Calculating -A^-1B * S^-1
	invertstatus=inplacerightmatmul(mat, orderd, pos+ordera, mb, bposm, bposn);

	return invertstatus;
}

int inplaceblocksbyd(int order, double** mat, int pos)
{
	//This function is based on invertability of d
	int invertstatus;

	int i,j;
	int ordera, orderd, mb, nb, mc, nc;
	int bposm, bposn, cposm, cposn;
	
	//step-1: Preparing the blocks A, B, C, D
	ordera=order/2;
	orderd=order-ordera;
	
	mb=ordera;
	nb=orderd;
	
	mc=orderd;
	nc=ordera;
	
	bposm=pos;
	bposn=pos+ordera;
	cposm=pos+ordera;
	cposn=pos;
	
	//step-2: Calculating D^-1
	invertstatus=invertinplace(orderd, mat, pos+ordera);
	
	//step-3: Calculating -1*D^-1*C
//int inplaceleftmatmul(double** mat, int ordera, int aposmn, int nb, int bposm, int bposn)
	invertstatus=inplaceleftmatmul(mat, orderd, pos+ordera, nc, cposm, cposn);
	
	//step-4: Calculating Schur complement S = D - B * D^-1C
//int schurcomplement(double** mat, int order, int matpos, int xposm, int xposn, int xn, int yposm, int yposn, int ym)
	invertstatus=schurcomplement(mat, ordera, pos, bposm, bposn, nb, cposm, cposn, mc);
	
	//step-5: Calculating B * D^-1
//int inplacerightmatmul(double** mat, int orderb, int bposmn, int ma, int aposm, int aposn)
	invertstatus=inplacerightmatmul(mat, orderd, pos+ordera, mb, bposm, bposn);
	
	//step-6: Calculating S^-1 
	invertstatus=invertinplace(ordera, mat, pos);
	
	//step-7: Calculatin S^-1 * BD^-1
	invertstatus=inplaceleftmatmul(mat, ordera, pos, nb, bposm, bposn);
	
	//step-8: Calculating Schur completed at D: D^-1 + D^-1C * S^-1BD^-1
	invertstatus=schurcomplement(mat, orderd, pos+ordera, cposm, cposn, nc, bposm, bposn, mb);
	
	//step-9: Calculating -D^-1C * S^-1
	invertstatus=inplacerightmatmul(mat, ordera, pos, mc, cposm, cposn);

	return invertstatus;
}


int inplaceleftmatmul(double** mat, int ordera, int aposmn, int nb, int bposm, int bposn)
{
	//This computes -1*mat A * mat B and stores it in mat B.
	//mat A is square matrix of order (ordera * ordera).
	//mat B is not square matrix.  It is of order (ordera * nb).
	//The result matrix will be stored at mat B's location with order (ordera * nb).
	//Since the required multiplications for block inversion are with sqare matrix, we use minimum variables.
	
	int i,j,k;

	double *btemp;
	btemp =(double *) malloc(ordera*sizeof(double));

	for(j=0;j<nb;j++)
		{
			//printf("-------");
			for(k=0; k<ordera; k++) 
			{
				btemp[k]=mat[bposm+k][bposn+j];  //printf("\n%lf",btemp[k]);
				//mat[bposm+i][bposn+k]=0;
			}
		
		for(i=0;i<ordera;i++) 
			{ 
				mat[bposm+i][bposn+j]=0;
				for(k=0; k<ordera; k++)
				mat[bposm+i][bposn+j]+=mat[aposmn+i][aposmn+k]*btemp[k];
				mat[bposm+i][bposn+j]*=-1;
			}
		}	
	free(btemp);
	return 1;
}

int inplacerightmatmul(double** mat, int orderb, int bposmn, int ma, int aposm, int aposn)
{
	//Left Multiplication we track -ve sign.  For Right Multiplication we track +ve sign.
	//This computes mat A * mat B and stores it in mat A.
	//mat A is not square matrix.  It is of order (ma * orderb).
	//mat B is square matrix of order (orderb * orderb).
	//The result matrix will be stored at mat A's location with order (ma * orderb).
	//Since the required multiplications for block inversion are with sqare matrix, we use minimum variables.
	
	int i,j,k;

	double *atemp;
	atemp =(double *) malloc(orderb*sizeof(double));

	for(i=0;i<ma;i++)
		{
			//printf("-------");
			for(k=0; k<orderb; k++) 
			{
				atemp[k]=mat[aposm+i][aposn+k];  //printf("\n%lf",atemp[k]);
				//mat[bposm+i][bposn+k]=0;
			}
		
		for(j=0;j<orderb;j++) 
			{ 
				mat[aposm+i][aposn+j]=0;
				for(k=0; k<orderb; k++)
				mat[aposm+i][aposn+j]+=atemp[k]*mat[bposmn+k][bposmn+j];
				//mat[aposm+i][aposn+j]*=-1;
			}
		}	
	free(atemp);
	return 1;
}

int schurcomplement(double** mat, int order, int matpos, int xposm, int xposn, int xn, int yposm, int yposn, int ym)
{
	//If block D is invertible, then Schur complement of the block D is
	//	M / D := A − B D^{-1} C 
	//If block A is invertible, then Schur complement of the block A is
	//	M / A := D − C A^{-1} B
	
	//In the first case, we are using following steps.
	//First we compute D^{-1}.  Then we calculate -D^{-1}C and stores at the place of C.
	//Using already stored B and -D^{-1}C, we can calculate the schur term and stores at the place of A.
	
	//In the Second case, we are using following steps.
	//First we compute A^{-1}.  Then we calculate -A^{-1}B and stores at the place of B.
	//Using already stored -A^{-1}B and C, we can calculate the schur term and stores at the place of D.

	//This function starts by assuming that we have already calculated and stored the following:-
	//case D is invertible: -D^{-1}C is available in place of C
	//case A is invertible: -A^{-1}B is available in place of A
	
	//mat is at the location i=matpos, j=matpos with order (order*order).
	//x matrix is at the location (xposm, xposn) with order (xm * order).
	//y matrix is at the location (yposm, yposn) with order (order * yn).
	//Here we calculate, mat = mat + xmatrix*ymatrix and stores at the location of mat.
	
	int i,j,k;
	
	for(i=0;i<order;i++)
		for(j=0;j<order;j++)
			for(k=0;k<xn;k++)  //xn == ym
			{
				mat[matpos+i][matpos+j]+=mat[xposm+i][xposn+k]*mat[yposm+k][yposn+j];
			}
	return 1;
}


