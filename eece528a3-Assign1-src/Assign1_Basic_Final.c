//This code is the base sequential code for Guassian Elimination Using Pivoting
// Name - Faizan Mansuri
// Student Number - 98753163
//****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define debugger 0
/*   Function Declaration  */

int  FunctionMatrixPrint (int num_of_rows, int num_of_cols, float **cooff);
int  FunctionVectorPrint (int num_of_rows, float *consts);
void Base_Guassian_Elimination(float **cooff, float *values, float *consts, int n);

int main (void)
{
	float  **cooff,*values,*consts;
	int    i,j,n;

	printf(" Enter the dimension n for A[n][n] and b[n]: ");
	scanf("%d",&n);

	cooff = calloc(n, sizeof(float *)); 
	--cooff;

	for(i = 1; i<=n; ++i) {
		cooff[i] =  calloc(n, sizeof(float)); 
		--cooff[i];
	}

	values = calloc(n, sizeof(float)); 
	--values;
	consts = calloc(n, sizeof(float)); 
	--consts;

	/*   Initialize Matrix A of AX = B */

	for (i=1;i<=n;i++){
		for (j=1;j<=n;j++) {
			cooff[i][j] = (float) (rand()%10);
			}
	}

	/*  Initialize the elements of B of AX = B */

	for (i=1;i<=n;i++){
	   values[i] = (float) (rand()%10);
	}

	printf("\nMatrices initiated\n");
	if(debugger){
        	printf("\nMatrix A\n\n");
        	FunctionMatrixPrint (n,n,cooff);

        	printf("\nVector b\n\n");
       		FunctionVectorPrint (n,values);
	}
/* Call the Gaussian elimination function */

        Base_Guassian_Elimination(cooff,values,consts,n);
	if(debugger){
        	printf("\nSolution x\n\n");
        	FunctionVectorPrint (n,consts);
	}
        return(0);
}

void Base_Guassian_Elimination(float **cooff, float *values, float *consts, int n)
{

        int   i,j,k,m;
        float max_val_pivot,temp,temp1,temp2,temp3,temp4,max_val_row;

//Start Timing of Code here //
        clock_t start = clock();
/*
_______________________________________

  Forward Elimination Code is Implemented below
_______________________________________

*/

for (k=1; k<=n-1; ++k) {


     max_val_row = (float) fabs(cooff[k][k]) ;
     m = k;
     for (i=k+1; i<=n; i++){   /* Find the row with largest pivot */
               max_val_pivot = (float) fabs(cooff[i][k]);
               if(max_val_pivot > max_val_row) {max_val_row = max_val_pivot; m=i;}
     }
     if(m != k) {  /* Row interchanges */
                 temp1 = values[k];
                 values[k]  = values[m];
                 values[m]  = temp1;
                 for(j=k; j<=n; j+=5) {
                       temp = cooff[k][j];
                       cooff[k][j] = cooff[m][j];
                       cooff[m][j] = temp;
                       temp1 = cooff[k][j+1];
                       cooff[k][j+1] = cooff[m][j+1];
                       cooff[m][j+1] = temp1;
                       temp2 = cooff[k][j+2];
                       cooff[k][j+2] = cooff[m][j+2];
                       cooff[m][j+2] = temp2;
                       temp3 = cooff[k][j+3];
                       cooff[k][j+3] = cooff[m][j+3];
                       cooff[m][j+3] = temp3;
                       temp4 = cooff[k][j+4];
                       cooff[k][j+4] = cooff[m][j+4];
                       cooff[m][j+4] = temp4;
                 }
      }
        for (i=k+1; i<=n; ++i) {
          max_val_pivot = cooff[i][k]/cooff[k][k];

               for (j=k+1; j<=n; ++j) {
                   cooff[i][j] = cooff[i][j]-max_val_pivot*cooff[k][j];
               }
          values[i] = values[i]-max_val_pivot*values[k];
       }

}

        clock_t stop = clock();
        printf("\n\nTime Elapsed for Forward Elimination of %d elements : %f seconds\n",n,(double)(stop-start)/CLOCKS_PER_SEC);

/*
_______________________________________

  Back Substitution Block
_______________________________________

*/

        for (j=1; j<=n; ++j) {
                k=n-j+1;
                consts[k] = values[k];
                for(i=k+1; i<=n; ++i) {
                        consts[k] = consts[k]-cooff[k][i]*consts[i];
                }
                consts[k] = consts[k]/cooff[k][k];
        }

}

int FunctionMatrixPrint (int num_of_rows, int num_of_cols, float **cooff)
{
        int i,j;

        if ( num_of_rows <= 0 ) return (-1);
        if ( num_of_cols <= 0 ) return (-2);

        for (i = 1; i <= num_of_rows; i++) {

                for (j = 1; j <= num_of_cols; j++) {
                        printf ("%9.4f  ", cooff[i][j]);
                }

                printf("\n"); 
        }
        return (0);
}


int FunctionVectorPrint (int num_of_rows, float *consts)
{
        int i;

        if ( num_of_rows <= 0 ) return (-1);

        for (i = 1; i <= num_of_rows; i++) {
                printf ("%9.4f  \n", consts[i]);
        }
        printf("\n");  
        return (0);
}

