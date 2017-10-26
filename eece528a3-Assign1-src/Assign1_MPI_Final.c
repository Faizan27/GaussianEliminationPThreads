#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <stdbool.h>
#include <sys/time.h>
 
void create2DMatrix(float ***a,int row,int col);

//int processCount; //number of tasks
//int taskid; //tasks identifier

int main(void) { 
	MPI_Init(NULL, NULL);
	int processCount; //number of tasks 
	int taskid; //tasks identifier
	//MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
        MPI_Comm_size(MPI_COMM_WORLD, &processCount);  
	//double first, last; //a variable to calculate elapsed time  
	int i, j, k,i1,i2;	//General variable to make operations 
	int left_overs,localNum_of_Rows; //left_overs is the piece of code that is not covered by any processor as code/processor counts left_overs.
	//localNum_of_Rows keeps count of rows each processor will execute or work with
	float **cooff,**localMat,**refRow,dividingFactor;
	clock_t start,stop;
	int Num_of_Rows1,Num_of_Rows;
	int Num_of_Cols1,Num_of_Cols;
	//MPI_Status status;
	//Initalize the processes task id and tasks identifier 
	//MPI_Init(&argc, &argv); 
	//MPI_Comm_rank(MPI_COMM_WORLD, &taskid); 
	//MPI_Comm_size(MPI_COMM_WORLD, &processCount); 
	
	/**************************** master task ************************************/ 
	/*Process 0 send matrix data to the slaves */ 
	if (taskid == 0) 
	{	
		printf("Enter the Size_of_Matrix : ");
		fflush(stdout);
		scanf("%d",&Num_of_Rows);
		Num_of_Cols = Num_of_Rows + 1; // For Augmented Matrix
		printf("Number of Processors used are %d\n",processCount);
        	printf("Augmented Data Matrix Size : Rows = %d and Cols = %d\n",Num_of_Rows,Num_of_Cols);
		
		//Create Augmented Matrix
		create2DMatrix(&cooff,Num_of_Rows,Num_of_Cols);
		
        	//Initalize Coefficients of Augmented A:B Matrix 		
		/*   Initialize Matrix A of AX = B */
        	for (i=0;i<=Num_of_Rows-1;i++){
                	for (j=0;j<=Num_of_Cols-1;j++) {
                        	cooff[i][j] = (float) (rand()%10);
                        	}
        	}
        	printf("************ Matrices initiated****************\n"); 
		printf("Moving into MPI World\n");
		start = clock();
		//Here is where the distribution takes place..
		// Send dimension info
		//This is done to make sure all processors know the size of matrix
		MPI_Bcast(&Num_of_Rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&Num_of_Cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
		//This is the piece of matrix that would stay untouched or undistributed
		left_overs = Num_of_Rows % processCount;
		localNum_of_Rows = Num_of_Rows / processCount + (left_overs > 0);
		printf("Process 0 is executing with %d number of local rows\n",localNum_of_Rows);

		// Distribute the rows to the workers
		create2DMatrix(&localMat,localNum_of_Rows, Num_of_Cols);		
		for (i = 0; i < Num_of_Rows / processCount; i++) {
			//Distribute the data all across all processors
			MPI_Scatter(&cooff[i * processCount][0], Num_of_Cols, MPI_FLOAT,&localMat[i][0], Num_of_Cols, MPI_FLOAT,0, MPI_COMM_WORLD);
		}
		if (left_overs > 0) {
			memcpy(&localMat[localNum_of_Rows - 1][0], &cooff[Num_of_Rows - left_overs][0], Num_of_Cols * sizeof(float));
			for ( i = 1; i < left_overs; i++) {
				MPI_Send(&cooff[Num_of_Rows - left_overs + i][0], Num_of_Cols, MPI_FLOAT, i % processCount, 0, MPI_COMM_WORLD);
			}
		}

		// Process the matrix
		for (i1 = 0; i1 < Num_of_Rows; i1++) {
			// Finalize the reference row - Set it to 0
			memset(&cooff[i1][0], 0, i1 * sizeof(float));
			if (i1 % processCount == 0) {
				memcpy(&cooff[i1][i1], &localMat[i1 / processCount][i1], (Num_of_Cols - i1) * sizeof(float));
			}
			// Send or receive the finalized reference row
			MPI_Bcast(&cooff[i1][i1], Num_of_Cols - i1, MPI_FLOAT, i1 % processCount, MPI_COMM_WORLD);
			// Update the local rows
			for (i2 = (i1 + 1) / processCount + ((i1 + 1) % processCount > 0);i2 < localNum_of_Rows;i2++) 
			{
				dividingFactor = localMat[i2][i1] / cooff[i1][i1];
				for (j = i1 + 1; j < Num_of_Cols; j++) {
					localMat[i2][j] -= dividingFactor * cooff[i1][j];
				}
			}
		}
		//*******************************************************************************/
		stop = clock();
		printf("Time taken for Parallel Execution : %f seconds\n",((double) (stop - start)) / CLOCKS_PER_SEC);
	}
	else
	{
		//*************************** Slave Machine ************************************************/
		// Receive dimension info
		MPI_Bcast(&Num_of_Rows1, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&Num_of_Cols1, 1, MPI_INT, 0, MPI_COMM_WORLD);
		left_overs = Num_of_Rows1 % processCount;
		localNum_of_Rows = Num_of_Rows1 / processCount + (left_overs > taskid);
		printf("Process %d is executing with %d number of local rows\n",taskid,localNum_of_Rows);
		// Receive the rows
		create2DMatrix(&localMat,localNum_of_Rows, Num_of_Cols1);
		for (i = 0; i < Num_of_Rows1 / processCount; i++) {
			MPI_Scatter(NULL, Num_of_Cols1, MPI_FLOAT,&localMat[i][0], Num_of_Cols1, MPI_FLOAT,0, MPI_COMM_WORLD);
		}
		if (left_overs > taskid) {
			MPI_Recv(&localMat[localNum_of_Rows - 1][0], Num_of_Cols1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		//create a referenceMatrix of Rowsize 1 for local access
		create2DMatrix(&refRow,1, Num_of_Cols1);
		// Process the matrix
		for (i1 = 0; i1 < Num_of_Rows1; i1++) {
			// Finalize the reference row
			if (i1 % processCount == taskid) {
				memcpy(&refRow[0][i1], &localMat[i1 / processCount][i1], (Num_of_Cols1 - i1) * sizeof(float));
			}
			// Send or receive the finalized reference row
			MPI_Bcast(&refRow[0][i1], Num_of_Cols1 - i1, MPI_FLOAT, i1 % processCount, MPI_COMM_WORLD);
			// Update the local rows
			for (i2 = (i1 + 1) / processCount + ((i1 + 1) % processCount > taskid);i2 < localNum_of_Rows;i2++) {
				dividingFactor = localMat[i2][i1] / refRow[0][i1];
				for (j = i1 + 1; j < Num_of_Cols1; j++) {
					localMat[i2][j] -= dividingFactor * refRow[0][j];
				}
			}
		}
	}
	MPI_Finalize();
    return 0;
}

void create2DMatrix(float ***a, int row, int col){
	int i;
	*a = calloc(row, sizeof(float *));
    for(i = 0; i<=row-1; ++i) {
        (*a)[i] =  calloc(col, sizeof(float));        
    }
	
}

