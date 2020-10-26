/* Gets the neighbors in a cartesian communicator
* Orginally written by Mary Thomas
* - Updated Mar, 2015
* Link: https://edoras.sdsu.edu/~mthomas/sp17.605/lectures/MPI-Cart-Comms-and-Topos.pdf
* Modifications to fix bugs, include an async send and receive and to revise print output
*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>
#include <memory.h>

#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1
#define THRESHOLD 80
#define TOL_RANGE 5
#define REQUEST_DATA_TAG 0
#define NO_REQUEST_DATA_TAG 1

int main(int argc, char *argv[]) {

	int ndims=2, size, my_rank, reorder, my_cart_rank, ierr;
	int nrows, ncols;
	int nbr_i_lo, nbr_i_hi;
	int nbr_j_lo, nbr_j_hi;
	MPI_Comm comm2D;
	int dims[ndims],coord[ndims];
	int wrap_around[ndims];
	
	/* start up initial MPI environment */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	/* process command line arguments*/
	if (argc == 3) {
		nrows = atoi (argv[1]);
		ncols = atoi (argv[2]);
		dims[0] = nrows; /* number of rows */
		dims[1] = ncols; /* number of columns */
		if( (nrows*ncols) != size-1) {
			if( my_rank ==0) printf("ERROR: nrows*ncols)=%d * %d = %d != %d\n", nrows, ncols, nrows*ncols,size);
			MPI_Finalize(); 
			return 0;
		}
	} else {
		nrows=ncols=(int)sqrt(size-1);
		dims[0]=dims[1]=0;
	}

    /*************************************************************/
    /* create cartesian topology for processes */
    /*************************************************************/
    MPI_Dims_create(size-1, ndims, dims);
    //if(my_rank==0)
        //printf("Root Rank: %d. Comm Size: %d: Grid Dimension = [%d x %d] \n",my_rank,size,dims[0],dims[1]);
    
    /* create cartesian mapping */
    wrap_around[0] = 0;
    wrap_around[1] = 0; /* periodic shift is .false. */
    reorder = 1;
    ierr = 0;
    ierr = MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, wrap_around, reorder, &comm2D);
    if(ierr != 0) printf("ERROR[%d] creating CART\n",ierr);
	
    if (my_rank == size-1) {
        // BASE STATION CODE
        if(my_rank==size-1)
            printf("Base Station Rank: %d. Comm Size: %d: Grid Dimension = [%d x %d] \n",my_rank,size,dims[0],dims[1]);
        
    } else {
        
        /* find my coordinates in the cartesian communicator group */
        MPI_Cart_coords(comm2D, my_rank, ndims, coord); // coordinated is returned into the coord array
        /* use my cartesian coordinates to find my rank in cartesian group*/
        MPI_Cart_rank(comm2D, coord, &my_cart_rank);
        /* get my neighbors; axis is coordinate dimension of shift */
        /* axis=0 ==> shift along the rows: P[my_row-1]: P[me] : P[my_row+1] */
        /* axis=1 ==> shift along the columns P[my_col-1]: P[me] : P[my_col+1] */
        
        MPI_Cart_shift( comm2D, SHIFT_ROW, DISP, &nbr_i_lo, &nbr_i_hi );
        MPI_Cart_shift( comm2D, SHIFT_COL, DISP, &nbr_j_lo, &nbr_j_hi );
        
        // Defines an array of neighbours for each node
        int neighbours[4] = {-1};
        neighbours[0] = nbr_i_lo;
        neighbours[1] = nbr_i_hi;
        neighbours[2] = nbr_j_lo;
        neighbours[3] = nbr_j_hi;

        // printf("Global rank: %d. Cart rank: %d. Coord: (%d, %d). Left: %d. Right: %d. Top: %d. Bottom: %d\n", my_rank, my_cart_rank, coord[0], coord[1], nbr_j_lo, nbr_j_hi, nbr_i_lo, nbr_i_hi);
        //fflush(stdout);
        
        /*
        int bCastVal = -1;
        if(my_cart_rank == 4){
            bCastVal = 400;
        }
        MPI_Bcast(&bCastVal, 1, MPI_INT, 4, comm2D);
        printf("Global rank: %d. Cart rank: %d. Coord: (%d, %d). BCast Value: %d\n", my_rank, my_cart_rank, coord[0], coord[1], bCastVal);
        fflush(stdout);
        */
    
        int timestamp = 0;
        int randomVal = -1;
        int noReques = 0;

        while (timestamp < 5) {
		    MPI_Request send_request[4];
		    MPI_Request receive_request[4];
		    MPI_Request receive_request2[4];
		    MPI_Request request;
		    MPI_Status send_status[4];
		    MPI_Status receive_status[4];
		    MPI_Status status;
	        unsigned int seed = time(NULL)*my_rank;
            randomVal = rand_r(&seed) % 150 + 1;
	
            for(int i=0; i<4; i++)
            {
                if (neighbours[i] != -2) {
                    MPI_Isend(&randomVal, 1, MPI_INT, neighbours[i], timestamp, comm2D, &send_request[i]);
                }
            }
	    //printf("Rank: %d. Value: %d. Timestmp: %d.\n", my_rank, randomVal, timestamp);
            if (randomVal > THRESHOLD) {
                int recVals[4] = {0};
                // send request msg to neighbour nodes and waits for a reply of requested msg
                
                for(int i=0; i<4; i++)
                {
                    if (neighbours[i] != -2) {
                        MPI_Irecv(&recVals[i], 1, MPI_INT, neighbours[i], timestamp, comm2D, &receive_request[i]);	
                        MPI_Wait(&receive_request[i], &status);
                        //printf("RECEIVE. Rank: %d, timestmp: %d, neighbour: %d.\n", my_rank, timestamp, neighbours[i]);
                    }
                }
                
                printf("Global rank: %d. Cart rank: %d. Coord: (%d, %d). Random Val: %d. Recv Top: %d. Recv Bottom: %d. Recv Left: %d. Recv Right: %d. Timestmp: %d\n", my_rank, my_cart_rank, coord[0], coord[1], randomVal, recVals[0], recVals[1], recVals[2], recVals[3], timestamp);

                int matchCnt = 0;
                int matchCntArr[4] = {0};
                for(int i=0; i<4; i++)
                {
                    int range = recVals[i] - randomVal;
                    if (recVals[i] > THRESHOLD && abs(range) <= 20) {
                        // Stores matched data rank
                        matchCntArr[matchCnt] = i;
                        matchCnt += 1;
                    }
                }

                bool toBaseFlag = false;
                // If more than 2 neighbours match, sends report to base station
                // Informs affected neighbour this node will sends report to base
                if (matchCnt >= 2) {
                    int maxReading = 0;
                    for (int i = 0; i < matchCnt; i++) {
                        if (randomVal > recVals[matchCntArr[i]]) {
                            toBaseFlag = true;
                        } else if (recVals[matchCntArr[i]] == randomVal)
                        {
                            if (my_rank > matchCntArr[i])
                                toBaseFlag = true;
                            else
                                toBaseFlag = false;
                        } else
                        {
                            toBaseFlag = false;
                        }                       
                    }

                    if (toBaseFlag) {
                        printf("TO BASE. Rank: %d", my_rank);
                        // Sends report to BASE station
                        // 1. temp. reading
                        // 2. timestamp
                        // 3. nodes rank
                    } 
                }
                
                // MPI_Waitall(4, send_request, send_status);
                // MPI_Waitall(4, receive_request, receive_status);


            } else {


            printf("Global rank: %d. Cart rank: %d. Coord: (%d, %d). Random Val: %d. Timestmp: %d.\n", my_rank, my_cart_rank, coord[0], coord[1], randomVal, timestamp);
            }
            //printf("MID. Rank: %d. Timestmp: %d.\n", my_rank, timestamp);
		
	    //printf("END. Rank: %d. Timestmp: %d.\n", my_rank, timestamp);
            timestamp += 1;
            sleep(2);
        }
        
        
        
        
        // printf("Global rank: %d. Cart rank: %d. Coord: (%d, %d). Random Val: %d. Recv Top: %d. Recv Bottom: %d. Recv Left: %d. Recv Right: %d.\n", my_rank, my_cart_rank, coord[0], coord[1], randomVal, recVals[0], recVals[1], recVals[2], recVals[3]);
        
        // // each process writes to it's own text file according to their rank
        // char filename[30];
        // sprintf(filename, "process_%d_task2b.txt", my_rank);
        // FILE *pFile = fopen(filename, "w");
        
        // // compare receive values who its randomVal
        // flag = true;
        // for (int i=0; i<4; i++) {
        //     if (randomVal != recVals[i] && recVals[i] != 0) {
        //         flag = false;
        //         break;
        //     }
        // }
        
        // // if compare flag is true, log to process_i.txt file
        // if (flag) {
        //     fprintf(pFile, "Global rank: %d. Cart rank: %d. Coord: (%d, %d). Random Val: %d. Recv Top: %d. Recv Bottom: %d. Recv Left: %d. Recv Right: %d.\n", my_rank, my_cart_rank, coord[0], coord[1], randomVal, recVals[0], recVals[1], recVals[2], recVals[3]);
        //     printf("[MATCHED] Global rank: %d. Cart rank: %d. Coord: (%d, %d). Random Val: %d. Recv Top: %d. Recv Bottom: %d. Recv Left: %d. Recv Right: %d.\n", my_rank, my_cart_rank, coord[0], coord[1], randomVal, recVals[0], recVals[1], recVals[2], recVals[3]);
        // }
        
        // fclose(pFile);
        MPI_Comm_free( &comm2D );
    }
	MPI_Finalize();
	return 0;
}
