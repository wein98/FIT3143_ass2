#include "mpi.h"
#include <time.h>
#include <pthread.h> 
#include <stdbool.h>

#define NUM_THREADS 2

int infraredReadings[100];	// shared array for infrared readings
int node_size;			// shared sensor nodes size

void *InfraredSim(); // infrared simulation function format
bool Comparison(int row, int col, int recv_rank, int recv_iteration);

int base_station(MPI_Comm comm_world, MPI_Comm comm, int row, int col, int _node_size) {
	
	pthread_t tid[NUM_THREADS];
	int threadNum[NUM_THREADS];
	node_size = _node_size;
	
	//int count, size;
	int flag = 0;
	int iteration = 0;
	//MPI_Comm_size(comm_world, &size);
	
	int true_alerts = 0;
	int false_alerts = 0;
	
	FILE *pFile = fopen("base_station_log.txt", "a");
	
	double main_timetaken, starttime, endtime;

	MPI_Status status, status2;
	
	// create a thread for infrared satellite simulation
	threadNum[0] = 0;
	pthread_create(&tid[0], 0, InfraredSim, &threadNum[0]);
	
	starttime = MPI_Wtime();
	
	// code for base station computation
	while (iteration < 100) 
	{
		int count;
		int x = 0;
		
		while (x < node_size/2)
		{
			int flag = 0;
			struct timespec start, end;
			double time_taken;
			clock_gettime(CLOCK_MONOTONIC, &start);
		    	
		    	// listens for any possible incoming messages
		    	while(flag==0)
		    	{
		    		clock_gettime(CLOCK_MONOTONIC, &end);
		    		time_taken = (end.tv_sec - start.tv_sec) * 1e9;
					time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;

					// break from this loop if it takes longer than 100ms
		    		if (time_taken > 0.1) {
		    			break;
		    		}
		    			
					MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm_world, &flag,&status);
		    	}
		    	
		    	if(flag)
		    	{
					MPI_Get_count(&status, MPI_INT, &count);
					
					// receive and unpack message
		    		if(status.MPI_TAG == 2) {
		    			int recv_rank = status.MPI_SOURCE;
		    			int recv_temp, recv_iteration;
		    			char recv_timestmp[50];
		    			char recv_buffer[150];
		    			double recv_ts, local_ts, time_taken;
		    			int position = 0;
		    			bool alert = false;
		    			
						MPI_Recv(recv_buffer, 150, MPI_PACKED, status.MPI_SOURCE, 2, comm_world, &status2);
		    			local_ts = MPI_Wtime(); 
		    			
		    			MPI_Unpack(recv_buffer, 150, &position, &recv_temp, 1, MPI_INT, comm_world);
		    			MPI_Unpack(recv_buffer, 150, &position, &recv_ts, 1, MPI_DOUBLE, comm_world);
		    			MPI_Unpack(recv_buffer, 150, &position, &recv_iteration, 1, MPI_INT, comm_world);	    
		    			MPI_Unpack(recv_buffer, 150, &position, &recv_timestmp, 50, MPI_CHAR, comm_world);
		    			
		    			time_taken = local_ts - recv_ts;
		    		    			
		    			// Do comparison between received alert with infrared readings
		    			alert = Comparison(row, col, recv_rank, recv_iteration);
		    			
						time_t now;
	                	struct tm ts;
						char loggedTime[50];
	                
	                	time(&now);
	                
	                	ts = *localtime(&now);
	                	strftime(loggedTime, sizeof(loggedTime), "%a %Y-%m-%d %H:%M:%S", &ts);
	    				// TODO: Log to file

						int infraredReadingPointer = (recv_iteration%50)*2; 
						int infraVal = infraredReadings[infraredReadingPointer];
						int infraRank = infraredReadings[infraredReadingPointer+1];
						
	    				fprintf(pFile, "------------------------------------------------------------------\n");
	    				fprintf(pFile, "Iteration: %d\n", recv_iteration);
	    				fprintf(pFile, "Alert reported time: \t%s\n", recv_timestmp);
	    				fprintf(pFile, "Logged time: \t\t\t%s\n", loggedTime);

		    			if (alert) 
		    			{
		    				fprintf(pFile, "Alert Type: True\n\n");
		    				true_alerts += 1;
		    			} else 
		    			{
		    				fprintf(pFile, "Alert Type: False\n\n");
		    				false_alerts += 1;
		    			}
	    				fprintf(pFile, "Activated node \nCoord: (%d, %d) \nTemp: %d \n\n", recv_rank/col, recv_rank%col, recv_temp);
	    				fprintf(pFile, "Infrared Satellite \n Reported Coord: (%d, %d) \nReported temp: %d\n", infraRank/col, infraRank%col, infraVal);
	    				fprintf(pFile, "Communication time (seconds): %3f\n", time_taken);
	    				fprintf(pFile, "------------------------------------------------------------------\n");
		    		}
		    	}
		    	x += 1;
		}

		sleep(2);
		iteration += 1;
	}
	
	fclose(pFile);
	
	// Calculate total time taken of the life of the program
	endtime = MPI_Wtime();
	main_timetaken = endtime - starttime;
	
	printf("Total true alerts: %d, false alerts: %d.\n", true_alerts, false_alerts);
	printf("Total time taken: %f.\n", main_timetaken);
	return 0;
}

void *InfraredSim() 
{
	for(int i=0; i<100; i+=2) 
	{
		int cyclicPointer = i%50;
		unsigned int seed = time(NULL)*i;
    	infraredReadings[cyclicPointer] = (rand_r(&seed) % (90 - 70 + 1)) + 70;
    	//infraredReadings[i] = 85;
    	infraredReadings[cyclicPointer+1] = rand_r(&seed) % (node_size);
    	//printf("INFRARED. Iteration: %d, value: %d, rank: %d.\n", i, infraredReadings[i], infraredReadings[i+1]);  
    	sleep(2);
	}
}

bool Comparison(int row, int col, int recv_rank, int recv_iteration) 
{
	int adjacent[8] = {0, 1, 0, -1, 1, 0, -1, 0};
	int recvCoor[2], infraCoor[2];	
	int x, y;
	int infraredReadingPointer = (recv_iteration%50)*2;		
	
	// retrieve simulated reading from infrared
	int infraVal = infraredReadings[infraredReadingPointer];
	int infraRank = infraredReadings[infraredReadingPointer+1];
		
	// received rank coord
	recvCoor[0] = recv_rank/col;
	recvCoor[1] = recv_rank%col;
	
	// infrared rank coord
	infraCoor[0] = infraRank/col;
	infraCoor[1] = infraRank%col;
	
	// if infrared value is > 80
	if (infraVal > 80) 
	{
		// returns true if received rank is same as infrared rank
		if (recv_rank == infraRank) 
		{
			return true;
		}
		
		// otherwise checks if infrared neighbourhood matches received rank
		for (int i=0; i<50; i+=2) 
		{
		    x = recvCoor[0] + adjacent[i];
		    y = recvCoor[1] + adjacent[i+1];
		    
		    if (x >= 0 && x < row && y >= 0 && y < col) 
		    	if (x == infraCoor[0] && y == infraCoor[1])
		    		return true;
		}
	}
	
	return false;
}
