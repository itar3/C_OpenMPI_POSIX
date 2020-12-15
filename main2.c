//
//  main2.c
//  assignment2
//
//  Created by Ivan Tarabykin on 30/10/20.
//

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>

#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1

#define SENTINEL 100
#define MAX_ITERATION 10
#define TERMINATION_SIGNAL 4
#define SATTELITE_ARRAY 1000

void* InfraredFunc(void *pArg);
int master_io(MPI_Comm world_comm, MPI_Comm comm);
int slave_io(MPI_Comm world_comm, MPI_Comm comm);

//structre used by the thread for simulating IRS
struct Log{
    time_t timeLog;
    int node;
    int temp;
  };

//structure that used for MPI_Send in order to send potential alerts
struct alert{
    
    int node;
    int temp;
    int nbs[4];
    int nbsTemp[4];
    char coords[10];
    
};

//structure used to keep track of time spend on 1)overall program 2)time taken to send potential alerts
struct timespec commStart, commEnd, totalStart, totalEnd;

//global vars used to generate user specified grid
int agc; // variable to store argc
int cols,rows;

struct Log satellitelogs[SATTELITE_ARRAY];


int main(int argc, char *argv[])
{
    int rank, size; //variables to keep track of size = number of processes, and rank = current process
    
    //initializing MPI
    MPI_Comm new_comm;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    agc = argc;
    //assigning extra arguments to rows and cols if provided (for user specified grid)
    if (argc == 3) {
        rows = atoi(argv[1]);
        cols = atoi(argv[2]);
    }
    
    
    clock_gettime(CLOCK_MONOTONIC, &totalStart); //start taking time for the whole program
    MPI_Comm_split( MPI_COMM_WORLD,rank == size-1, 0, &new_comm);
    if (rank == size-1) //size- 1 = last process, which is used as a base station (master node)
        master_io(MPI_COMM_WORLD, new_comm );
    else                                        //all the other processes used by workers - WSN nodes
        slave_io( MPI_COMM_WORLD, new_comm );
    MPI_Finalize();
    clock_gettime(CLOCK_MONOTONIC, &totalEnd);// end taking time for the whole program
    
    //calculate the total time running and print the results out
    double time_taken = (totalEnd.tv_sec - totalStart.tv_sec)* 1e9;
    time_taken = (time_taken + (totalEnd.tv_nsec - totalEnd.tv_nsec)) * 1e-9;
    if(rank == 0){
        printf("\n\nTotal time taken %f\n", time_taken);
    }
    return 0;
}


int master_io(MPI_Comm world_comm, MPI_Comm comm)
{
    int size, nslaves;
    MPI_Comm_size(world_comm, &size );
    
    //structure and all the necessary stuff to recieve and store potential alets
    struct alert recv;
    MPI_Datatype Valuetype;
    MPI_Datatype type[5] = {MPI_INT, MPI_INT,MPI_INT,MPI_INT,MPI_CHAR};
    int blocklen[5] = {1,1,4,4,10};
    MPI_Aint disp[5];
    MPI_Status status;

    MPI_Get_address(&recv.node, &disp[0]);
    MPI_Get_address(&recv.temp, &disp[1]);
    MPI_Get_address(&recv.nbs, &disp[2]);
    MPI_Get_address(&recv.nbsTemp, &disp[3]);
    MPI_Get_address(&recv.coords, &disp[4]);
    
    //displacement of each vars from 0
    disp[4] = disp[4] - disp[0];
    disp[3] = disp[3] - disp[0];
    disp[2] = disp[2] - disp[0];
    disp[1] = disp[1] - disp[0];
    disp[0] = 0;
    
    //create struct and commit the custom data type
    MPI_Type_create_struct(5, blocklen, disp, type, &Valuetype);
    MPI_Type_commit(&Valuetype);
    
    //initialize and create the thread that runs simulated IRS
    pthread_t tid;
    pthread_create(&tid, NULL, InfraredFunc, NULL); // Create the thread
    
    
    
    int val = 0; //sentinel value checker
    FILE* input;    //file for user input
    FILE* log;      //file used for logs
    int iteration = 0;      //keep track of iterations
    int time_treshold = 5;  //time threshold between potential alert and IRS reading
    int temp_threshold = 10;    // temperature threshold between potential alert and IRS reading
    int dims[2];                //dimensions of the grid
    
    nslaves = size - 1; // number of slave
    
    MPI_Recv(dims, 8, MPI_INT, 0, 0, world_comm, MPI_STATUS_IGNORE);  //recieve the grid dimensions
    
    //main loop of the base station with MAX_ITERATION specified in the #define section
    while(iteration <  MAX_ITERATION){
        printf("[base] iteration %d\n", iteration);
        
        if (iteration == MAX_ITERATION-1){      //if on the last iteration  - send the termination signal to worker nodes in order to exit
            for(int i = 0; i< nslaves; i++){
                MPI_Send(0,0,MPI_INT, i ,TERMINATION_SIGNAL, world_comm);
                
            }
            break; //break from the while loop
        }
        //check for the user sentinel input
        input = fopen("user.txt", "r");
        fscanf(input, "%d", &val);
        fclose(input);
        
        //if user put sentinel value, and it matches the sentinel value defined above - early break from the program via sending termination signal to worker nodes
        if(val == SENTINEL){
            
            for(int i = 0; i< nslaves; i++){
                MPI_Send(0,0,MPI_INT, i ,TERMINATION_SIGNAL, world_comm);
            }
            break;
        }
        
        
        //recieve the potential alert from worker nodes and assign to the struct recv
        MPI_Recv(&recv, 1, Valuetype, MPI_ANY_SOURCE, 1, world_comm, MPI_STATUS_IGNORE);
        
        clock_gettime(CLOCK_MONOTONIC, &commEnd); //finish calculating time of communication
        
        // check the time of recieving
        time_t r_time;
        time_t recv_time = time(&r_time);
        
        //values to check for previous iteration logs
        int prev_false_iter;
        int prev_true_iter;
        int prev_true_node;
        
        //open file for logging
        log = fopen("log.txt", "a");
        
        //check for the readings of the IRS from the back (most recent reading at the back)
        for(int i = SATTELITE_ARRAY - 1 ; i>0; i--){
            if(recv.node == satellitelogs[i].node){ //if the node of the alert matches the IRS reading
                if(fabs(difftime(recv_time, satellitelogs[i].timeLog)) < time_treshold){    //if the time within threshold
                    
                    //true alert logging
                    if(abs(recv.temp - satellitelogs[i].temp) < temp_threshold || satellitelogs[i].temp > recv.temp ){ //if the temperature within threshold or satellite temp > node temp (since node is >80, and IRS is bigger than that - its a match)
                        
                        //check to prevent double logging
                        if(prev_true_iter != iteration && prev_true_node != recv.node){
                            
                            //basic printf and fprintf of a log
                            fprintf(log,"-----------------------------------------\n");
                            printf("-----------------------------------------\n");
                            
                            time_t ltime;
                            time(&ltime);
                            
                            fprintf(log, "Iteration: %d\n", iteration);
                            printf("Iteration: %d\n", iteration);
                            
                            fprintf(log, "Logged time: %s", ctime(&ltime));
                            printf("Logged time: %s", ctime(&ltime));
                            
                            fprintf(log, "Reported time: %s", ctime(&recv_time));
                            printf("Reported time: %s", ctime(&recv_time));
                            
                            fprintf(log, "Aler time: true\n\n");
                            printf("Aler time: true\n\n");
                            
                            fprintf(log,"Reporting node\t\t Coord\t\tTemp\n");
                            printf("Reporting node\t\t Coord\t\tTemp\n");
                            
                            fprintf(log,"\t%d \t\t %s \t\t %d\n", recv.node, recv.coords, recv.temp  );
                            printf("\t%d \t\t %s \t\t %d\n", recv.node, recv.coords, recv.temp  );
                            
                            int count = 0; //count number of matching neigbhours
                            fprintf(log, "\nAdjusting nodes\t\t Coord\t\tTemp\n");
                            printf("\nAdjusting nodes\t\t Coord\t\tTemp\n");
                            for (int i = 0; i<4; i++){
                                if(recv.nbs[i] != -10){
                                    //convert node number to coordinates format
                                    int x = recv.nbs[i]/dims[1];
                                    int y = recv.nbs[i]%dims[1];
                                    fprintf(log,"\t%d \t\t (%d,%d) \t\t %d\n", recv.nbs[i], x,y, recv.nbsTemp[i]);
                                    printf("\t%d \t\t (%d,%d) \t\t %d\n", recv.nbs[i], x,y, recv.nbsTemp[i]);
                                    count++;
                                }
                            }
                            //convert node number to coordinates format
                            int x = satellitelogs[i].node/dims[1];
                            int y = satellitelogs[i].node%dims[1];
                            fprintf(log,"\nInfrared Satellite Reporting Time: %s", ctime(&satellitelogs[i].timeLog));
                            printf("\nInfrared Satellite Reporting Time: %s", ctime(&satellitelogs[i].timeLog));
                            fprintf(log,"Infrared Satellite Reporting (Celsius): %d\n", satellitelogs[i].temp);
                            printf("Infrared Satellite Reporting (Celsius): %d\n", satellitelogs[i].temp);
                            fprintf(log,"Infrared Satellite Reporting Coords: (%d,%d)\n\n", x,y);
                            printf("Infrared Satellite Reporting Coords: (%d,%d)\n\n", x,y);
                            
                            //time taken of the communication
                            double time_taken = (commEnd.tv_sec - commStart.tv_sec)* 1e9;
                            time_taken = (time_taken + (commEnd.tv_nsec - commStart.tv_nsec)) * 1e-16;
                            
                            fprintf(log,"Communication Time(secods): %lf\n", time_taken);
                            printf("Communication Time(secods): %lf\n", time_taken);
                            
                            fprintf(log,"Number of adjacent matches to reporting node: %d\n", count);
                            printf("Number of adjacent matches to reporting node: %d\n", count);
                            
                            fprintf(log,"Time difference between WSN and Satellite reading: %f\n", difftime(recv_time, satellitelogs[i].timeLog));
                            printf("Time difference between WSN and Satellite reading: %f\n", difftime(recv_time, satellitelogs[i].timeLog));
                            
                            fprintf(log,"-----------------------------------------\n");
                            printf("-----------------------------------------\n");
                            
                            //prevents double logging of the same alert
                            prev_true_iter = iteration;
                            prev_true_node = recv.node;
                            
                        }
                        
                        
                        
                        
                        
                        
                        
                    //false alert logging
                    }else{
                        //prevents double logging of the same alert
                        if (prev_false_iter != iteration){
                            
                            fprintf(log,"-----------------------------------------\n");
                            printf("-----------------------------------------\n");
                            
                            time_t ltime;
                            time(&ltime);
                            
                            fprintf(log, "Iteration: %d\n", iteration);
                            printf("Iteration: %d\n", iteration);
                            
                            fprintf(log, "Logged time: %s", ctime(&ltime));
                            printf("Logged time: %s", ctime(&ltime));
                            
                            fprintf(log, "Reported time: %s", ctime(&recv_time));
                            printf("Reported time: %s", ctime(&recv_time));
                            
                            fprintf(log, "Aler time: false\n\n");
                            printf("Aler time: false\n\n");
                            
                            fprintf(log,"Reporting node\t\t Coord\t\tTemp\n");
                            printf("Reporting node\t\t Coord\t\tTemp\n");
                            
                            fprintf(log,"\t%d \t\t %s \t\t %d\n", recv.node, recv.coords, recv.temp  );
                            printf("\t%d \t\t %s \t\t %d\n", recv.node, recv.coords, recv.temp  );
                            
                            int count = 0;
                            fprintf(log, "\nAdjusting nodes\t\t Coord\t\tTemp\n");
                            printf("\nAdjusting nodes\t\t Coord\t\tTemp\n");
                            for (int i = 0; i<4; i++){
                                if(recv.nbs[i] != -10){
                                    int x = recv.nbs[i]/dims[1];
                                    int y = recv.nbs[i]%dims[1];
                                    fprintf(log,"\t%d \t\t (%d,%d) \t\t %d\n", recv.nbs[i], x,y, recv.nbsTemp[i]);
                                    printf("\t%d \t\t (%d,%d) \t\t %d\n", recv.nbs[i], x,y, recv.nbsTemp[i]);
                                    count++;
                                }
                            }
                            
                            int x = satellitelogs[i].node/dims[1];
                            int y = satellitelogs[i].node%dims[1];
                            fprintf(log,"\nInfrared Satellite Reporting Time: %s", ctime(&satellitelogs[i].timeLog));
                            printf("\nInfrared Satellite Reporting Time: %s", ctime(&satellitelogs[i].timeLog));
                            fprintf(log,"Infrared Satellite Reporting (Celsius): %d\n", satellitelogs[i].temp);
                            printf("Infrared Satellite Reporting (Celsius): %d\n", satellitelogs[i].temp);
                            fprintf(log,"Infrared Satellite Reporting Coords: (%d,%d)\n\n", x,y);
                            printf("Infrared Satellite Reporting Coords: (%d,%d)\n\n", x,y);
                            
                            double time_taken = (commEnd.tv_sec - commStart.tv_sec)* 1e9;
                            time_taken = (time_taken + (commEnd.tv_nsec - commStart.tv_nsec)) * 1e-16;
                            
                            fprintf(log,"Communication Time(secods): %lf\n", time_taken);
                            printf("Communication Time(secods): %lf\n", time_taken);
                            
                            fprintf(log,"Number of adjacent matches to reporting node: %d\n", count);
                            printf("Number of adjacent matches to reporting node: %d\n", count);
                            
                            fprintf(log,"Time difference between WSN and Satellite reading: %f\n", difftime(recv_time, satellitelogs[i].timeLog));
                            printf("Time difference between WSN and Satellite reading: %f\n", difftime(recv_time, satellitelogs[i].timeLog));
                            
                            fprintf(log,"-----------------------------------------\n");
                            printf("-----------------------------------------\n");
                            
                            //prevents double logging of the same alert
                            prev_false_iter = iteration;
                        
                        }
                        
                    }

                }

            }
        }

 
        sleep(1);
        iteration++;
    }
    
    printf("\n[BASE STATION] exit\n");
    printf("\n\n\n--- End of session ---");
    return 0;
}

//worker nodes (WSN)
int slave_io(MPI_Comm world_comm, MPI_Comm comm){
    int ndims=2, size, my_rank, reorder, my_cart_rank, ierr, worldSize;
    int nrows, ncols; //keep track of dimension
    
    //neighbours
    int nbr_i_lo, nbr_i_hi;
    int nbr_j_lo, nbr_j_hi;

    
    int dims[ndims],coord[ndims];
    int wrap_around[ndims];
    
    /* start up initial MPI environment */
    MPI_Comm_size(world_comm, &worldSize); // size of the master communicator
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm comm2D;
    
    
    /* process command line arguments*/
    if (agc == 3) {
        dims[0] = rows; /* number of rows */
        dims[1] = cols; /* number of columns */
        if( (rows*cols) != size) {
            if( my_rank ==0) printf("ERROR: nrows*ncols = %d * %d = %d != %d\n", rows, cols, rows*cols,size);
            return 0;
        }
    } else {
        nrows=ncols=(int)sqrt(size);
        dims[0]=dims[1]=0;
    }
    
    
    /*************************************************************/
    /* create cartesian topology for processes */
    /*************************************************************/
    MPI_Dims_create(size, ndims, dims);
    
    if(my_rank == 0){
        MPI_Send(dims,8, MPI_INT, worldSize-1, 0, world_comm);
    }

    if(my_rank==0) printf("Root Rank: %d. Comm Size: %d: Grid Dimension = [%d x %d] \n",my_rank,size,dims[0],dims[1]);
    fflush(stdout);
    /* create cartesian mapping */
    wrap_around[0] = 0;
    wrap_around[1] = 0; /* periodic shift is .false. */
    reorder = 1;
    //printf("%d\n", my_rank);
    ierr =0;
    ierr = MPI_Cart_create(comm, ndims, dims, wrap_around, reorder, &comm2D);
    if(ierr != 0) printf("ERROR[%d] creating CART\n",ierr);
    /* find my coordinates in the cartesian communicator group */
    MPI_Cart_coords(comm2D, my_rank, ndims, coord); // coordinated is returned into the coord array
    /* use my cartesian coordinates to find my rank in cartesian group*/
    MPI_Cart_rank(comm2D, coord, &my_cart_rank);
    /* get my neighbors; axis is coordinate dimension of shift */
    /* axis=0 ==> shift along the rows: P[my_row-1]: P[me] : P[my_row+1] */
    /* axis=1 ==> shift along the columns P[my_col-1]: P[me] : P[my_col+1] */
    

    MPI_Cart_shift( comm2D, SHIFT_ROW, DISP, &nbr_i_lo, &nbr_i_hi);
    MPI_Cart_shift( comm2D, SHIFT_COL, DISP, &nbr_j_lo, &nbr_j_hi );
    
    //used for Isend and Irecv
    MPI_Request send_request[4];
        MPI_Request receive_request;
        MPI_Status send_status[4];
        MPI_Status receive_status;

    //keeps track of neibghour for each of the WSN nodes
    int nbs[4];
    nbs[0] = nbr_i_lo;
    nbs[1] = nbr_i_hi;
    nbs[2] = nbr_j_lo;
    nbs[3] = nbr_j_hi;
    
    MPI_Status exitStatus;  //used for checking for termination signal tag
    MPI_Status status;
    int randomVal;
    
    
    
    unsigned int seed = time(NULL)+ my_rank;    //seed used for rand function that generates temperature of the node
    
    //structure and all the necessary stuff to send and store potential alets
    struct alert send;
    MPI_Datatype Valuetype;
    MPI_Datatype type[5] = {MPI_INT, MPI_INT,MPI_INT,MPI_INT,MPI_CHAR};
    int blocklen[5] = {1,1,4,4,10};
    MPI_Aint disp[5];

    MPI_Get_address(&send.node, &disp[0]);
    MPI_Get_address(&send.temp, &disp[1]);
    MPI_Get_address(&send.nbs, &disp[2]);
    MPI_Get_address(&send.nbsTemp, &disp[3]);
    MPI_Get_address(&send.coords, &disp[4]);

    disp[4] = disp[4] - disp[0];
    disp[3] = disp[3] - disp[0];
    disp[2] = disp[2] - disp[0];
    disp[1] = disp[1] - disp[0];
    disp[0] = 0;

    MPI_Type_create_struct(5, blocklen, disp, type, &Valuetype);
    MPI_Type_commit(&Valuetype);
    
    //flagof the Iprobe
    int flag = 0;
    
    //tolerance between WSN temp and its neighbours temp
    int tolerance = 10;
    
    //main loop runs until base station sends a signal to stop
    while(1){
        //checking if termination signal was sent by the base station
        MPI_Iprobe(worldSize - 1, TERMINATION_SIGNAL, world_comm,&flag, &exitStatus);
        if(flag){
            if(exitStatus.MPI_TAG == TERMINATION_SIGNAL){ //if it was - terminate the WSN nodes
                break;
            }
        }
        
        
        int i, dest, source;
        
        //temp values of the neighbouring nodes
        int recvValues[4];
        recvValues[0] = -1;
        recvValues[1] = -1;
        recvValues[2] = -1;
        recvValues[3] = -1;
        
        //get the number of a WSN node in the grid and its temp
        send.node = my_rank;
        send.temp = randomVal;
        //get the coords of the WSN
        sprintf( send.coords, "(%d,%d)", coord[0], coord[1]);
        
        //random value (temperature) of the WSN node
        randomVal = (rand_r(&seed) % (100 - 40 + 1)) + 40;
        
        //if the temperature greater than 80 degrees
        if (randomVal > 80){
            //for every neighbour of the node
            for (i = 0; i < 4; i++){
                dest = nbs[i];
                source = nbs[i];
                
                //send the value to the neighbour and recieve their values
                MPI_Isend(&randomVal, 1, MPI_INT, dest, 4, comm2D, &send_request[i]);
                MPI_Iprobe(source, 4, comm2D, &flag, &status);
                if(flag){
                    MPI_Irecv(&recvValues[i], 1, MPI_INT, source, 4, comm2D, &receive_request);
                    MPI_Wait(&receive_request, &receive_status);
                }
            }
            
            //wait to send
            MPI_Waitall(1, send_request, send_status);
            
            //count of the neighbours that have temp value withing threshold to the main WSN node
            int count = 0;
            for (i = 0; i < 4; i++){
                if(abs(randomVal - recvValues[i]) < tolerance){ //if withing threshold - record the neighbour as well as its value and increase count by 1
                    send.nbs[i] = nbs[i];
                    send.nbsTemp[i] = recvValues[i];
                    count+=1;
                    
                }else{ //else give negative value to know it is not accepted
                    send.nbs[i] = -10;
                    send.nbsTemp[i] = -10;
                }
            }
            //if number of neighbours that passed threshold 2 and more - potential alert - report to the base station
            if(count>=2){
                clock_gettime(CLOCK_MONOTONIC, &commStart);     //check the time of communication start
                MPI_Send(&send, 1, Valuetype, worldSize-1, 1, world_comm );
               
            }
            

        }
        sleep(1);
    }
    printf("[WSN %d] exit\n", my_rank);
    return 0;
}


//IRS thread function
void* InfraredFunc(void *pArg)
{
    
    int lowerTemp = 50, upperTemp = 100; // Define upper bound and lower bound of temp
    int node;
    int temp;
    int size;
   // srand(time(NULL));
    unsigned int seed = time(NULL);
    
    MPI_Comm_size(MPI_COMM_WORLD, &size ); //get the size (number) of processes
    size = size - 1; //get the number of WSN nodes
    
    // Loop to populate the array
    for (int i = 0; i < SATTELITE_ARRAY; i++){
        temp = (rand_r(&seed) % (upperTemp - lowerTemp + 1)) + lowerTemp;
        node = (rand_r(&seed) % (size - 0 + 1)) + 0;
        time_t ltime; // calendar time
        
        
        satellitelogs[i].timeLog = time(&ltime);;
        satellitelogs[i].node = node;
        satellitelogs[i].temp = temp;
        
        sleep(1);
    }
    printf("[SATELLITE] exit\n");
    return 0;
}

