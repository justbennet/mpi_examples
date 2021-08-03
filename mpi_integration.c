#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#ifndef PI
  #define PI 3.14159265358979323846
#endif

int timeval_subtract(struct timeval *result, struct timeval *t2,
                     struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec)
                     - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}

int gettimeofday(struct timeval *tv, struct timezone *tz);


int main(int argc, char **argv)
{
  /*****************************************************************/
  /* Problem description: Numerically integrate the normal distribution
     N(mu,t) from a to b.
     \int_a^b 1/sqrt(2*\pi t^2) * exp(-(x-mu)^2/2t^2) dx = \int_a^b f(x) dx  
     by the midpoint rule:
     ~ \sum_{i=0}^{n-1} f(x_i) \Delta x,
     x_i = a + (i+1/2)*dx
     For the purposes of this demonstration, we assume that only
     the process with rank 0 has access to the parameters of the 
     problem: a, b, t, n, and the specific form of f(x). 
  */
  /*****************************************************************/


  /*****************************************************************/
  /* MPI setup */
  /* Assign to each process a rank, in the range [0,1,...,numprocs-1]
     numprocs is the total number of processors used, X, selected by
     calling this program as:
     
     mpirun -np X mpicode
  */
  MPI_Init(&argc, &argv);
  int rank, numprocs, i;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  // Some timing variables
  struct timeval tvStart, tvFinish, tvSub, tvBegin, tvEnd, tvDiff;
  /*****************************************************************/


  /*****************************************************************/
  gettimeofday(&tvStart, NULL);
  /*****************************************************************/


  /*****************************************************************/
  /* Allow processor 0 to read in the parameters. Allocate memory
     on all processors to store the parameters. Recall, any
     commands *not* encased in a "if(rank == XXX)" statement will
     be executed by *all* processes.
  */
  double *pd = malloc(4 * sizeof(double));
  int n;
  if(rank == 0)
  {
    //  integration lower bound
    pd[0] = -3.0;
    //  integration upper bound
    pd[1] =  3.0;
    //  distribution mean
    pd[2] =  0.0;
    //  distribution std deviation
    pd[3] =  2.0;
    //  number of intervals for numerical integration
    n = 1000000000;
  }
    /* Broadcast this information to all processes. 
       Arguments: (starting memory address, number of values, datatype,
       sending rank, communicator (generally, MPI_COMM_WORLD))
     */
    for (i = 0; i < 50; i++)
    /* for (i = 0; i < 150; i++) */
  {
    gettimeofday(&tvBegin, NULL);
    /* gettimeofday(&tvBegin, NULL); */
    MPI_Bcast(pd, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    double a, b, mu, t;
    a  = pd[0];
    b  = pd[1];
    mu = pd[2];
    t  = pd[3];
    /*****************************************************************/
    /* Determine which part of the interval belongs to each processor: partition n */
    int ii, nmin, nmax;
    double sum;
    double xi,zi,dx;
    /* n/numprocs is the number of intervals per processor
       take the floor gives just the integer part, multiplying
       by rank gives 0*ints, 1*ints, 2*ints, 3*ints, the starting
       point, then rank+1*ints gives the upper limit  */
    nmin = floor(rank * ((double)n / (double)numprocs));
    nmax = floor((rank+1) * ((double)n / (double)numprocs));
  
    /* Compute the terms of the sum belonging to this processor */
    dx = (b - a) / (double)n;
    sum = 0.0; /* sum is a variable local to each processor */
    for(ii = nmin; ii<nmax; ii++)
    {
      xi = (double)(a + (0.5 + (double)ii) * dx);
      zi = (xi - mu) / t;
      sum += exp(-(zi * zi) / 2.0);
    }
    sum *= dx / sqrt(2.0 * PI * t * t);
    /*****************************************************************/
    /* Collect the terms of sum on processor 0 (totalsum is not valid 
       on other processes) and output
     */
    double totalsum;
    /* Sum up information in "sum" on each processor, store in "totalsum" at rank 0
       Arguments: (input memory address, output memory address, number of arguments, datatype, operation, root processor, communicator)
    */
    MPI_Reduce(&sum, &totalsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0) {
      gettimeofday(&tvEnd, NULL);
      printf("The sum = %.6lf\n", totalsum);
      timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
      printf("Elapsed time is:  %ld.%06ld\n", tvDiff.tv_sec, tvDiff.tv_usec);
    }
  }
  free(pd);
  if (rank == 0)
  {
    gettimeofday(&tvFinish, NULL);
    timeval_subtract(&tvSub, &tvFinish, &tvStart);
    printf("Total time is:  %ld.%06ld\n", tvSub.tv_sec, tvSub.tv_usec);
  }
  /*****************************************************************/
  /* End MPI */
  MPI_Finalize();
  return(0);
}
