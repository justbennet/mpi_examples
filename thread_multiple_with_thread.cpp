#include "math.h"
#include "mpi.h"
#include "omp.h"

#include <iostream>
#include "stdlib.h"

int main(int argc,char* args[])
{
    int required = MPI_THREAD_MULTIPLE;
      int provided = 0;
        MPI_Init_thread(NULL,NULL,required,&provided);
          if(provided!=required)
                {
                      std::cout << "Error: MPI thread support insufficient! required " << required << " provided " << provided;
                        abort();
                              
                }
    long n=1000000;
    double start = MPI_Wtime();
    double *d = new double[n];
    double *d2 = new double[n];

#pragma omp parallel for
    for(long i=0;i<n;i++)
    {
        
        d2[i] = cos(d[i])*pow(d[i],3.0);


    }
    delete[] d;
    delete[] d2;


    double end1 = MPI_Wtime();
    std::cout << "MPI_THREAD_MULTIPLE w/o rank: " << end1-start << std::endl;

    d = new double[n];
    d2 = new double[n];

#pragma omp parallel for
    for(long i=0;i<n;i++)
    {
        int myProcID;
        for(int j=0;j<10;j++)
              MPI_Comm_rank(MPI_COMM_WORLD,&myProcID);
        d2[i] = cos(d[i])*pow(d[i],3.0);


    }
    
    double end2 = MPI_Wtime();
    std::cout << "MPI_THREAD_MULTIPLE w rank: " << end2-end1 << std::endl;
    MPI_Finalize();

    return 0;

}

