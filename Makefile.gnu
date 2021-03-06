# Makefile variables
# GCC

CC = mpicc
CFLAGS = -O2 -g

CXX = mpicxx
CXXFLAGS = -O2 -g

# Entries to bring executable up to date

all: mpi_integration prompted_mpi_integration print_sin \
     thread_multiple_no_thread thread_multiple_with_thread ioda

mpi_integration: mpi_integration.o
	$(CC) $(CFLAGS) -o mpi_integration mpi_integration.o -lm

mpi_integration.o: mpi_integration.c
	$(CC) $(CFLAGS) -c mpi_integration.c

prompted_mpi_integration: prompted_mpi_integration.o
	$(CC) $(CFLAGS) -o prompted_mpi_integration prompted_mpi_integration.o -lm

prompted_mpi_integration.o: prompted_mpi_integration.c
	$(CC) $(CFLAGS) -c prompted_mpi_integration.c

print_sin: print_sin.c
	$(CC) $(CFLAGS) -o print_sin print_sin.c -lm

thread_multiple_with_thread: thread_multiple_with_thread.cpp
	$(CXX) $(CXXFLAGS) -openmp -o thread_multiple_with_thread thread_multiple_with_thread.cpp

thread_multiple_no_thread: thread_multiple_no_thread.cpp
	$(CXX) $(CXXFLAGS) -openmp -o thread_multiple_no_thread thread_multiple_no_thread.cpp

ioda: ioda.c
	$(CC) $(CFLAGS) -o ioda ioda.c

clean:
	rm -f mpi_integration mpi_integration.o prompted_mpi_integration \
          prompted_mpi_integration.o thread_multiple_with_thread \
          thread_multiple_no_thread print_sin slurm-*.out ioda
