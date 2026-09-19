#if defined(USE_MPI)
#include <mpi.h>
#endif
#include <stdlib.h>
#include <iostream>


#if defined(__cplusplus)
extern "C" {
#endif

#if defined(__cplusplus)
}
#endif

// C++ function
extern int SimpleAdvect_cpp();


int main(int argc, char* argv[])
{
	SimpleAdvect_cpp();
	return EXIT_SUCCESS;
}
