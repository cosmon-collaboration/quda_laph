
void QL_abort(int status)
{
#ifdef ARCH_SERIAL
 exit(status);
#elif ARCH_PARALLEL
 MPI_abort(
#else
#error "Invalid architecture"
#endif
} 



bool isPrimaryRank()


broadcast
