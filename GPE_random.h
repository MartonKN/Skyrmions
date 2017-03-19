double GPE_ran(long *idum);
// Produces random number between 0 and 1 using multiplicative method (see Numerical Recipies p. 280). Cycle around 1e8.
// Recommended instead of the built in random number generators of the C language implementations.
// idum: shall be initialized as a negative long variable once in the program and not touched afterwards. (This is the seed of ran1.)