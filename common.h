#ifndef COMMON_H
#define COMMON_H

#include <vector>
#include <string>

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_FLOAT
typedef float real_t;
#else

typedef double real_t;
#endif

typedef long int_t;

enum Side { L = 0, R = 1, LR = 2 };

enum Boundary {
	FREE_BOUNDARY,
	OPEN_BOUNDARY
};

struct InterpPoint {
	real_t x, y;
	real_t c[LR][LR]; // coefficient for bilinear interpolation
	int_t i, j; // coordinate of the left down bot node
};

struct Source {
	real_t y;
	real_t c[LR];
	real_t j;
	std::vector<real_t> val;
};

struct Receiver : public InterpPoint {
};

template <typename T>
inline T sqr(T x) {
	return x * x;
}

#endif // COMMON_H
