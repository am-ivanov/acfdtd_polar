#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>
#include <map>

#include "kutils/config.h"
#include "fdcoeff.h"
#include "common.h"
#include "rgrid/darrayscatter.h"

class Config {
public:
	void readConfig(std::string file);
	void readRhoX(std::string file);
	void readRhoY(std::string file);
	void readK(std::string file);
	void readSources(std::string file);
	void readReceivers(std::string file);

	int_t nx, ny; // number of nodes
	real_t ox, oy; // origin, start coordinate
	real_t dx, dy; // step size
	real_t ex, ey; // end of computational area

	int_t steps; // number of time steps
	real_t dt; // time step

	int_t lx, ly; // local partitioning
	int_t gx, gy; // global partitioning
	int_t saveStep; // how often save wavefield (if 0 then no saving)

	bool isPmlRad, isPmlTop, isPmlBot; // presense of pml layer on each side in each direction
	int_t pml_len; // number of nodes which are used for pml layer
	real_t max_pml; // pml grows as x^2 from inner area to outer area, so at the boundary it's value is max_pml

	std::string rcvsOut; // where to store receivers

// Ready to work data
	std::vector<Source> src; // coordinates of sources with values at each time step
	std::vector<Receiver> rcv; // coordiantes of receivers

	rgrid::DArrayScatter<real_t, int_t> K; // bulk modulus
	rgrid::DArrayScatter<real_t, int_t> rhox, rhoy; // density

	int_t ho; // half order
	FDCoeff fdc;

	kutils::Config kconf;
};

#endif // CONFIG_H
