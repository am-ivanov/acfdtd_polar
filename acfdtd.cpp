#include "config.h"

#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include "rgrid/vtksaver.h"
#include "rgrid/darrayscatter.h"

#include <fenv.h>

using namespace std;
using namespace rgrid;

// first and second derivatives in pml layer
enum {
	PHI1X,
	PHI1Y,
	PHI2X,
	PHI2Y,
	PHI1R,
	PHINUM
};

struct PMLParams {
	real_t pmlVal;
	real_t a;
	real_t b;
};

int main(int argc, char** argv) {
	feenableexcept(FE_INVALID | FE_OVERFLOW);  // Enable all floating point exceptions but FE_INEXACT
	rgmpi::init(&argc, &argv);
	{
	if (argc != 2) {
		cout << "Usage: " << argv[0] << " config" << endl;
		return 0;
	}
	Config cfg;
	cfg.readConfig(argv[1]);

	for (int i = -cfg.ho; i <= cfg.ho; ++i) {
		if (i != 0)
			cout << cfg.fdc.sc1(i) << " ";
	}
	cout << endl;
	for (int i = -cfg.ho; i <= cfg.ho; ++i) {
		cout << cfg.fdc.c1(i) << " ";
	}
	cout << endl;
	for (int i = -cfg.ho; i <= cfg.ho; ++i) {
		cout << cfg.fdc.c2(i) << " ";
	}
	cout << endl;

	if (rgmpi::worldRank() == 0) {
		cout << "order = " << cfg.ho * 2 << endl;
		cout << "nodes = " << cfg.nx << ", " << cfg.ny << endl;
		cout << "origin = " << cfg.ox << ", " << cfg.oy << endl;
		cout << "space_step = " << cfg.dx << ", " << cfg.dy << endl;
		cout << "time_steps = " << cfg.steps << endl;
		cout << "time_step = " << cfg.dt << endl;
		cout << "save_every = " << cfg.saveStep << endl;
		cout << "pml_max = " << cfg.max_pml << endl;
		cout << "pml_nodes = " << cfg.pml_len << endl;
		cout << "rad_boundary = " << cfg.isPmlRad << endl;
		cout << "top_boundary = " << cfg.isPmlTop << endl;
		cout << "bot_boundary = " << cfg.isPmlBot << endl;
		cout << "global_parts = " << cfg.gx << ", " << cfg.gy << endl;
		cout << "local_parts = " << cfg.lx << ", " << cfg.ly << endl;

		cout << "Sources: " << cfg.src.size() << endl;
		for (vector<Source>::size_type i = 0; i < cfg.src.size(); ++i) {
			cout << "source " << i << " position:" << endl;
			cout
				<< "index: "
				<< cfg.src.at(i).j
				<< ", coord: "
				<< cfg.src.at(i).y
				<< endl;
		}

		cout << "Receivers: " << cfg.src.size() << endl;
		for (vector<Receiver>::size_type i = 0; i < cfg.src.size(); ++i) {
			cout << "receiver " << i << " position:" << endl;
			cout
				<< "index: "
				<< cfg.rcv.at(i).j
				<< ", coord: "
				<< cfg.rcv.at(i).y
				<< endl;
		}
	}

	DArrayScatter<real_t, int_t> das;
	das.setSizes(
		Dim3D<int_t>(cfg.nx, cfg.ny, 1),
		Dim3D<int_t>(cfg.gx, cfg.gy, 1),
		Dim3D<int_t>(cfg.lx, cfg.ly, 1),
		Dim3D<int_t>(cfg.ho, cfg.ho, 0),
		1);
	DArrayScatter<real_t, int_t> dasNext;
	dasNext.setSizes(
		Dim3D<int_t>(cfg.nx, cfg.ny, 1),
		Dim3D<int_t>(cfg.gx, cfg.gy, 1),
		Dim3D<int_t>(cfg.lx, cfg.ly, 1),
		Dim3D<int_t>(cfg.ho, cfg.ho, 0),
		1);
	DArrayScatter<float, int_t> dasSave;
	dasSave.setSizes(
		Dim3D<int_t>(cfg.nx, cfg.ny, 1),
		Dim3D<int_t>(cfg.gx, cfg.gy, 1),
		Dim3D<int_t>(1, 1, 1),
		Dim3D<int_t>(0, 0, 0),
		1);

	VTKSaver<float, int_t> vs;
	vs.setHeaderStructPoints(
		"Acoustic FDTD with CPML",
		Dim3D<int_t>(cfg.nx, cfg.ny, 1),
		Dim3D<float>(0, 0, 0),
		Dim3D<float>(cfg.dx, cfg.dy, 1));
	vs.appendData("u", dasSave, VTKSaver<float, int_t>::POINT_DATA);

	DArrayScatter<real_t, int_t> x1das, y1das;
	x1das.setSizes(
		cfg.rhox.getWidth(),
		cfg.rhox.getLocalContainer().getWidth(),
		Dim3D<int_t>(cfg.ho, 0, 0),
		1);
	y1das.setSizes(
		cfg.rhoy.getWidth(),
		cfg.rhoy.getLocalContainer().getWidth(),
		Dim3D<int_t>(1, cfg.ho, 0),
		1);

	DArrayScatter<real_t, int_t>* u = &das;
	DArrayScatter<real_t, int_t>* un = &dasNext;

	DArray<real_t, int_t> pmlTop;
	DArray<real_t, int_t> pmlBot;
	DArray<real_t, int_t> pmlRad;

	//pmlTop.resize(cfg.nx, cfg.pml_len, 1, cfg.nx, cfg.pml_len, 1, 0, 0, 0, 1, 0, 0, PHINUM);
	//pmlBot.resize(cfg.nx, cfg.pml_len, 1, cfg.nx, cfg.pml_len, 1, 0, 0, 0, 1, 0, 0, PHINUM);
	pmlTop.resize(cfg.nx, cfg.pml_len, 1, cfg.nx, cfg.pml_len, 1, 0, 0, 0, 0, 0, 0, PHINUM);
	pmlBot.resize(cfg.nx, cfg.pml_len, 1, cfg.nx, cfg.pml_len, 1, 0, 0, 0, 0, 0, 0, PHINUM);
	pmlRad.resize(cfg.pml_len, cfg.ny, 1, cfg.pml_len, cfg.ny, 1, 0, 0, 0, 0, 0, 0, PHINUM);

	pmlTop.fill(0);
	pmlBot.fill(0);
	pmlRad.fill(0);

	vector<PMLParams> pmlParams1; // first derivative
	vector<PMLParams> pmlParams2; // second derivative
	for (int_t i = 0; i != cfg.pml_len; ++i) {
		PMLParams p;
		p.pmlVal = cfg.max_pml * sqr(1 - (i + 0.5)/cfg.pml_len);
		p.b = exp(-p.pmlVal*cfg.dt);
		p.a = p.b - 1;
		pmlParams1.push_back(p);
	}
	for (int_t i = 0; i != cfg.pml_len; ++i) {
		PMLParams p;
		p.pmlVal = cfg.max_pml * sqr(1 - (i * 1.0)/cfg.pml_len);
		p.b = exp(-p.pmlVal * cfg.dt);
		p.a = p.b - 1;
		pmlParams2.push_back(p);
	}

#ifdef USE_MPI
	MPI_File fh;
	if (0 != MPI_File_open(MPI_COMM_WORLD, cfg.rcvsOut.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh))
		throw logic_error("Can't open receivers file for writing");
	if (0 != MPI_File_set_size(fh, 0))
		throw logic_error("Can't truncate receivers file");
#else
	fstream recvsout;
	if (rgmpi::worldRank() == 0) {
		recvsout.open(cfg.rcvsOut.c_str(), ios_base::out | ios_base::trunc);
		if (!recvsout.is_open()) throw logic_error("Can't open receivers file for writing");
	}
#endif

	real_t t = 0.0;
	for (int_t step = 0; step != cfg.steps; ++step) {

		DArrayContainer<real_t, int_t>& dac = u->getLocalContainer();
		DArrayContainer<real_t, int_t>& dacNext = un->getLocalContainer();

		if (step % cfg.saveStep == 0) {
			// save to VTK
			if (rgmpi::worldRank() == 0)
				cout << "Step " << step << endl;
			char name[100];
			sprintf(name, "out-%06ld.vtk", (long)step);
			DArray<float, int_t>& ds = dasSave.getDArrayPart(0, 0, 0);
			for (int_t k = 0; k != dac.numParts(Z); ++k)
			for (int_t j = 0; j != dac.numParts(Y); ++j)
			for (int_t i = 0; i != dac.numParts(X); ++i) {
				DArray<real_t, int_t>& d = dac.getDArrayPart(i, j, k);
				for (int_t k2 = 0; k2 != d.localSize(Z); ++k2)
				for (int_t j2 = 0; j2 != d.localSize(Y); ++j2)
				for (int_t i2 = 0; i2 != d.localSize(X); ++i2) {
					ds(i2 + dac.partOrigin(X, i), j2 + dac.partOrigin(Y, j), k2 + dac.partOrigin(Z, k), 0) = d(i2, j2, k2, 0);
				}
			}
			ds.inverseBytes();
			vs.save(std::string(name));
		}

		u->externalSyncStart();
		u->internalSync();
		u->externalSyncEnd();

		// save receivers
		for (vector<Receiver>::size_type i = 0; i != cfg.rcv.size(); ++i) {
			const Dim3D<int_t> rind(cfg.rcv.at(i).i, cfg.rcv.at(i).j, 0);
			if (u->isPresentGhostGlobal(rind)) {
				float val[LR][LR] = {{0}};
				val[L][L] = u->getNodeGhostLocal(Dim3D<int_t>(rind[X],rind[Y],rind[Z]),0) * cfg.rcv.at(i).c[L][L];
				val[L][R] = u->getNodeGhostLocal(Dim3D<int_t>(rind[X],rind[Y]+1,rind[Z]),0) * cfg.rcv.at(i).c[L][R];
				val[R][L] = u->getNodeGhostLocal(Dim3D<int_t>(rind[X]+1,rind[Y],rind[Z]),0) * cfg.rcv.at(i).c[R][L];
				val[R][R] = u->getNodeGhostLocal(Dim3D<int_t>(rind[X]+1,rind[Y]+1,rind[Z]),0) * cfg.rcv.at(i).c[R][R];
				float rval =
					val[L][L] +
					val[R][L] +
					val[L][R] +
					val[R][R];
#ifdef USE_MPI
				if (0 != MPI_File_write_at(fh, (step * cfg.rcv.size() + i) * sizeof(float), &rval, 1, MPI_FLOAT, MPI_STATUS_IGNORE))
					throw logic_error("Can't write in receivers file");
#else
				recvsout.write(reinterpret_cast<const char*>(&rval), sizeof(float));
#endif
			}
		}

		// recalculate own DArrays
		for (int_t gk = 0; gk != dac.numParts(Z); ++gk)
			for (int_t gj = 0; gj != dac.numParts(Y); ++gj)
				for (int_t gi = 0; gi != dac.numParts(X); ++gi) {
					DArray<real_t, int_t>& p = dac.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& x1 = x1das.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& y1 = y1das.getDArrayPart(gi, gj, gk);

					DArray<real_t, int_t>& rhox = cfg.rhox.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& rhoy = cfg.rhoy.getDArrayPart(gi, gj, gk);

					// inner area
					for (int_t j = 0; j != rhox.localSize(Y); ++j) {
						for (int_t i = 0; i != rhox.localSize(X); ++i) {
							int_t i2 = i + rhox.origin(X);
							int_t j2 = j + rhox.origin(Y);
							x1(i, j, 0, 0) = 0;
							for (int_t n = 0; n < cfg.ho; ++n)
								x1(i, j, 0, 0) += cfg.fdc.sc1(n+1) * (p(i+1+n, j, 0, 0) - p(i-n, j, 0, 0));
							x1(i,j,0,0) /= cfg.dx;
							if (i2 < cfg.ho) {
								// copy axis nodes
								x1(-i-1,j,0,0) = -x1(i,j,0,0);
							}
							if (i2 > cfg.nx-cfg.pml_len-2 && cfg.isPmlRad) {
								pmlRad(cfg.nx-i2-2,j2,0,PHI1X) = pmlParams1.at(cfg.nx-i2-2).b * pmlRad(cfg.nx-i2-2,j2,0,PHI1X) + pmlParams1.at(cfg.nx-i2-2).a * x1(i,j,0,0);
								x1(i,j,0,0) = x1(i,j,0,0) + pmlRad(cfg.nx-i2-2,j2,0,PHI1X);
							}
						}
					}

					for (int_t j = 0; j != rhoy.localSize(Y); ++j) {
						//for (int_t i = -1; i != rhoy.localSize(X); ++i) {
						for (int_t i = 0; i != rhoy.localSize(X); ++i) {
							int_t i2 = i + rhoy.origin(X);
							int_t j2 = j + rhoy.origin(Y);
							y1(i,j,0,0) = 0;
							for (int_t n = 0; n < cfg.ho; ++n)
								y1(i,j,0,0) += cfg.fdc.sc1(n+1) * (p(i,j+1+n,0,0) - p(i,j-n,0,0));
							y1(i,j,0,0) /= cfg.dy;
							if (j2 < cfg.pml_len && cfg.isPmlTop) {
								pmlTop(i2,j2,0,PHI1Y) = pmlParams1.at(j2).b * pmlTop(i2,j2,0,PHI1Y) + pmlParams1.at(j2).a * y1(i,j,0,0);
								y1(i,j,0,0) = y1(i,j,0,0) + pmlTop(i2,j2,0,PHI1Y);
							} else if (j2 > cfg.ny-cfg.pml_len-2 && cfg.isPmlBot) {
								pmlBot(i2,cfg.ny-j2-2,0,PHI1Y) = pmlParams1.at(cfg.ny-j2-2).b * pmlBot(i2,cfg.ny-j2-2,0,PHI1Y) + pmlParams1.at(cfg.ny-j2-2).a * y1(i,j,0,0);
								y1(i,j,0,0) = y1(i,j,0,0) + pmlBot(i2,cfg.ny-j2-2,0,PHI1Y);
							}
							y1(i,j,0,0) /= rhoy(i,j,0,0);
						}
					}

				}

		x1das.externalSyncStart();
		y1das.externalSyncStart();

		x1das.internalSync();
		y1das.internalSync();

		x1das.externalSyncEnd();
		y1das.externalSyncEnd();

		for (int_t gk = 0; gk != dac.numParts(Z); ++gk)
			for (int_t gj = 0; gj != dac.numParts(Y); ++gj)
				for (int_t gi = 0; gi != dac.numParts(X); ++gi) {

					DArray<real_t, int_t>& p = dac.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& pn = dacNext.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& K = cfg.K.getDArrayPart(gi, gj, gk);

					DArray<real_t, int_t>& x1 = x1das.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& y1 = y1das.getDArrayPart(gi, gj, gk);

					DArray<real_t, int_t>& rhox = cfg.rhox.getDArrayPart(gi, gj, gk);

					for (int_t j = 0; j != p.localSize(Y); ++j) {
						for (int_t i = 0; i != p.localSize(X); ++i) {
							real_t x2 = 0, y2 = 0;
							real_t x1r = 0;
							int_t i2 = i + p.origin(X);
							int_t j2 = j + p.origin(Y);

							// assume rhox constant
							for (int_t n = 0; n < cfg.ho; ++n)
								x2 += cfg.fdc.sc1(n+1) * (x1(i+n,j,0,0) - x1(i-1-n,j,0,0));
							x2 /= cfg.dx;
							if (i2 > cfg.nx - cfg.pml_len - 1 && cfg.isPmlRad) {
								pmlRad(cfg.nx-i2-1,j2,0,PHI2X) = pmlParams2.at(cfg.nx-i2-1).b * pmlRad(cfg.nx-i2-1,j2,0,PHI2X) + pmlParams2.at(cfg.nx-i2-1).a * x2;
								x2 = x2 + pmlRad(cfg.nx-i2-1,j2,0,PHI2X);
							}
							x2 /= rhox(i,j,0,0);

							for (int_t n = 0; n < cfg.ho; ++n)
								y2 += cfg.fdc.sc1(n+1) * (y1(i,j+n,0,0) - y1(i,j-1-n,0,0));
							y2 /= cfg.dy;
							if (j2 < cfg.pml_len && cfg.isPmlTop) {
								pmlTop(i2,j2,0,PHI2Y) = pmlParams2.at(j2).b * pmlTop(i2,j2,0,PHI2Y) + pmlParams2.at(j).a * y2;
								y2 = y2 + pmlTop(i2,j2,0,PHI2Y);
							} else if (j2 > cfg.ny - cfg.pml_len - 1 && cfg.isPmlBot) {
								pmlBot(i2,cfg.ny-j2-1,0,PHI2Y) = pmlParams2.at(cfg.ny-j2-1).b * pmlBot(i2,cfg.ny-j2-1,0,PHI2Y) + pmlParams2.at(cfg.ny-j2-1).a * y2;
								y2 = y2 + pmlBot(i2,cfg.ny-j2-1,0,PHI2Y);
							}

							if (i2 == 0) x1r = x2;
							else {
								// assume rhox constant
								for (int_t n = 1; n <= cfg.ho; ++n)
									x1r += cfg.fdc.c1(n) * (p(i+n,j,0,0) - p(i-n,j,0,0));
								x1r /= cfg.dx;
								if (i2 > cfg.nx - cfg.pml_len - 1 && cfg.isPmlRad) {
									pmlRad(cfg.nx-i2-1,j2,0,PHI1R) = pmlParams2.at(cfg.nx-i2-1).b * pmlRad(cfg.nx-i2-1,j2,0,PHI1R) + pmlParams2.at(cfg.nx-i2-1).a * x1r;
									x1r = x1r + pmlRad(cfg.nx-i2-1,j2,0,PHI1R);
								}
								//x1r /= (i2 + 1) * cfg.dx; // divide by radius
								x1r /= i2 * cfg.dx; // divide by radius
								x1r /= rhox(i,j,0,0); // divide by density
							}

							pn(i, j, 0, 0) =
								2.0 * p.val(i, j, 0, 0) - pn(i, j, 0, 0)
								+ K(i, j, 0, 0) * sqr(cfg.dt) * (x2 + y2 + x1r);

							// fill axis ghost nodes
							// calculate them in 3d cartesian coordinates
							if (i2 < cfg.ho) {
							//if (i2 < cfg.ho-1) {
								//pn(-i-2, j, 0, 0) = pn(i, j, 0, 0);
								pn(-i-1, j, 0, 0) = pn(i+1, j, 0, 0);
							}
							//if (i2 == 0) {
							//	real_t x2 = cfg.fdc.c2(0) * p(-1, j, 0, 0);
							//	for (int_t n = 1; n <= cfg.ho; ++n)
							//		x2 += cfg.fdc.c2(n) * 2 * p(-1+n,j,0,0);
							//	x2 /= sqr(cfg.dx);
							//	x2 /= rhox(0, j, 0, 0);

							//	real_t y2 = 0;
							//	for (int_t n = 0; n < cfg.ho; ++n)
							//		y2 += cfg.fdc.sc1(n+1) * (y1(-1,j+n,0,0) - y1(-1,j-1-n,0,0));
							//	y2 /= cfg.dy;

							//	if (j2 < cfg.pml_len && cfg.isPmlTop) {
							//		pmlTop(-1,j2,0,PHI2Y) = pmlParams2.at(j2).b * pmlTop(-1,j2,0,PHI2Y) + pmlParams2.at(j).a * y2;
							//		y2 = y2 + pmlTop(-1,j2,0,PHI2Y);
							//	} else if (j2 > cfg.ny - cfg.pml_len - 1 && cfg.isPmlBot) {
							//		pmlBot(-1,cfg.ny-j2-1,0,PHI2Y) = pmlParams2.at(cfg.ny-j2-1).b * pmlBot(-1,cfg.ny-j2-1,0,PHI2Y) + pmlParams2.at(cfg.ny-j2-1).a * y2;
							//		y2 = y2 + pmlBot(-1,cfg.ny-j2-1,0,PHI2Y);
							//	}

							//	pn(-1, j, 0, 0) =
							//		2.0 * p.val(-1, j, 0, 0) - pn(-1, j, 0, 0)
							//		+ K(-1, j, 0, 0) * sqr(cfg.dt) * (x2 + x2 + y2);
							//}
						}
					}

				}

		// insert source
		for (vector<Source>::iterator it = cfg.src.begin(); it != cfg.src.end(); ++it) {
			real_t val = sqr(cfg.dt) * it->val.at(step) / (sqr(cfg.dx) * cfg.dy * M_PI / 4.0);
			if (un->isPresentGhostGlobal(Dim3D<int_t>(0,it->j,0))) {
				//un->getNodeGhostGlobal(Dim3D<int_t>(-1,it->j,0),0)
				un->getNodeGhostGlobal(Dim3D<int_t>(0,it->j,0),0)
					+= it->c[L] * val
					* cfg.K.getNode(Dim3D<int_t>(0,it->j,0),0);
			}
			if (un->isPresentGhostGlobal(Dim3D<int_t>(0,it->j+1,0))) {
				//un->getNodeGhostGlobal(Dim3D<int_t>(-1,it->j+1,0),0)
				un->getNodeGhostGlobal(Dim3D<int_t>(0,it->j+1,0),0)
					+= it->c[R] * val
					* cfg.K.getNode(Dim3D<int_t>(0,it->j+1,0),0);
			}
		}

		swap(un, u);

		t += cfg.dt;
	}

#ifdef USE_MPI
	MPI_File_close(&fh);
#else
	if (rgmpi::worldRank() == 0) {
		recvsout.close();
	}
#endif
	}
	rgmpi::forceFinalize();
}
