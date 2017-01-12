#include "config.h"
#include "rgrid/darrayscatter.h"
#include "kutils/config.h"

#include <fstream>

using namespace rgrid;

using namespace std;

void Config::readConfig(string file) {
#ifdef USE_MPI
	MPI_File fh;
	MPI_Offset filesize;
	if (0 != MPI_File_open(MPI_COMM_WORLD, file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
		throw logic_error("(MPI) Couldn't open \"" + file + "\"");
	if (0 != MPI_File_get_size(fh, &filesize))
		throw logic_error("(MPI) Couldn't get \"" + file + "\" size");
	vector<char> buf(filesize);
	MPI_File_read_all(fh, &buf[0], filesize, MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);

	std::istringstream ss(std::string(&buf[0],static_cast<size_t>(filesize)));
	istream& is = ss;
#else
	fstream fs;
	fs.open(file.c_str(), ios_base::in);
	if (!fs.is_open()) throw logic_error("Couldn't open \"" + file + "\"");
	istream& is = fs;
#endif

	kconf.parseStream(is);

	int_t order = kconf.get<int_t>("order");
	if (!(order > 0 && order % 2 == 0)) throw logic_error("Wrong order");
	ho = order / 2;
	fdc.calc(ho);

	std::vector<int> a;
	kconf.getEntry("nodes").getIntArray(a);
	// ghost nodes -1 and nx are boundary nodes on the grid
	// in input nx and ny is number of nodes without ghost
	nx = a.at(0);
	ny = a.at(1);

	std::vector<double> b;
	b.clear();
	kconf.getEntry("origin").getDoubleArray(b);
	oy = b.at(0);

	b.clear();
	kconf.getEntry("space_step").getDoubleArray(b);
	dx = b.at(0);
	dy = b.at(1);

	steps = kconf.get<int_t>("time_steps");
	dt = kconf.get<real_t>("time_step");

	ex = nx * dx;
	ey = (ny + 1) * dy + oy;

	saveStep = kconf.get<int_t>("save_every");

	max_pml = kconf.get<real_t>("pml_max");
	pml_len = kconf.get<real_t>("pml_nodes");

	std::vector<std::string> s;
	std::string boundary;
	boundary = kconf.get<std::string>("rad_boundary");
	isPmlRad = (boundary == "absorb");
	boundary = kconf.get<std::string>("top_boundary");
	isPmlTop = (boundary == "absorb");
	boundary = kconf.get<std::string>("bot_boundary");
	isPmlBot = (boundary == "absorb");

	a.clear();
	kconf.getEntry("global_parts").getIntArray(a);
	gx = a.at(0);
	gy = a.at(1);

	a.clear();
	kconf.getEntry("local_parts").getIntArray(a);
	lx = a.at(0);
	ly = a.at(1);

	string sourcesFile = kconf.get<std::string>("sources_position");
	string receiversFile = kconf.get<std::string>("receivers_position");
	string KFile = kconf.get<std::string>("bulk_modulus");

	s.clear();
	kconf.getEntry("density").getStringArray(s);
	string rhoXFile = s.at(0);
	string rhoYFile = s.at(1);

	rcvsOut = kconf.get<std::string>("receivers_output");
#ifndef USE_MPI
	fs.close();
#endif
	readSources(sourcesFile);
	readReceivers(receiversFile);
	readK(KFile);
	readRhoX(rhoXFile);
	readRhoY(rhoYFile);
}

void Config::readSources(string file) {
	std::istringstream ss;
#ifdef USE_MPI
	{
		MPI_File fh;
		MPI_Offset filesize;
		if (0 != MPI_File_open(MPI_COMM_WORLD, file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
			throw logic_error("(MPI) Couldn't open \"" + file + "\"");
		if (0 != MPI_File_get_size(fh, &filesize))
			throw logic_error("(MPI) Couldn't get \"" + file + "\" size");
		vector<char> buf(filesize);
		MPI_File_read_all(fh, &buf[0], filesize, MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_close(&fh);
		ss.str(std::string(&buf[0],static_cast<size_t>(filesize)));
	}
#else
	{
		fstream fs;
		fs.open(file.c_str(), ios_base::in);
		if (!fs.is_open()) throw logic_error("Couldn't open \"" + file + "\"");
		std::string contents((std::istreambuf_iterator<char>(fs)), std::istreambuf_iterator<char>());
		fs.close();
		ss.str(contents);
	}
#endif
	while (1) {
		Source e;
		string file_src;
		ss >> file_src >> e.y;
		if (ss.eof()) break;
		if (!(oy + dy <= e.y && e.y <= ey - dy)) throw logic_error("Wrong source position Y");

		// find nodes
		e.j = static_cast<int_t>((e.y - oy) / dy) - 1;

		// distance from node with indexes i,j,k to current node
		real_t yl = e.y - oy - (e.j + 1) * dy;

		real_t yr = dy - yl;

		real_t volume = dy;
		e.c[L] = yr / volume;
		e.c[R] = yl / volume;

		istringstream ss_src;
#ifdef USE_MPI
	{
		string& file = file_src;
		MPI_File fh;
		MPI_Offset filesize;
		if (0 != MPI_File_open(MPI_COMM_WORLD, file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
			throw logic_error("(MPI) Couldn't open \"" + file + "\"");
		if (0 != MPI_File_get_size(fh, &filesize))
			throw logic_error("(MPI) Couldn't get \"" + file + "\" size");
		vector<char> buf(filesize);
		MPI_File_read_all(fh, &buf[0], filesize, MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_close(&fh);
		ss_src.str(std::string(&buf[0],static_cast<size_t>(filesize)));
	}
#else
	{
		string& file = file_src;
		fstream fs;
		fs.open(file.c_str(), ios_base::in);
		if (!fs.is_open()) throw logic_error("Couldn't open \"" + file + "\"");
		std::string contents((std::istreambuf_iterator<char>(fs)), std::istreambuf_iterator<char>());
		fs.close();
		ss_src.str(contents);
	}
#endif
		for (int_t i = 0; i != steps; ++i) {
			real_t val;
			ss_src >> val;
			e.val.push_back(val);
		}
		src.push_back(e);
	}
}

void Config::readReceivers(string file) {
	std::istringstream ss;
#ifdef USE_MPI
	{
		MPI_File fh;
		MPI_Offset filesize;
		if (0 != MPI_File_open(MPI_COMM_WORLD, file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
			throw logic_error("(MPI) Couldn't open \"" + file + "\"");
		if (0 != MPI_File_get_size(fh, &filesize))
			throw logic_error("(MPI) Couldn't get \"" + file + "\" size");
		vector<char> buf(filesize);
		MPI_File_read_all(fh, &buf[0], filesize, MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_close(&fh);
		ss.str(std::string(&buf[0],static_cast<size_t>(filesize)));
	}
#else
	{
		fstream fs;
		fs.open(file.c_str(), ios_base::in);
		if (!fs.is_open()) throw logic_error("Couldn't open \"" + file + "\"");
		std::string contents((std::istreambuf_iterator<char>(fs)), std::istreambuf_iterator<char>());
		fs.close();
		ss.str(contents);
	}
#endif
	while (1) {
		Receiver e;
		ss >> e.x >> e.y;
		if (ss.eof()) break;

		if (!(0 <= e.x && e.x <= ex)) throw logic_error("Wrong receiver position X");
		if (!(oy <= e.y && e.y <= ey)) throw logic_error("Wrong receiver position Y");

		//e.i = static_cast<int_t>((e.x - ox) / dx) - 1;
		e.i = static_cast<int_t>(e.x / dx);
		e.j = static_cast<int_t>((e.y - oy) / dy) - 1;

		// distance from node with indexes i,j,k to current node
		real_t xl = e.x - e.i * dx;
		real_t yl = e.y - oy - (e.j + 1) * dy;

		if (e.i == nx) { xl = dx; } // e.x == ex
		if (e.j == ny) { e.j -= 1; yl = dy; } // e.y == ey

		real_t xr = dx - xl;
		real_t yr = dy - yl;

		real_t volume = dx * dy;
		e.c[L][L] = xr * yr / volume;
		e.c[L][R] = xr * yl / volume;
		e.c[R][L] = xl * yr / volume;
		e.c[R][R] = xl * yl / volume;

		rcv.push_back(e);
	}
}

void Config::readK(string file) {
	K.setSizes(
		Dim3D<int_t>(nx, ny, 1),
		Dim3D<int_t>(gx, gy, 1),
		Dim3D<int_t>(lx, ly, 1),
		Dim3D<int_t>(ho, ho, 0),
		1);
	K.loadData(file);
	K.externalSyncStart();
	K.internalSync();
	K.fillGhost();
	K.externalSyncEnd();
}

void Config::readRhoX(string file) {
	// last block in global partitioning takes one additional node
	Dim3D<vector<int_t> > globalWidth = K.getWidth();
	globalWidth[X].at(globalWidth[X].size()-1) -= 1;
	Dim3D<vector<int_t> > localWidth = K.getLocalContainer().getWidth();
	// so if last block on current process we add node here
	if (K.getInternalPos(X) == gx - 1) {
		localWidth[X].at(localWidth[X].size()-1) -= 1;
	}
	rhox.setSizes(
		globalWidth,
		localWidth,
		Dim3D<int_t>(ho, ho, 0),
		1);
	rhox.loadData(file);
	rhox.externalSyncStart();
	rhox.internalSync();
	rhox.fillGhost();
	rhox.externalSyncEnd();
}

void Config::readRhoY(string file) {
	Dim3D<vector<int_t> > globalWidth = K.getWidth();
	globalWidth[Y].at(globalWidth[Y].size()-1) -= 1;
	Dim3D<vector<int_t> > localWidth = K.getLocalContainer().getWidth();
	if (K.getInternalPos(Y) == gy - 1) {
		localWidth[Y].at(localWidth[Y].size()-1) -= 1;
	}
	rhoy.setSizes(
		globalWidth,
		localWidth,
		Dim3D<int_t>(ho, ho, 0),
		1);
	rhoy.loadData(file);
	rhoy.externalSyncStart();
	rhoy.internalSync();
	rhoy.fillGhost();
	rhoy.externalSyncEnd();
}
