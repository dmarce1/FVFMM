#include "grid.hpp"
#include <silo.h>
#include <fstream>
#include <list>
#include <set>
#include <cmath>

typedef std::array<real, NDIM> xpoint;

inline bool float_eq(real a, real b) {
	return std::abs(a - b) < 1.0e-10;
}

inline bool xpoint_eq(const xpoint& a, const xpoint& b) {
	bool rc = true;
	for (integer d = 0; d != NDIM; ++d) {
		if (!float_eq(a[d], b[d])) {
			rc = false;
			break;
		}
	}
	return rc;
}

struct node_point {
	xpoint pt;
	integer index;
	bool operator==(const node_point& other) const {
		return xpoint_eq(other.pt, pt);
	}
	bool operator<(const node_point& other) const {
		bool rc = false;
		for (integer d = 0; d != NDIM; ++d) {
			if (!float_eq(pt[d], other.pt[d])) {
				rc = (pt[d] < other.pt[d]);
				break;
			}
		}
		return rc;
	}
};

void grid::output(const char* filename) {

	if (!system("cp hello.chk goodbye.chk\n")) {
	}
	save("hello.chk");

	const integer vertex_order[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };
	std::set<node_point> node_list;

	std::list<integer> zone_list;
	const integer this_bw = HBW;
	for (integer i = this_bw; i != HNX - this_bw; ++i) {
		for (integer j = this_bw; j != HNX - this_bw; ++j) {
			for (integer k = this_bw; k != HNX - this_bw; ++k) {
				for (integer ci = 0; ci != NVERTEX; ++ci) {
					const integer iii = i * DNX + j * DNY + k * DNZ;
					const integer vi = vertex_order[ci];
					const integer xi = (vi >> 0) & 1;
					const integer yi = (vi >> 1) & 1;
					const integer zi = (vi >> 2) & 1;
					node_point this_x;
					this_x.pt[XDIM] = X[XDIM][iii] + (real(xi) - HALF) * dx;
					this_x.pt[YDIM] = X[YDIM][iii] + (real(yi) - HALF) * dx;
					this_x.pt[ZDIM] = X[ZDIM][iii] + (real(zi) - HALF) * dx;
					auto iter = node_list.find(this_x);
					integer index;
					if (iter != std::end(node_list)) {
						index = iter->index;
					} else {
						index = node_list.size();
						this_x.index = index;
						node_list.insert(this_x);
					}
					zone_list.push_back(index);
				}
			}
		}
	}

	const int nzones = zone_list.size() / NVERTEX;
	std::vector<int> zone_nodes(nzones * NVERTEX);
	integer index = 0;
	for (auto iter = std::begin(zone_list); iter != std::end(zone_list); ++iter) {
		zone_nodes[index] = *iter;
		++index;
	}

	const int nnodes = node_list.size();
	std::vector<double> x_coord(nnodes);
	std::vector<double> y_coord(nnodes);
	std::vector<double> z_coord(nnodes);
	std::array<double*, NDIM> node_coords = { x_coord.data(), y_coord.data(), z_coord.data() };
	for (auto iter = std::begin(node_list); iter != std::end(node_list); ++iter) {
		x_coord[iter->index] = iter->pt[0];
		y_coord[iter->index] = iter->pt[1];
		z_coord[iter->index] = iter->pt[2];
	}

	constexpr
	int nshapes = 1;
	int shapesize[1] = { NVERTEX };
	int shapetype[1] = { DB_ZONETYPE_HEX };
	int shapecnt[1] = { nzones };
	const char* coord_names[NDIM] = { "x", "y", "z" };

	DBfile *db = DBCreateReal(filename, DB_CLOBBER, DB_LOCAL, "Euler Mesh", DB_PDB);
	DBPutZonelist2(db, "zones", nzones, int(NDIM), zone_nodes.data(), nzones * NVERTEX, 0, 0, 0, shapetype, shapesize,
			shapecnt, nshapes, nullptr);
	DBPutUcdmesh(db, "mesh", int(NDIM), coord_names, node_coords.data(), nnodes, nzones, "zones", nullptr, DB_DOUBLE,
			nullptr);

	const char* field_names[] = { "rho", "egas", "sx", "sy", "sz", "tau", "pot", "phi", "gx", "gy", "gz" };
	const integer I3 = std::pow(HNX - 2 * this_bw, 3);
	std::array<double*, NF + NGF> u_data;
	for (integer field = 0; field != NF + NGF; ++field) {
		u_data[field] = new double[I3];
		index = 0;
		for (integer i = this_bw; i != HNX - this_bw; ++i) {
			for (integer j = this_bw; j != HNX - this_bw; ++j) {
				for (integer k = this_bw; k != HNX - this_bw; ++k) {
					const integer iii = DNX * i + DNY * j + DNZ * k;
					if (field < NF) {
						u_data[field][index] = U[field][iii];
					} else {
						u_data[field][index] = G[field - NF][iii];
					}
					++index;
				}
			}
		}
	}
	for (int field = 0; field != NF + NGF; ++field) {
		DBPutUcdvar1(db, field_names[field], "mesh", u_data[field], nzones, nullptr, 0, DB_DOUBLE, DB_ZONECENT,
				nullptr);
	}
	for (int field = 0; field != NF + NGF; ++field) {
		delete[] u_data[field];
	}
	DBClose(db);
}

void grid::save(const char* filename) const {
	std::ofstream ofs(filename);
	boost::archive::binary_oarchive arc(ofs);
	arc << *this;
	ofs.close();
}

void grid::load(const char* filename) {
	std::ifstream ifs(filename);
	boost::archive::binary_iarchive arc(ifs);
	arc >> *this;
	ifs.close();
}

