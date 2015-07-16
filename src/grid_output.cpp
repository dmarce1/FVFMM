#include "grid.hpp"
#include <silo.h>
#include <fstream>
#include <cmath>
#include <thread>

bool grid::float_eq(real a, real b) {
	return std::abs(a - b) < 1.0e-10;
}

bool grid::xpoint_eq(const xpoint& a, const xpoint& b) {
	bool rc = true;
	for (integer d = 0; d != NDIM; ++d) {
		if (!float_eq(a[d], b[d])) {
			rc = false;
			break;
		}
	}
	return rc;
}

bool grid::node_point::operator==(const node_point& other) const {
	return xpoint_eq(other.pt, pt);
}

bool grid::node_point::operator<(const node_point& other) const {
	bool rc = false;
	for (integer d = 0; d != NDIM; ++d) {
		if (!float_eq(pt[d], other.pt[d])) {
			rc = (pt[d] < other.pt[d]);
			break;
		}
	}
	return rc;
}

void grid::merge_output_lists(grid::output_list_type& l1, const grid::output_list_type& l2) {
	std::map<integer, integer> index_map;

	for (auto i = l2.nodes.begin(); i != l2.nodes.end(); ++i) {
		integer index, oindex;
		auto this_x = *i;
		oindex = this_x.index;
		auto j = l1.nodes.find(this_x);
		if (j != l1.nodes.end()) {
			index = j->index;
		} else {
			index = l1.nodes.size();
			this_x.index = index;
			l1.nodes.insert(this_x);
		}
		index_map[oindex] = index;
	}
	for (auto i = l2.zones.begin(); i != l2.zones.end(); ++i) {
		l1.zones.push_back(index_map[*i]);
	}
	for (integer field = 0; field < NF + NGF; ++field) {
		std::vector<real> data(l1.data[field].size() + l2.data[field].size());
		std::copy(l1.data[field].begin(), l1.data[field].end(), data.begin());
		std::copy(l2.data[field].begin(), l2.data[field].end(), data.begin() + l1.data[field].size());
		l1.data[field] = std::move(data);
	}
}

grid::output_list_type grid::get_output_list() const {
	output_list_type rc;
	const integer vertex_order[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };

	std::set<node_point>& node_list = rc.nodes;
	std::list<integer>& zone_list = rc.zones;
	std::array<std::vector<real>, NF + NGF>& data = rc.data;

	for (integer field = 0; field != NF + NGF; ++field) {
		data[field].resize(INX * INX * INX);
	}
	const integer this_bw = HBW;
	integer di = 0;
	for (integer i = this_bw; i != HNX - this_bw; ++i) {
		for (integer j = this_bw; j != HNX - this_bw; ++j) {
			for (integer k = this_bw; k != HNX - this_bw; ++k) {
				const integer iii = i * DNX + j * DNY + k * DNZ;
				for (integer ci = 0; ci != NVERTEX; ++ci) {
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
					if (iter != node_list.end()) {
						index = iter->index;
					} else {
						index = node_list.size();
						this_x.index = index;
						node_list.insert(this_x);
					}
					zone_list.push_back(index);
				}
				for (integer field = 0; field != NF; ++field) {
					data[field][di] = U[field][iii];

				}
				for (integer field = 0; field != NGF; ++field) {
					data[field + NF][di] = G[field][iii];
				}
				++di;
			}
		}
	}

	return rc;
}

void grid::output(const output_list_type& olists, const char* filename) {

	std::thread(
			[=]() {
				const std::set<node_point>& node_list = olists.nodes;
				const std::list<integer>& zone_list = olists.zones;

				const int nzones = zone_list.size() / NVERTEX;
				std::vector<int> zone_nodes(nzones * NVERTEX);
				integer index = 0;
				for (auto iter = zone_list.begin(); iter != zone_list.end(); ++iter) {
					zone_nodes[index] = *iter;
					++index;
				}

				const int nnodes = node_list.size();
				std::vector<double> x_coord(nnodes);
				std::vector<double> y_coord(nnodes);
				std::vector<double> z_coord(nnodes);
				std::array<double*, NDIM> node_coords = {x_coord.data(), y_coord.data(), z_coord.data()};
				for (auto iter = node_list.begin(); iter != node_list.end(); ++iter) {
					x_coord[iter->index] = iter->pt[0];
					y_coord[iter->index] = iter->pt[1];
					z_coord[iter->index] = iter->pt[2];
				}

				constexpr
				int nshapes = 1;
				int shapesize[1] = {NVERTEX};
				int shapetype[1] = {DB_ZONETYPE_HEX};
				int shapecnt[1] = {nzones};
				const char* coord_names[NDIM] = {"x", "y", "z"};

				DBfile *db = DBCreateReal(filename, DB_CLOBBER, DB_LOCAL, "Euler Mesh", DB_PDB);
				DBPutZonelist2(db, "zones", nzones, int(NDIM), zone_nodes.data(), nzones * NVERTEX, 0, 0, 0, shapetype, shapesize,
						shapecnt, nshapes, nullptr);
				DBPutUcdmesh(db, "mesh", int(NDIM), coord_names, node_coords.data(), nnodes, nzones, "zones", nullptr, DB_DOUBLE,
						nullptr);

				const char* field_names[] = {"rho", "egas", "sx", "sy", "sz", "tau", "pot", "phi", "gx", "gy", "gz"};
				for (int field = 0; field != NF + NGF; ++field) {
					DBPutUcdvar1(db, field_names[field], "mesh", olists.data[field].data(), nzones, nullptr, 0, DB_DOUBLE, DB_ZONECENT,
							nullptr);
				}
				DBClose(db);
			}).join();
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

