#ifndef BUILD_GRAPH_MAIN
#define BUILD_GRAPH_MAIN

#include <map>
#include <vector>
#include <bitset>

#include <iomanip>
#include <algorithm>

#include <omp.h>

#define PI 3.14159265358979323846

#include "build_graph_structures.h"
#include "build_graph_prj_config.h"
#include "build_graph_read_write.h"
#include "build_graph_calculation.h"

#ifdef USE_MPI
#include "mpi.h"
#define RETURN(a) return a;
#define MPI_END MPI_Finalize();
#define MPI_RETURN(a) { MPI_END RETURN(a) }
#else
#define MPI_RETURN(a) return a;
#endif //USE_MPI


template<typename Str>
int ReadStartSettings(Str name_file_settings, Str& name_file_vtk, Str& name_file_sphere_direction, Str& name_file_graph,
	Str& name_file_normals, Str& name_file_pairs, Str& name_file_inner_boundary, Str& name_file_init_boundary, Str& name_file_inner_sphere) {

	std::ifstream ifile;
	ifile.open(name_file_settings);
	if (!ifile.is_open()) {
		std::cerr << " Error : file settings build graph is not open !\n";
		return 1;
	}

	std::getline(ifile, name_file_vtk);
	std::getline(ifile, name_file_sphere_direction);
	std::getline(ifile, name_file_graph);
	std::getline(ifile, name_file_normals);
	std::getline(ifile, name_file_pairs);
	std::getline(ifile, name_file_inner_boundary);
	std::getline(ifile, name_file_init_boundary);
	std::getline(ifile, name_file_inner_sphere);

	ifile.close();
	return 0;
}

#endif