#ifndef BUILD_GRAPH_READ_WRITE
#define BUILD_GRAPH_READ_WRITE

#include "build_graph_prj_config.h"
#include "build_graph_structures.h"

#include <math.h>

#include<fstream>
#include<iostream>

#include<set>
#include<string>
#include<vector>

#include <chrono>
#include <cstdio>
#include <inttypes.h>
#include <memory>
#include< stdio.h >

typedef std::string Str_Type;
using namespace std::chrono;


#ifdef _MSC_VER
#define fwrite_unlocked _fwrite_nolock
#define fread_unlocked  _fread_nolock
#endif


#ifdef OnlyWriteFiles
#include<vtk-9.0/vtkUnstructuredGrid.h>
#include <vtk-9.0\vtkGenericDataObjectReader.h>
#include <vtk-9.0\vtkIdList.h>
#include <vtk-9.0\vtkCellData.h>
#include <vtk-9.0\vtkGenericDataObjectWriter.h>
#include<vtk-9.0\vtkPoints.h>

int BuildSetForClaster(const std::string name_file_vtk, const std::string name_file_pairs,
	const std::string name_file_boundary, const std::string name_file_normals, const std::string name_file_boundary_inner,
	const std::string name_file_face_and_id, const std::string name_file_squares, const std::string name_file_volume, const std::string name_file_centers);

int NormalAndSquareFace(size_t NumberCell, size_t NumberFace, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Type& S, Eigen::Vector3d& n);
int WriteNormalAndSquaresFile(const std::string name_file_normals, const std::string name_file_squares, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid);
int WriteVolume(const std::string name_file_volume, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid);

int WritePairsId(const std::string name_file_pairs, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid);

int WriteInitBoundarySetAndInnerBoundaryFace(const std::string name_file_boundary, const std::string name_file_boundary_inner,
	const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid);

int WriteInnerCellOfSphere(const std::string name_file_inner_sphere, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid,
	const std::string name_file_boundary_inner);

int WriteInnerCellAndIdFaces(const std::string name_file_boundary_inner, const std::string name_file_face_and_id,
	const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid);

#endif

#ifdef OnlyReadFiles
#include "build_graph_calculation.h"
#include <map>

int ReadPairsId(const std::string name_file_pairs, std::vector<IntId>& all_pairs_face);

int ReadInitBoundarySet(const std::string name_file_boundary, std::set<IntId>& boundary_cells);

int ReadInnerBoundary(const std::string name_file_boundary_inner, std::set<IntId>& id_inner_boundary_face);

int ReadInnerCellBoundary(const std::string name_file_boundary_inner, std::set<IntId>& id_inner_boundary_cell);

int ReadNormalFile(const std::string name_file_normals, std::vector<Normals>& normals);

int ReadInnerCellOfSphere(const std::string name_file_inner_sphere, std::vector<Face>& inner_faces);

int ReadInnerCellOfSphereAndId(const std::string name_file_face_and_id, std::map<IntId, Face>& inner_faces);

int ReadInnerCellOfSphereAndId(const std::string name_file_face_and_id, std::map<IntId, FaceCell>& inner_cells);

int ReadSphereDirectionDecartToSpherical(const std::string name_file_sphere_direction, std::vector<Type>& directions_all);

int ReadSphereDirectionDecart(const std::string name_file_sphere_direction, std::vector<Vector3>& directions_all);

int WriteFileGraph(const int i, const std::string& name_file_graph, const std::vector<IntId>& graph);

int WriteFileGraph(std::unique_ptr<FILE, int(*)(FILE*)>& file_graph, std::unique_ptr<FILE, int(*)(FILE*)>& file_id,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_dist, std::unique_ptr<FILE, int(*)(FILE*)>& file_x,
	const int i, const int n, const std::vector<IntId>& graph);
#ifdef USE_VTK
int WriteFileBoundary(const std::string name_file_out, const std::string name_file_graph, const std::string name_file_grid);
#endif

#endif

#ifdef USE_VTK
#include<vtk-9.0/vtkUnstructuredGrid.h>
#include <vtk-9.0\vtkCellData.h>
#include <vtk-9.0\vtkGenericDataObjectReader.h>
#include <vtk-9.0\vtkGenericDataObjectWriter.h>
template<typename Str>
int ReadFileVtk(const Str name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {
	vtkSmartPointer<vtkGenericDataObjectReader> reader_vtk =
		vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader_vtk->ReadAllScalarsOn();
	reader_vtk->SetFileName(name_file_vtk);
	reader_vtk->Update();

	if (reader_vtk->IsFileUnstructuredGrid()) {
		unstructured_grid = reader_vtk->GetUnstructuredGridOutput();
		unstructured_grid->Modified();
	}
	else {
		std::cout << "Error read file\n";
		std::cout << "file_vtk is not UnstructuredGrid\n";
		return 1;
	}

	std::cout << "Grid has Cell: " << unstructured_grid->GetNumberOfCells() << '\n';
	return 0;
}
#endif //USE_VTK

#endif