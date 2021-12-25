#ifndef BUILD_GRAPH_CALCULATION
#define BUILD_GRAPH_CALCULATION

#ifdef  USE_OMP
#include <omp.h>
#endif 

#include <map>
#include <vector>
#include<set>
#include <bitset>
#include<list>

#include <iomanip>
#include <algorithm>

#include "build_graph_prj_config.h"
#include "build_graph_read_write.h"
#include "build_graph_structures.h"

#define PI 3.14159265358979323846


int FromDecartToSphere(const Type* decart, Type& fi, Type& theta);
int FromSphericalToDecart(const int number_cur, const std::vector<Type>& all_directions, Vector3& direction);
// ----------------------------------------------------------------------------------------------------------------

size_t SetBasis(const Type * start_point, Vector3 & normal, Eigen::Matrix3d & basis);
size_t Make2dPoint(const Type* start, const Eigen::Matrix3d& local_basis, const Type* point, Vector3& new_point);
int IntersectionWithPlane(const Face& face, const Vector3& start_point, const Vector3& direction, Vector3& result);
int InTriangle(int number_face, const Face& cell_face, const Normals& normals_cell, const Eigen::Vector3d& XX);

#ifdef 0
int FindInAndOutFaces(const Vector3& direction, const int number_cell, const std::vector<Normals>& normals, int* face_state);

int FindInAndOutFaces(const Vector3& direction, const int number_cell, const std::vector<Normals>& normals, uint8_t& face_state);

#endif
int FindInAndOutFaces(const Vector3& direction, const int number_cell, const std::vector<Normals>& normals, std::bitset<4>& face_state);

int InitFacesState(const std::vector<IntId>& all_pairs_id, std::vector<State>& faces_state, const std::set<IntId>& inner_boundary_faces);

#ifdef 0
int FindIdCellInBoundary(const Vector3& direction, const std::set<IntId>& inner_bound, const std::map<IntId, Face>& inner_faces,
	const std::vector<Normals>& normals, const int cur_face);

int FindIdCellInBoundary(const Vector3& direction, const std::set<IntId>& inner_bound, const std::map<IntId, Face>& inner_faces,
	const std::vector<Normals>& normals, const int cur_face, int* id);
bool CheckCell(const IntId id_cell, const Vector3& direction, const std::vector<Normals>& normals,
	std::vector<State>& faces_state, const std::vector<IntId>& all_pairs_id, const std::set<IntId>& inner_boundary_faces,
	const std::map<IntId, Face>& inner_faces);

bool CheckCell3(const IntId id_cell, const Vector3& direction, const std::vector<Normals>& normals,
	std::vector<State>& faces_state, const std::vector<IntId>& all_pairs_id, const std::set<IntId>& inner_boundary_faces,
	const std::map<IntId, Face>& inner_faces);

IntId GetNextCell(const Vector3& direction, const std::set<IntId>& boundary_id, const std::vector<Normals>& normals,
	std::vector<State>& faces_state, const std::vector<IntId>& all_pairs_id, const std::set<IntId>& inner_boundary_faces,
	const std::map<IntId, Face>& inner_faces);

IntId GetNextCell2(const Vector3& direction, std::vector<IntId>& boundary_id, const std::vector<Normals>& normals,
	std::vector<State>& faces_state, const std::vector<IntId>& all_pairs_id, const std::set<IntId>& inner_boundary_faces,
	const std::map<IntId, Face>& inner_faces);

IntId GetNextCellNew(const Vector3& direction, const std::set<IntId>& boundary_id, const std::vector<Normals>& normals,
	std::vector<State>& faces_state, const std::vector<IntId>& all_pairs_id, const std::set<IntId>& inner_boundary_faces,
	const std::map<IntId, Face>& inner_faces, std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::vector<IntId>& cur_el);

int ReBuildSetBondary(const IntId id_cell, const Vector3& direction, std::set<IntId>& boundary_id, const std::vector<Normals>& normals,
	std::vector<State>& faces_state, const std::vector<IntId>& all_pairs_id, std::set<IntId>& set_graph);

int ReBuildSetBondaryNew(const Vector3& direction, std::set<IntId>& boundary_id, const std::vector<Normals>& normals,
	std::vector<State>& faces_state, const std::vector<IntId>& all_pairs_id, std::set<IntId>& set_graph, const std::vector<IntId>& cur_el);
//---------------------------------------------------------------------------------------------------
#endif

int FindCurCellWithHole(const std::set<IntId>& next_step_el, const std::vector<IntId>& count_in_face,  std::vector<IntId>& count_knew_face,
	std::vector<IntId>& cur_el,
	const std::set<IntId>& inner_part, std::set<IntId>& outter_part, const std::map<IntId, FaceCell>& inner_cells,
	const std::vector<IntId>& all_pairs_id, const Vector3& direction, const std::vector<Normals>& normals);
#ifdef 0
int FindNumberOfAllInnerFace(const Vector3& dir, const std::vector<Normals>& normals, std::vector<IntId>& count_in_face);
int FindNumberOfAllKnewInnerDFace(const std::vector<IntId>& bound_vect, const std::vector<State>& faces_state, std::vector<IntId>& count_knew_face);

int ReduceOuterBoundary(std::vector<IntId>& bound_vect, const Vector3& direction, const std::vector<Normals>& normals,
	const std::vector<State>& faces_state);

#endif
int InitFacesState(const std::vector<IntId>& all_pairs_id, std::vector<State>& faces_state, const  std::map<IntId, FaceCell>& inner_boundary_faces);
//=======================================OMP=======================================
#ifdef USE_OMP

int NewStep(const std::vector<IntId>& all_pairs_id, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::list<IntId>& cur_el,
	std::vector<IntId>& next_step_el);
int FindCurCell(const std::vector<IntId>& next_step_el, const std::vector<IntId>& count_in_face, const std::vector<IntId>& count_knew_face,
	std::list<IntId>& cur_el);
int FindNumberOfAllInnerFaceAndKnew(const Vector3& dir, const std::vector<Normals>& normals, const std::vector<State>& faces_state,
	std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::vector<IntId>& next_step_el);
int FindCurCellWithHole(const std::vector<IntId>& next_step_el, const std::vector<IntId>& count_in_face, const std::vector<IntId>& count_knew_face,
	std::list<IntId>& cur_el,
	const std::set<IntId>& inner_boundary_cell, const std::map<IntId, FaceCell>& inner_cells,
	const std::vector<IntId>& all_pairs_id, const Vector3& direction, const std::vector<Normals>& normals);
#else
int FindNumberOfAllInnerFaceAndKnew(const Vector3& dir, const std::vector<Normals>& normals, const std::vector<State>& faces_state,
	std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::set<IntId>& next_step_el);

int FindCurCell(const std::set<IntId>& next_step_el, const std::vector<IntId>& count_in_face, const std::vector<IntId>& count_knew_face, std::vector<IntId>& cur_el);
int NewStep(const std::vector<IntId>& all_pairs_id, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, const std::vector<IntId>& cur_el,
	std::set<IntId>& next_step_el);


#endif //USE_OMP

#endif