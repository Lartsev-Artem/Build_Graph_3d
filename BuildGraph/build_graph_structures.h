#ifndef BUILD_GRAPH_STRUCT
#define BUILD_GRAPH_STRUCT

#include <eigen3/Eigen/Dense>
#include <vector>

typedef Eigen::Vector3d Vector3;

extern std::vector<int> id_try_surface;		 // id граней, определяющих внутренюю границу
extern std::vector<double> dist_try_surface; // расстояния между точками (через полость внутри) 
extern std::vector<Vector3> x_try_surface;   // x точка выхода

struct Face {
	Vector3 A;
	Vector3 B;
	Vector3 C;
	Face& operator=(const Face& face) {
		A = face.A;
		B = face.B;
		C = face.C;
	}
};
struct Normals {
	std::vector<Vector3> n;
	Normals() {

	}
	Normals(const int size) {
		n.resize(size);
	}
};
struct FaceCell {
	int face_id;
	Face face;
	FaceCell(const int id = 0, const Face& face_init = Face()) {
		face_id = id;
		face = face_init;
	}
};

struct TrySolve {
	int id_1;
	int id_2;
	int id_3;

	double s_1;
	double s_2;
	double s_3;

	Vector3 x1;
	Vector3 x2;
	Vector3 x3;
};

extern TrySolve buf_try;

#endif
