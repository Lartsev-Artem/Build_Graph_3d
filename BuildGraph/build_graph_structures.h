#ifndef BUILD_GRAPH_STRUCT
#define BUILD_GRAPH_STRUCT

#include <eigen3/Eigen/Dense>
#include <vector>

typedef Eigen::Vector3d Vector3;

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

#endif
