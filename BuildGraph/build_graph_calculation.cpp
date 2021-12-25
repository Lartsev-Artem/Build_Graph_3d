#include "build_graph_calculation.h"

//extern std::vector<IntId> bound_vect;

//extern std::set<IntId> inner_part;
//extern std::set<IntId> outter_part;

int FromDecartToSphere(const Type* decart, Type& fi, Type& theta) {
	Type x = decart[0];
	Type y = decart[1];

	theta = atan(sqrt(x * x + y * y) / decart[2]) + PI / 2;

	if (x <= 0)
		fi = atan(y / x) + PI;
	else if (x > 0 && y < 0)
		fi = atan(y / x) + 2 * PI;
	else if (x > 0 && y >= 0)
		fi = atan(y / x);

	return 0;
}
int FromSphericalToDecart(const int number_cur, const std::vector<Type>& all_directions, Vector3& direction) {
	Type theta = all_directions[number_cur];
	Type fi = all_directions[all_directions.size() / 2 + number_cur];

	direction[0] = sin(theta) * cos(fi);
	direction[1] = sin(theta) * sin(fi);
	direction[2] = cos(theta);
	return 0;
}

//-------------------------------------------------------------------------------------
size_t SetBasis(const Type* start_point, Vector3& normal, Eigen::Matrix3d& basis) {
	/*�� ��������� ����� � ������� ������ ��������� ����� ��������� ��������� (vec1, vec2).
	  ������� ����. ������ ���� ������ �����������(������������ normal). ������ ������ �� ���������� ������������*/
	Vector3 vec_1;
	Vector3 vec_2;

	if (abs(normal[1]) < 1e-20) {
		vec_1[0] = 0;
		vec_1[1] = 1;
		vec_1[2] = 0;
	}
	else {
		vec_1[0] = 1;
		vec_1[2] = 0;
		vec_1[1] = -(normal[0] * vec_1[0] + normal[2] * vec_1[2]) / normal[1];  //��-�� ���������� ������������ (N, vec1)==0
	}

	// ���������� ���������� ������ ���������
	if (normal[1] < 0)
		for (int i = 0; i < 3; ++i)
			vec_1[i] *= -1;

	// ������� ��������� ���������. Eigen �������� � �� �����!!!
	Eigen::Vector3d c = normal.cross(vec_1);

	for (size_t i = 0; i < 3; ++i)
		vec_2[i] = -c(i);

	vec_1.normalize();
	vec_2.normalize();

	basis.row(0) = vec_1;
	basis.row(1) = vec_2;
	basis.row(2) = normal;

	return 0;
}
size_t Make2dPoint(const Type* start, const Eigen::Matrix3d& local_basis, const Type* point, Vector3& new_point) {



	for (size_t i = 0; i < 3; i++)
		new_point[i] = 0;

	//������� 3d ����� � 2d (� ��������� ������ {start, local_basis}) 
	for (size_t k = 0; k < 3; k++) {
		new_point[0] += (point[k] - start[k]) * local_basis(0, k);
		new_point[1] += (point[k] - start[k]) * local_basis(1, k);
	}
	return 0;
}
int IntersectionWithPlane(const Face& face, const Vector3& start_point, const Vector3& direction, Vector3& result) {

	//������� ������������


	Type a, b, c, d;  // ��������� ��������� ���������
	Type t;

	a = face.A[1] * (face.B[2] - face.C[2]) + face.B[1] * (face.C[2] - face.A[2]) + face.C[1] * (face.A[2] - face.B[2]);
	b = face.A[0] * (face.C[2] - face.B[2]) + face.B[0] * (face.A[2] - face.C[2]) + face.C[0] * (face.B[2] - face.A[2]);
	c = face.A[0] * (face.B[1] - face.C[1]) + face.B[0] * (face.C[1] - face.A[1]) + face.C[0] * (face.A[1] - face.B[1]);
	d = face.A[0] * (face.C[1] * face.B[2] - face.B[1] * face.C[2]) + face.B[0] * (face.A[1] * face.C[2] -
		face.C[1] * face.A[2]) + face.C[0] * (face.B[1] * face.A[2] - face.A[1] * face.B[2]);

	t = -(a * start_point[0] + b * start_point[1] + c * start_point[2] + d) / (a * direction[0] + b * direction[1] + c * direction[2]);

	for (size_t i = 0; i < 3; ++i)
		result[i] = (direction[i] * t + start_point[i]);  // ����� ����������� ����  (start->direction) � ����������!!! face

	return 0;
}
int InTriangle(int number_face, const Face& cell_face, const Normals& normals_cell, const Eigen::Vector3d& XX) {
	/*face --- �����������, X --- ����� ��� ��������*/

	// ������� ������������
	const Type* AA = cell_face.A.data();
	const Type* BB = cell_face.B.data();
	const Type* CC = cell_face.C.data();

	Vector3 A, B, C, X;  // ����� ����� �� ���������
	{
		Eigen::Matrix3d basis;
		Vector3 n = normals_cell.n[number_face % 4];
		SetBasis(AA, n, basis);
		Make2dPoint(AA, basis, AA, A);
		Make2dPoint(AA, basis, BB, B);
		Make2dPoint(AA, basis, CC, C);
		Make2dPoint(AA, basis, XX.data(), X);
	}

	// �������� �������
	Type r1 = (A[0] - X[0]) * (B[1] - A[1]) - (B[0] - A[0]) * (A[1] - X[1]);
	Type r2 = (B[0] - X[0]) * (C[1] - B[1]) - (C[0] - B[0]) * (B[1] - X[1]);
	Type r3 = (C[0] - X[0]) * (A[1] - C[1]) - (A[0] - C[0]) * (C[1] - X[1]);

	if (r1 < 0 && r2 < 0 && r3 < 0)
		return true;
	else if (r1 > 0 && r2 > 0 && r3 > 0)
		return true;
	else return false;
}

#ifdef 0

int FindInAndOutFaces(const Vector3& direction, const int number_cell, const std::vector<Normals>& normals, int* face_state) {
	//face_state  -0=> ��������� �����,  1=> ��������  face_state.size=4!!!  

	Normals cell = normals[number_cell];
	for (size_t i = 0; i < 4; ++i) {

		const Type eps = 1e-6;
		if (cell.n[i].dot(direction) < -eps)
			face_state[i] = 1;
		else if (cell.n[i].dot(direction) > eps)
			face_state[i] = 0;
		else {
			face_state[i] = 0;  // ���� ����� �����������, �������, ��� ��� �� �������� ������������
			//std::cout << "FindInAndOutFaces: error directions\n";
		}
	}

	return 0;
}

int FindInAndOutFaces(const Vector3& direction, const int number_cell, const std::vector<Normals>& normals, uint8_t& face_state) {
	//face_state  -0=> ��������� �����,  1=> ��������  face_state.size=4!!!  

	Normals cell = normals[number_cell];
	const Type eps = 1e-6;
	for (size_t i = 0; i < 4; ++i) {

		if (cell.n[i].dot(direction) < -eps)
			face_state |= (1 << i);
	}

	return 0;
}

#endif // 0

int FindInAndOutFaces(const Vector3& direction, const int number_cell, const std::vector<Normals>& normals, std::bitset<4>& face_state) {
	//face_state  -0=> ��������� �����,  1=> ��������  face_state.size=4!!!  

	Normals cell = normals[number_cell];
	const Type eps = 1e-6;
	for (size_t i = 0; i < 4; ++i) {

		if (cell.n[i].dot(direction) < -eps)
			face_state.set(i);
		else
			face_state.set(i, false);
	}

	return 0;
}

#ifdef 0
int InitFacesState(const std::vector<IntId>& all_pairs_id, std::vector<State>& faces_state, const std::set<IntId>& inner_boundary_faces) {

	const int n = all_pairs_id.size();
	faces_state.assign(n, 0);

	for (size_t i = 0; i < n; i++) {
		if (all_pairs_id[i] == -1)
			faces_state[i] = 1;
	}

	for (auto& el : inner_boundary_faces)
		faces_state[el] = 0;

	return 0;
}
#endif

int InitFacesState(const std::vector<IntId>& all_pairs_id, std::vector<State>& faces_state, const std::map<IntId,FaceCell>& inner_boundary_faces) {

	const int n = all_pairs_id.size();
	faces_state.assign(n, 0);

	for (size_t i = 0; i < n; i++) {
		if (all_pairs_id[i] == -1)
			faces_state[i] = 1;
	}

	for (auto& el : inner_boundary_faces)
		faces_state[el.second.face_id] = 0;

	return 0;
}
#ifdef 0
int FindIdCellInBoundary(const Vector3& direction, const std::set<IntId>& inner_bound, const std::map<IntId, Face>& inner_faces,
	const std::vector<Normals>& normals, const int cur_face) {

	// ������� ������������.
	Vector3 P1(inner_faces.find(cur_face)->second.A.data());
	Vector3 P2(inner_faces.find(cur_face)->second.B.data());
	Vector3 P3(inner_faces.find(cur_face)->second.C.data());

	Vector3 vertex;
	{
		// ���� ���� ����� (�����)

		Type a = (P1 - P2).norm();
		Type b = (P3 - P2).norm();
		Type c = (P1 - P3).norm();
		Type S = a + b + c;

		vertex[0] = (a * P3[0] + b * P1[0] + c * P2[0]) / S;
		vertex[1] = (a * P3[1] + b * P1[1] + c * P2[1]) / S;
		vertex[2] = (a * P3[2] + b * P1[2] + c * P2[2]) / S;
	}

	// ���� ����������� vectrex->direction c ������� ���������� �������
	Vector3 intersect_point;
	int count = 0;
	int result = -1;
	for (auto& in_id : inner_bound) {
		IntersectionWithPlane(inner_faces.find(in_id)->second, vertex, direction, intersect_point);
		if (InTriangle(in_id, inner_faces.find(in_id)->second, normals[in_id / 4], intersect_point))
			if (in_id != cur_face) {
				count++;
				result = in_id;
				//break;
			}
	}

	if (count > 1)
		std::cout << "Inner intersections: " << count << " ??\n";
	if (!count)
		std::cout << "No intersections: " << "\n";

	return result;
}

int FindIdCellInBoundary(const Vector3& direction, const std::set<IntId>& inner_bound, const std::map<IntId, Face>& inner_faces,
	const std::vector<Normals>& normals, const int cur_face, int* id) {

	// ������� ������������.
	Vector3 P1(inner_faces.find(cur_face)->second.A.data());
	Vector3 P2(inner_faces.find(cur_face)->second.B.data());
	Vector3 P3(inner_faces.find(cur_face)->second.C.data());

	// �������� ������ (�������������� ����� P1,P2,P3 ��������������)
	Vector3 P11 = (P3 + P2) / 2;
	Vector3 P22 = (P3 + P1) / 2;
	Vector3 P33 = (P2 + P1) / 2;

	// ����� �� ��������
	Vector3 vertex1 = P1 + (P11 - P1) / 3;
	Vector3 vertex2 = P2 + (P22 - P2) / 3;
	Vector3 vertex3 = P3 + (P33 - P3) / 3;

	// ���� ����������� vectrex->direction c ������� ���������� �������
	Vector3 intersect_point;
	Face plane_face;
	for (auto& in_id : inner_bound) {
		plane_face = inner_faces.find(in_id)->second;
		IntersectionWithPlane(plane_face, vertex1, direction, intersect_point);
		if (InTriangle(in_id, plane_face, normals[in_id / 4], intersect_point))
			if (in_id != cur_face) {
				id[0] = in_id;
				break;
			}
	}

	for (auto& in_id : inner_bound) {
		plane_face = inner_faces.find(in_id)->second;
		IntersectionWithPlane(plane_face, vertex2, direction, intersect_point);
		if (InTriangle(in_id, plane_face, normals[in_id / 4], intersect_point))
			if (in_id != cur_face) {
				id[1] = in_id;
				break;
			}
	}

	for (auto& in_id : inner_bound) {
		plane_face = inner_faces.find(in_id)->second;
		IntersectionWithPlane(plane_face, vertex3, direction, intersect_point);
		if (InTriangle(in_id, plane_face, normals[in_id / 4], intersect_point))
			if (in_id != cur_face) {
				id[2] = in_id;
				break;
			}
	}

	return 0;
}
bool CheckCell(const IntId id_cell, const Vector3& direction, const std::vector<Normals>& normals,
	std::vector<State>& faces_state, const std::vector<IntId>& all_pairs_id, const std::set<IntId>& inner_boundary_faces,
	const std::map<IntId, Face>& inner_faces) {


	//IntId face_state[4];
	uint8_t state = 0;

	//FindInAndOutFaces(direction, id_cell, normals, face_state);
	FindInAndOutFaces(direction, id_cell, normals, state);

	for (size_t i = 0; i < 4; i++)
		//if (face_state[i]) {// �������� �����
		if (state & (1 << i)) {
			int id_face = id_cell * 4 + i;
			if (faces_state[id_face] == 0 && all_pairs_id[id_face] == -1) {//  �������������� ���������� �������
				if (inner_boundary_faces.count(id_face) != 0) {
					int try_id = FindIdCellInBoundary(direction, inner_boundary_faces, inner_faces, normals, id_face);
					if (try_id == -1) {
						std::cout << "NewCheckCell: wtf??\n";
						continue;
					}

					if (faces_state[try_id] == 1) { // ���� ����� �� ������ ����� ����������
						faces_state[id_face] = 1;
						continue;
					}
					else
						return false;
				}
			}

			if (faces_state[id_cell * 4 + i] == 0) {// && faces_state[all_pairs_id[id_cell * 4 + i]] == 0) {
				// ����� �� ����������
				return false;
			}
		}
	return true;
}

bool CheckCell3(const IntId id_cell, const Vector3& direction, const std::vector<Normals>& normals,
	std::vector<State>& faces_state, const std::vector<IntId>& all_pairs_id, const std::set<IntId>& inner_boundary_faces,
	const std::map<IntId, Face>& inner_faces) {


	//IntId face_state[4];
	uint8_t state = 0;

	//FindInAndOutFaces(direction, id_cell, normals, face_state);
	FindInAndOutFaces(direction, id_cell, normals, state);
	IntId try_id[3] = { -1,-1,-1 };
	for (size_t i = 0; i < 4; i++)
		//if (face_state[i]) {// �������� �����
		if (state & (1 << i)) {
			int id_face = id_cell * 4 + i;
			if (faces_state[id_face] == 0 && all_pairs_id[id_face] == -1) {//  �������������� ���������� �������
				if (inner_boundary_faces.count(id_face) != 0) {
					FindIdCellInBoundary(direction, inner_boundary_faces, inner_faces, normals, id_face, try_id);
					if (try_id[0] == -1 || try_id[1] == -1 || try_id[2] == -1) {
						std::cout << "NewCheckCell: wtf??\n";
						continue;
					}

					if (faces_state[try_id[0]] == 1 && faces_state[try_id[1]] == 1 && faces_state[try_id[2]] == 1) { // ���� ����� �� ������ ����� ����������
						faces_state[id_face] = 1;
						continue;
					}
					else
						return false;
				}
			}

			if (faces_state[id_cell * 4 + i] == 0) {// && faces_state[all_pairs_id[id_cell * 4 + i]] == 0) {
				// ����� �� ����������
				return false;
			}
		}
	return true;
}

IntId GetNextCell(const Vector3& direction, const std::set<IntId>& boundary_id, const std::vector<Normals>& normals,
	std::vector<State>& faces_state, const std::vector<IntId>& all_pairs_id, const std::set<IntId>& inner_boundary_faces,
	const std::map<IntId, Face>& inner_faces) {
	// ���������� (���������)������ ������ �� ����� ������������� ��������� �������

	for (auto el : boundary_id) {
		if (CheckCell3(el, direction, normals, faces_state, all_pairs_id, inner_boundary_faces, inner_faces)) {
			for (size_t i = 0; i < 4; i++) {
				faces_state[el * 4 + i] = 1;
				if (all_pairs_id[el * 4 + i] != -1)
					faces_state[all_pairs_id[el * 4 + i]] = 1;
			}
			return el;
		}//if
	}//for

	std::cout << "NextCell is -1\n";

	return -1;
}

IntId GetNextCell2(const Vector3& direction, std::vector<IntId>& boundary_id, const std::vector<Normals>& normals,
	std::vector<State>& faces_state, const std::vector<IntId>& all_pairs_id, const std::set<IntId>& inner_boundary_faces,
	const std::map<IntId, Face>& inner_faces) {
	// ���������� (���������)������ ������ �� ����� ������������� ��������� �������

	static int num_start_threads = omp_get_num_threads();
	int ret_el = -1;
	int ind = -1;
	bool flag = true;
	//if(boundary_id.size()< num_start_threads+1)omp_set_num_threads(1);
#pragma omp parallel default(none) shared(boundary_id,direction, normals, faces_state, all_pairs_id, inner_boundary_faces, inner_faces,flag, ret_el,ind)  
	{
		int el = -1;

#pragma omp for
		for (int i = 0; i < boundary_id.size(); i++) {
			if (flag) {
				el = boundary_id[i];
				if (CheckCell3(el, direction, normals, faces_state, all_pairs_id, inner_boundary_faces, inner_faces)) {
#pragma omp critical
					{
						if (CheckCell3(el, direction, normals, faces_state, all_pairs_id, inner_boundary_faces, inner_faces))
						{
							for (size_t i = 0; i < 4; i++) {
								faces_state[el * 4 + i] = 1;
								if (all_pairs_id[el * 4 + i] != -1)
									faces_state[all_pairs_id[el * 4 + i]] = 1;
							}
							ret_el = el;
							ind = i;

							flag = false;
							//break; //return el;
						}//if_critical
					}//critical
				}//if
			}//if_flag
		}//for

	}//parallel
	if (flag) {
		std::cout << "NextCell is -1\n";
		return -1;
	}

	std::swap(boundary_id[boundary_id.size() - 1], boundary_id[ind]); // ��������� ������� � ����� ������
	return ret_el;
}

IntId GetNextCellNew2(const Vector3& direction, std::vector<IntId>& boundary_id, const std::vector<Normals>& normals,
	std::vector<State>& faces_state, const std::vector<IntId>& all_pairs_id, const std::set<IntId>& inner_boundary_faces,
	const std::map<IntId, Face>& inner_faces, std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face) {
	// ���������� (���������)������ ������ �� ����� ������������� ��������� �������

	static int num_start_threads = omp_get_num_threads();
	int ret_el = -1;
	int ind = -1;
	bool flag = true;
	//if(boundary_id.size()< num_start_threads+1)omp_set_num_threads(1);
#pragma omp parallel default(none) shared(boundary_id,direction, normals, faces_state, all_pairs_id, inner_boundary_faces, inner_faces,flag, ret_el,ind)  
	{
		int el = -1;

#pragma omp for
		for (int i = 0; i < boundary_id.size(); i++) {
			if (flag) {
				el = boundary_id[i];
				if (CheckCell3(el, direction, normals, faces_state, all_pairs_id, inner_boundary_faces, inner_faces)) {
#pragma omp critical
					{
						if (CheckCell3(el, direction, normals, faces_state, all_pairs_id, inner_boundary_faces, inner_faces))
						{
							for (size_t i = 0; i < 4; i++) {
								faces_state[el * 4 + i] = 1;
								if (all_pairs_id[el * 4 + i] != -1)
									faces_state[all_pairs_id[el * 4 + i]] = 1;
							}
							ret_el = el;
							ind = i;

							flag = false;
							//break; //return el;
						}//if_critical
					}//critical
				}//if
			}//if_flag
		}//for

	}//parallel
	if (flag) {
		std::cout << "NextCell is -1\n";
		return -1;
	}

	std::swap(boundary_id[boundary_id.size() - 1], boundary_id[ind]); // ��������� ������� � ����� ������
	return ret_el;
}

IntId GetNextCellNew(const Vector3& direction, const std::set<IntId>& boundary_id, const std::vector<Normals>& normals,
	std::vector<State>& faces_state, const std::vector<IntId>& all_pairs_id, const std::set<IntId>& inner_boundary_faces,
	const std::map<IntId, Face>& inner_faces, std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::vector<IntId>& cur_el) {
	// ���������� (���������)������ ������ �� ����� ������������� ��������� �������
	cur_el.clear();

	for (auto el : boundary_id) {
		if (CheckCell3(el, direction, normals, faces_state, all_pairs_id, inner_boundary_faces, inner_faces)) {
			for (size_t i = 0; i < 4; i++) {
				faces_state[el * 4 + i] = 1;
				if (all_pairs_id[el * 4 + i] != -1)
					faces_state[all_pairs_id[el * 4 + i]] = 1;
			}
			if (el == 359)
				int a = 0;
			cur_el.push_back(el);  //return el;
		}//if
	}//for

	if (cur_el.size() == 0) {
		std::cout << "NextCell is -1\n";
		return -1;
	}

	return cur_el.size();
}

int ReBuildSetBondary(const IntId id_cell, const Vector3& direction, std::set<IntId>& boundary_id, const std::vector<Normals>& normals,
	std::vector<State>& faces_state, const std::vector<IntId>& all_pairs_id, std::set<IntId>& set_graph) {

	IntId face_state[4];
	FindInAndOutFaces(direction, id_cell, normals, face_state);
	bound_vect.pop_back();
	for (size_t i = 0; i < 4; i++)
		if (face_state[i] == 0) {
			int num_face = all_pairs_id[id_cell * 4 + i] / 4;


			if (set_graph.count(num_face) == 0)
				if (boundary_id.count(num_face) == 0) {
					boundary_id.emplace(num_face);
					bound_vect.push_back(num_face);
				}

			//int buf_s = set_graph.size();
			//set_graph.emplace(num_face / 4);
			//if (set_graph.size() != buf_s) {
			//	boundary_id.emplace(num_face / 4); 
			//	set_graph.erase(num_face / 4);
			//}
		}
	boundary_id.erase(id_cell);

	return 0;
}

int ReBuildSetBondaryNew( const Vector3& direction, std::set<IntId>& boundary_id, const std::vector<Normals>& normals,
	std::vector<State>& faces_state, const std::vector<IntId>& all_pairs_id, std::set<IntId>& set_graph, const std::vector<IntId>& cur_el) {


	for (int hh = 0; hh < cur_el.size(); hh++) {
		int	id_cell = cur_el[hh];

		std::bitset<4> face_state;  //IntId face_state[4];
		FindInAndOutFaces(direction, id_cell, normals, face_state);
		//bound_vect.pop_back();
		for (size_t i = 0; i < 4; i++)
			if (face_state[i] == 0) {
				int num_face = all_pairs_id[id_cell * 4 + i] / 4;


				if (set_graph.count(num_face) == 0)
					if (boundary_id.count(num_face) == 0) {
						boundary_id.emplace(num_face);
						//bound_vect.push_back(num_face);
					}

				//int buf_s = set_graph.size();
				//set_graph.emplace(num_face / 4);
				//if (set_graph.size() != buf_s) {
				//	boundary_id.emplace(num_face / 4); 
				//	set_graph.erase(num_face / 4);
				//}
			}
		boundary_id.erase(id_cell);
	}
	return 0;
}

//---------------------------------------------------------------------------------------------------

// ������������� �������(���������) ������� ��� �������� ����������� (�.�. ������������ "������ ���������")
int ReduceOuterBoundary(std::vector<IntId>& bound_vect, const Vector3& direction, const std::vector<Normals>& normals,
	const std::vector<State>& faces_state) {
	uint8_t state = 0;
	std::list<IntId> new_bound;

	for (int i = 0; i < bound_vect.size(); i++) {

		int cell = bound_vect[i];

		std::bitset<4> state;
		FindInAndOutFaces(direction, cell, normals, state);

		for (int k = 0; k < 4; k++)
			if (state[k] == 1 && faces_state[cell * 4 + k] == 1)  // ���� ����� �������� � ����������
			{
				new_bound.push_back(cell); // ��������� ������
				break;
			}
	}

	bound_vect.assign(new_bound.begin(), new_bound.end());

	return 0;
}

#endif
#ifdef 0
// ����� �������� ������ ��� ������ ������ 
int FindNumberOfAllInnerFace(const Vector3& dir, const std::vector<Normals>& normals, std::vector<IntId>& count_in_face) {

	count_in_face.resize(normals.size(), 0);  // ����� �����

	int state[4] = { -1 };
	for (int i = 0; i < count_in_face.size(); i++)
	{
		// ����� ������� ����� ������������
		FindInAndOutFaces(dir, i, normals, state);
		for (int j = 0; j < 4; j++)
		{
			if (state[j])
				count_in_face[i]++;
		}
	}
	return 0;
}
#endif

#ifdef 0
// ����� ��������� (������������) �������� ������ ��� ������ ������ 
int FindNumberOfAllKnewInnerDFace(const std::vector<IntId>& bound_vect, const std::vector<State>& faces_state, std::vector<IntId>& count_knew_face) {

	count_knew_face.resize(faces_state.size() / 4, 0);  // ����� �����

	for (int i = 0; i < faces_state.size() / 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (faces_state[i * 4 + j])
				count_knew_face[i]++;

		}

	}
	return 0;
}

#endif
int FindIdCellInBoundary(const Vector3& direction, const std::set<IntId>& inner_bound, const std::map<IntId, FaceCell>& inner_cells,
	const std::vector<Normals>& normals, const int cur_cell, int* id) {

	// ������� ������������.
	Face face = inner_cells.find(cur_cell)->second.face;
	Vector3 P1(face.A.data());
	Vector3 P2(face.B.data());
	Vector3 P3(face.C.data());

	// �������� ������ (�������������� ����� P1,P2,P3 ��������������)
	Vector3 P11 = (P3 + P2) / 2;
	Vector3 P22 = (P3 + P1) / 2;
	Vector3 P33 = (P2 + P1) / 2;

	// ����� �� ��������
	Vector3 vertex1 = P1 + (P11 - P1) / 3;
	Vector3 vertex2 = P2 + (P22 - P2) / 3;
	Vector3 vertex3 = P3 + (P33 - P3) / 3;

	// ���� ����������� vectrex->direction c ������� ���������� �������
	Vector3 intersect_point;
	FaceCell plane_face;
	for (auto& in_id : inner_bound) {
		plane_face = inner_cells.find(in_id)->second;
		IntersectionWithPlane(plane_face.face, vertex1, direction, intersect_point);
		if (InTriangle(plane_face.face_id, plane_face.face, normals[in_id], intersect_point))
			//if (in_id != cur_cell && ((intersect_point - vertex1).dot(direction) < 0)) 
		{
			/*	std::bitset<4>id;
				FindInAndOutFaces(direction, in_id, normals, id);
				if (id[plane_face.face_id % 4] == 0) break;*/

			id[0] = in_id;
			break;
		}
	}

	for (auto& in_id : inner_bound) {
		plane_face = inner_cells.find(in_id)->second;
		IntersectionWithPlane(plane_face.face, vertex2, direction, intersect_point);
		if (InTriangle(plane_face.face_id, plane_face.face, normals[in_id], intersect_point))
			//if (in_id != cur_cell && ((intersect_point - vertex2).dot(direction) < 0)) 
		{
			/*	std::bitset<4>id;
				FindInAndOutFaces(direction, in_id, normals, id);
				if (id[plane_face.face_id % 4] == 0) break;*/
			id[1] = in_id;
			break;
		}
	}

	for (auto& in_id : inner_bound) {
		plane_face = inner_cells.find(in_id)->second;
		IntersectionWithPlane(plane_face.face, vertex3, direction, intersect_point);
		if (InTriangle(plane_face.face_id, plane_face.face, normals[in_id], intersect_point))
			//if (in_id != cur_cell && ((intersect_point - vertex3).dot(direction) < 0)) 
		{
			/*	std::bitset<4>id;
				FindInAndOutFaces(direction, in_id, normals, id);
				if (id[plane_face.face_id % 4] == 0) break;*/
			id[2] = in_id;
			break;
		}
	}

	return 0;
}

#ifdef USE_OMP
int FindCurCell(const std::vector<IntId>& next_step_el, const std::vector<IntId>& count_in_face, const std::vector<IntId>& count_knew_face,
	std::list<IntId>& cur_el) {

	cur_el.clear();

	const int N = next_step_el.size();
#pragma omp parallel default(none) shared(next_step_el, count_in_face, count_knew_face,cur_el,  N)  
	{
#pragma omp for
		for (int i = 0; i < N; i++)
		{
			int cell = next_step_el[i];
			if (count_in_face[cell] == count_knew_face[cell]) {
#pragma omp critical 
				{
					if (count_in_face[cell] == count_knew_face[cell])
						cur_el.push_back(cell);
				}
			}
		}
	}

	if (cur_el.size() == 0) {
		std::cout << "NextCell is -1\n";
		return -1;
	}

	return 0;
}
int NewStep(const std::vector<IntId>& all_pairs_id, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::list<IntId>& cur_el,
	std::vector<IntId>& next_step_el) {

	int buf_size = next_step_el.size();
	next_step_el.clear();
	next_step_el.reserve(buf_size); // ����������� ������ �� ������ ����������� ����(�����������, ��� ������� �������� ����� ��� ��)

	int N = cur_el.size();

	std::set<IntId> next_step;

	for (auto cell: cur_el)
//	for (size_t i = 0; i < N; i++)
	{
	//	int cell = cur_el[i];
		// �� ���� �������
		for (size_t j = 0; j < 4; j++)
		{
			int neighbour = all_pairs_id[cell * 4 + j];
			if (neighbour == -1) continue;
			neighbour /= 4;

			if (count_in_face[neighbour] > count_knew_face[neighbour]) {
				count_knew_face[neighbour]++;  // ������ �� ��� ����� ����� �������� (�������� �� ��������??)
				if (next_step.count(neighbour) == 0) {
					next_step.emplace(neighbour);
					next_step_el.push_back(neighbour);  // ������ ���� ��������, ��������� �� ���������� �� ��������� ����
				}
			}
		}

	}

	return 0;
}

// ����� �������� ������ ��� ������ ������ + ����� ��������� �� ��� + ��������� �������
int FindNumberOfAllInnerFaceAndKnew(const Vector3& dir, const std::vector<Normals>& normals, const std::vector<State>& faces_state,
	std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::vector<IntId>& next_step_el) {

	const int N = normals.size();  // ����� �����

	count_in_face.assign(N, 0);
	count_knew_face.assign(N, 0);
	next_step_el.clear();

	std::bitset<4> state;
	for (int i = 0; i < N; i++)
	{
		FindInAndOutFaces(dir, i, normals, state);
		for (int j = 0; j < 4; j++)
		{
			if (state[j]) {// ��� i-�� ������ ���������:
				count_in_face[i]++;  // ������� �����
				if (faces_state[i * 4 + j]) {
					count_knew_face[i]++;  //����������� ������� ����� (�.�. ���������)
					next_step_el.push_back(i);  // ��������� ����� ����������
				}
			}
		}
	}
	return 0;
}


int FindCurCellWithHole(const std::vector<IntId>& next_step_el, const std::vector<IntId>& count_in_face, const std::vector<IntId>& count_knew_face,
	std::list<IntId>& cur_el,
	const std::set<IntId>& inner_boundary_cell, const std::map<IntId, FaceCell>& inner_cells,
	const std::vector<IntId>& all_pairs_id, const Vector3& direction, const std::vector<Normals>& normals) {

	cur_el.clear();

	const int N = next_step_el.size();
	for (auto cell : next_step_el)
	{
		if (inner_boundary_cell.count(cell) != 0) { // ���� ������ �� ���������� �������
			if (count_in_face[cell] > count_knew_face[cell]) {  // ������� �� ����������
				IntId try_id[3] = { -1,-1,-1 };

				FindIdCellInBoundary(direction, inner_boundary_cell, inner_cells, normals, cell, try_id);
				if (try_id[0] == -1 || try_id[1] == -1 || try_id[2] == -1) continue;
				if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
					count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
					count_in_face[try_id[2]] == count_knew_face[try_id[2]]) { // ���� ����� �� ������ ����� ����������
					cur_el.push_back(cell);
					continue;
				}
			}
			else  // ������� ����������
				cur_el.push_back(cell);
		}
		else if (count_in_face[cell] == count_knew_face[cell])
			cur_el.push_back(cell);
	}

	if (cur_el.size() == 0) {

		// �����, �� ��� ����. ���� �� ������ ����� �� ���� ������ ��������� \
		����������� ������ �������� �� ���������� �������, 

		for (auto cell : inner_boundary_cell) {
			if (count_in_face[cell] > count_knew_face[cell]) {  // ������� �� ����������
				IntId try_id[3] = { -1,-1,-1 };
				FindIdCellInBoundary(direction, inner_boundary_cell, inner_cells, normals, cell, try_id);
				if (try_id[0] == -1 || try_id[1] == -1 || try_id[2] == -1) continue;
				if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
					count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
					count_in_face[try_id[2]] == count_knew_face[try_id[2]]) { // ���� ����� �� ������ ����� ����������
					cur_el.push_back(cell);
					continue;
				}
			}
			else  // ������� ����������
				cur_el.push_back(cell);
		}

		if (cur_el.size() == 0) {
			std::cout << "NextCell is -1\n";
			return -1;
		}
	}

	return 0;
}

#else

// ����� �������� ������ ��� ������ ������ + ����� ��������� �� ��� + ��������� �������
int FindNumberOfAllInnerFaceAndKnew(const Vector3& dir, const std::vector<Normals>& normals, const std::vector<State>& faces_state,
	std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::set<IntId>& next_step_el) {

	const int N = normals.size();  // ����� �����

	count_in_face.assign(N, 0);
	count_knew_face.assign(N, 0);
	next_step_el.clear();

	std::bitset<4> state;
	for (int i = 0; i < N; i++)
	{
		FindInAndOutFaces(dir, i, normals, state);
		for (int j = 0; j < 4; j++)
		{
			if (state[j]) {// ��� i-�� ������ ���������:
				count_in_face[i]++;  // ������� �����
				if (faces_state[i * 4 + j]) {
					count_knew_face[i]++;  //����������� ������� ����� (�.�. ���������)
					next_step_el.emplace(i);  // ��������� ����� ����������
				}
			}
		}
	}
	return 0;
}


int FindCurCell(const std::set<IntId>& next_step_el, const std::vector<IntId>& count_in_face, const std::vector<IntId>& count_knew_face, std::vector<IntId>& cur_el) {

	cur_el.clear();

	const int N = next_step_el.size();
	for (auto cell : next_step_el)
	{
		if (count_in_face[cell] == count_knew_face[cell])
			cur_el.push_back(cell);
	}

	if (cur_el.size() == 0) {
		std::cout << "NextCell is -1\n";
		return -1;
	}

	return 0;
}

int NewStep(const std::vector<IntId>& all_pairs_id, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face,
	const std::vector<IntId>& cur_el,
	std::set<IntId>& next_step_el) {

	next_step_el.clear();
	int N = cur_el.size();

	for (size_t i = 0; i < N; i++)
	{
		int cell = cur_el[i];
		//	count_knew_face[cell]++;
			// �� ���� �������
		for (size_t j = 0; j < 4; j++)
		{
			int neighbour = all_pairs_id[cell * 4 + j];
			if (neighbour == -1) continue;
			neighbour /= 4;

			if (count_in_face[neighbour] > count_knew_face[neighbour]) {
				next_step_el.emplace(neighbour);  // ������ ���� ��������, ��������� �� ���������� �� ��������� ����
				count_knew_face[neighbour]++;  // ������ �� ��� ����� ����� �������� (�������� �� ��������??)
			}
		}

	}

	return 0;
}





int FindCurCellWithHole(const std::set<IntId>& next_step_el, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face,
	std::vector<IntId>& cur_el,
	const std::set<IntId>& inner_part, std::set<IntId>& outter_part,
	const std::map<IntId, FaceCell>& inner_cells,
	const std::vector<IntId>& all_pairs_id, const Vector3& direction, const std::vector<Normals>& normals) {

	cur_el.clear();
	std::set<IntId> test;
	const int N = next_step_el.size();
	for (auto cell : next_step_el)
	{
		if (outter_part.count(cell) != 0) {
			//if (inner_boundary_cell.count(cell) != 0) { // ���� ������ �� ���������� �������
			if (count_in_face[cell] == count_knew_face[cell] + 1) {  // ������� �� ����������
				IntId try_id[3] = { -1,-1,-1 };

				FindIdCellInBoundary(direction, inner_part, inner_cells, normals, cell, try_id);
				//FindIdCellInBoundary(direction, inner_boundary_cell, inner_cells, normals, cell, try_id);
				if (try_id[0] == -1 || try_id[1] == -1 || try_id[2] == -1) continue;
				if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
					count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
					count_in_face[try_id[2]] == count_knew_face[try_id[2]]) { // ���� ����� �� ������ ����� ����������
					cur_el.push_back(cell);
					test.emplace(cell);
					outter_part.erase(cell);
					count_knew_face[cell]++;
					continue;
				}
			}
			else if (count_in_face[cell] == count_knew_face[cell]) {  // ������� ����������
				cur_el.push_back(cell);
				//count_knew_face[cell]++;
				test.emplace(cell);
				outter_part.erase(cell);
			}
		}
		else if (count_in_face[cell] == count_knew_face[cell]) {
			cur_el.push_back(cell);
			test.emplace(cell);
		}
	}


	if (cur_el.size() != test.size())
		int a = 0;

	if (cur_el.size() == 0) {

		// �����, �� ��� ����. ���� �� ������ ����� �� ���� ������ ��������� \
		����������� ������ �������� �� ���������� �������, 

		std::list<IntId> buf_erase;

		for (auto cell : outter_part) {
			if (count_in_face[cell] == count_knew_face[cell] + 1) {  // ������� �� ����������
				IntId try_id[3] = { -1,-1,-1 };
				FindIdCellInBoundary(direction, inner_part, inner_cells, normals, cell, try_id);
				if (try_id[0] == -1 || try_id[1] == -1 || try_id[2] == -1) continue;
				if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
					count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
					count_in_face[try_id[2]] == count_knew_face[try_id[2]]) { // ���� ����� �� ������ ����� ����������
					cur_el.push_back(cell);
					count_knew_face[cell]++;
					buf_erase.push_back(cell);//outter_part.erase(cell);
					continue;
				}
			}
			else if (count_in_face[cell] == count_knew_face[cell]) {
				buf_erase.push_back(cell);//outter_part.erase(cell);
				cur_el.push_back(cell);   //��������� ��� ��������� ����� �������!!!
			}
		}

		for (auto el : buf_erase)
		{
			outter_part.erase(el);
		}

		if (cur_el.size() == 0) {
			std::cout << "NextCell is -1\n";
			return -1;
		}
	}

	return 0;
}

#endif // USE_OMP
