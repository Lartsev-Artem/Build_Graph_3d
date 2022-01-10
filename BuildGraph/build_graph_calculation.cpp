#include "build_graph_calculation.h"

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

int IntersectionWithPlane(const Face& face, const Vector3& start_point, const Vector3& direction, Vector3& result) {

	//вершины треугольника


	Type a, b, c, d;  // параметры уравнения плоскости
	Type t;

	a = face.A[1] * (face.B[2] - face.C[2]) + face.B[1] * (face.C[2] - face.A[2]) + face.C[1] * (face.A[2] - face.B[2]);
	b = face.A[0] * (face.C[2] - face.B[2]) + face.B[0] * (face.A[2] - face.C[2]) + face.C[0] * (face.B[2] - face.A[2]);
	c = face.A[0] * (face.B[1] - face.C[1]) + face.B[0] * (face.C[1] - face.A[1]) + face.C[0] * (face.A[1] - face.B[1]);
	d = face.A[0] * (face.C[1] * face.B[2] - face.B[1] * face.C[2]) + face.B[0] * (face.A[1] * face.C[2] -
		face.C[1] * face.A[2]) + face.C[0] * (face.B[1] * face.A[2] - face.A[1] * face.B[2]);

	t = -(a * start_point[0] + b * start_point[1] + c * start_point[2] + d) / (a * direction[0] + b * direction[1] + c * direction[2]);

	for (size_t i = 0; i < 3; ++i)
		result[i] = (direction[i] * t + start_point[i]);  // точка пересечения луча  (start->direction) с плоскостью!!! face

	return 0;
}
int InTriangle(int number_face, const Face& cell_face, const Normals& normals_cell, const Eigen::Vector3d& XX) {
	/*face --- треугольник, X --- точка для проверки*/

	// вершины треугольника
	const Type* AA = cell_face.A.data();
	const Type* BB = cell_face.B.data();
	const Type* CC = cell_face.C.data();

	Vector3 A, B, C, X;  // новые точки на плоскости
	{
		Eigen::Matrix3d basis;
		Vector3 n = normals_cell.n[number_face % 4];
		SetBasis(AA, n, basis);
		Make2dPoint(AA, basis, AA, A);
		Make2dPoint(AA, basis, BB, B);
		Make2dPoint(AA, basis, CC, C);
		Make2dPoint(AA, basis, XX.data(), X);
	}

	// линейная алгебра
	Type r1 = (A[0] - X[0]) * (B[1] - A[1]) - (B[0] - A[0]) * (A[1] - X[1]);
	Type r2 = (B[0] - X[0]) * (C[1] - B[1]) - (C[0] - B[0]) * (B[1] - X[1]);
	Type r3 = (C[0] - X[0]) * (A[1] - C[1]) - (A[0] - C[0]) * (C[1] - X[1]);

	if (r1 < 0 && r2 < 0 && r3 < 0)
		return true;
	else if (r1 > 0 && r2 > 0 && r3 > 0)
		return true;
	else return false;
}

size_t SetBasis(const Type* start_point, Vector3& normal, Eigen::Matrix3d& basis) {
	/*по начальной точке и нормале строит локальный базис картинной плоскости (vec1, vec2).
	  нормаль дана. задаем один вектор произвольно(ортагонально normal). третий вектор из векторного произведения*/
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
		vec_1[1] = -(normal[0] * vec_1[0] + normal[2] * vec_1[2]) / normal[1];  //св-во скалярного произведения (N, vec1)==0
	}

	// правельная ориентация базиса плоскости
	if (normal[1] < 0)
		for (int i = 0; i < 3; ++i)
			vec_1[i] *= -1;

	// обычное векторное умножение. Eigen временно и не нужен!!!
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

	//перевод 3d точки в 2d (в локальном базисе {start, local_basis}) 
	for (size_t k = 0; k < 3; k++) {
		new_point[0] += (point[k] - start[k]) * local_basis(0, k);
		new_point[1] += (point[k] - start[k]) * local_basis(1, k);
	}
	return 0;
}

int FindIdCellInBoundary(const Vector3& direction, const std::set<IntId>& inner_bound, const std::map<IntId, FaceCell>& inner_cells,
	const std::vector<Normals>& normals, const int cur_cell, int* id) {

	// вершины треугольника.
	Face face = inner_cells.find(cur_cell)->second.face;
	Vector3 P1(face.A.data());
	Vector3 P2(face.B.data());
	Vector3 P3(face.C.data());

	// середины сторон (противолежащий точке P1,P2,P3 соответственно)
	Vector3 P11 = (P3 + P2) / 2;
	Vector3 P22 = (P3 + P1) / 2;
	Vector3 P33 = (P2 + P1) / 2;

	// точки на медианах
	Vector3 vertex1 = P1 + (P11 - P1) / 3;
	Vector3 vertex2 = P2 + (P22 - P2) / 3;
	Vector3 vertex3 = P3 + (P33 - P3) / 3;

	// ищем пересечения vectrex->direction c гранями внутренней границы
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

int FindInAndOutFaces(const Vector3& direction, const int number_cell, const std::vector<Normals>& normals, std::bitset<4>& face_state) {
	//face_state  -0=> выходящая грань,  1=> входящая  face_state.size=4!!!  

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

int FractionInnerBoundary(const Vector3& direction, const std::vector<Normals>& normals, const std::map<IntId, FaceCell>& inner_cells,
	const std::set<IntId>& full_boundary, std::set<IntId>& inner_part, std::set<IntId>& outter_part) {

	inner_part.clear();
	outter_part.clear();

	std::bitset<4> state;

	int bound_face_id = -1;
	for (auto el : full_boundary)
	{
		FindInAndOutFaces(direction, el, normals, state);

		bound_face_id = inner_cells.find(el)->second.face_id;

		if (state[bound_face_id % 4] != 0)  // выходящая грань
			outter_part.emplace(el);
		else
			inner_part.emplace(el);
	}

	return 0;
}

//=======================================OMP=======================================
#ifdef USE_OMP

#ifdef GRID_WITH_INNER_BOUNDARY	
int FindCurCellWithHole(const std::vector<IntId>& next_step_el, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face,
	std::list<IntId>& cur_el,
	const std::set<IntId>& inner_part, std::set<IntId>& outter_part, const std::map<IntId, FaceCell>& inner_cells,
	const std::vector<IntId>& all_pairs_id, const Vector3& direction, const std::vector<Normals>& normals) {

	cur_el.clear();

	const int N = next_step_el.size();
#pragma omp parallel default(none) shared(next_step_el, count_in_face, count_knew_face,cur_el,inner_part,outter_part, all_pairs_id, direction,normals,inner_cells, N)  
	{
		int cell;
#pragma omp for
		for (int i = 0; i < N; ++i)
		{
			cell = next_step_el[i];
			if (outter_part.count(cell) != 0) {
				if (count_in_face[cell] == count_knew_face[cell] + 1) {  // граница не определена
					IntId try_id[3] = { -1,-1,-1 };

					FindIdCellInBoundary(direction, inner_part, inner_cells, normals, cell, try_id);

					if (try_id[0] == -1 || try_id[1] == -1 || try_id[2] == -1)
						int a;
					else if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
						count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
						count_in_face[try_id[2]] == count_knew_face[try_id[2]]) { // если грань на другом конце определена
#pragma omp critical 
							{
								if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
									count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
									count_in_face[try_id[2]] == count_knew_face[try_id[2]]) {
									cur_el.push_back(cell);
									outter_part.erase(cell);
									count_knew_face[cell]++;
								}
							}
					}
				}
				else if (count_in_face[cell] == count_knew_face[cell]) {  // граница определена
#pragma omp critical 
					{
						if (count_in_face[cell] == count_knew_face[cell]) {
							cur_el.push_back(cell);
							outter_part.erase(cell);
						}
					}
				}
			}
			else if (count_in_face[cell] == count_knew_face[cell]) {
#pragma omp critical 
				{
					if (count_in_face[cell] == count_knew_face[cell])
						cur_el.push_back(cell);
				}
			}
		}


	}

	if (cur_el.size() == 0) {

		// плохо, но как есть. Если не смогли найти ни одну ячейку кандидата \
		попробовать пройти отдельно по внутренней границе, 

		std::list<IntId> buf_erase;

		for (auto cell : outter_part) {
			if (count_in_face[cell] == count_knew_face[cell] + 1) {  // граница не определена
				IntId try_id[3] = { -1,-1,-1 };
				FindIdCellInBoundary(direction, inner_part, inner_cells, normals, cell, try_id);
				if (try_id[0] == -1 || try_id[1] == -1 || try_id[2] == -1) continue;
				if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
					count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
					count_in_face[try_id[2]] == count_knew_face[try_id[2]]) { // если грань на другом конце определена

					cur_el.push_back(cell);
					count_knew_face[cell]++;
					buf_erase.push_back(cell);//outter_part.erase(cell);
					continue;
				}
			}
			else if (count_in_face[cell] == count_knew_face[cell]) {
				buf_erase.push_back(cell);//outter_part.erase(cell);
				cur_el.push_back(cell);   //добавляет все найденные ранее границы!!!
			}
		}

		for (auto el : buf_erase)
			outter_part.erase(el);

		if (cur_el.size() == 0) {
			std::cout << "NextCell is -1\n";
			return -1;
		}
	}

	return 0;
}
#else
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
#endif

// число входящих граней для каждой ячейки + число известных из них + начальная граница
int FindNumberOfAllInnerFaceAndKnew(const Vector3& dir, const std::vector<Normals>& normals, const std::vector<State>& faces_state,
	std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::vector<IntId>& next_step_el) {

	const int N = normals.size();  // число ячеек

	count_in_face.assign(N, 0);
	count_knew_face.assign(N, 0);
	next_step_el.clear();

	std::bitset<4> state;
	for (int i = 0; i < N; i++)
	{
		FindInAndOutFaces(dir, i, normals, state);
		for (int j = 0; j < 4; j++)
		{
			if (state[j]) {// для i-ой ячейки добавляем:
				count_in_face[i]++;  // входную грань
				if (faces_state[i * 4 + j]) {
					count_knew_face[i]++;  //определнную входную грань (т.е. граничную)
					next_step_el.push_back(i);  // начальный набор кандидатов
				}
			}
		}
	}
	return 0;
}

int NewStep(const std::vector<IntId>& all_pairs_id, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::list<IntId>& cur_el,
	std::vector<IntId>& next_step_el) {

	int buf_size = next_step_el.size();
	next_step_el.clear();
	next_step_el.reserve(buf_size); // резервируем память на основе предыдущего шага(предпологая, что порядок величины будет тот же)

	int N = cur_el.size();

	std::set<IntId> next_step;

	for (auto cell : cur_el)
	{
		// по всем соседям
		for (size_t j = 0; j < 4; j++)
		{
			int neighbour = all_pairs_id[cell * 4 + j];
			if (neighbour == -1) continue;
			neighbour /= 4;

			if (count_in_face[neighbour] > count_knew_face[neighbour]) {
				count_knew_face[neighbour]++;  // всегда ли эта грань будет входящей (проверка по нормалям??)
				if (next_step.count(neighbour) == 0) {
					next_step.emplace(neighbour);
					next_step_el.push_back(neighbour);  // ячейка была изменена, проверить ее готовность на следующем шаге
				}
			}
		}

	}

	return 0;
}

#else

// число входящих граней для каждой ячейки + число известных из них + начальная граница
int FindNumberOfAllInnerFaceAndKnew(const Vector3& dir, const std::vector<Normals>& normals, const std::vector<State>& faces_state,
	std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face, std::set<IntId>& next_step_el) {

	const int N = normals.size();  // число ячеек

	count_in_face.assign(N, 0);
	count_knew_face.assign(N, 0);
	next_step_el.clear();

	std::bitset<4> state;
	for (int i = 0; i < N; i++)
	{
		FindInAndOutFaces(dir, i, normals, state);
		for (int j = 0; j < 4; j++)
		{
			if (state[j]) {// для i-ой ячейки добавляем:
				count_in_face[i]++;  // входную грань
				if (faces_state[i * 4 + j]) {
					count_knew_face[i]++;  //определнную входную грань (т.е. граничную)
					next_step_el.emplace(i);  // начальный набор кандидатов
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

int FindCurCellWithHole(const std::set<IntId>& next_step_el, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face,
	std::vector<IntId>& cur_el, const std::set<IntId>& inner_part, std::set<IntId>& outter_part,
	const std::map<IntId, FaceCell>& inner_cells, const std::vector<IntId>& all_pairs_id, const Vector3& direction, const std::vector<Normals>& normals) {

	cur_el.clear();

	const int N = next_step_el.size();
	for (auto cell : next_step_el)
	{
		if (outter_part.count(cell) != 0) {
			if (count_in_face[cell] == count_knew_face[cell] + 1) {  // граница не определена
				IntId try_id[3] = { -1,-1,-1 };

				FindIdCellInBoundary(direction, inner_part, inner_cells, normals, cell, try_id);

				if (try_id[0] == -1 || try_id[1] == -1 || try_id[2] == -1) continue;

				if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
					count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
					count_in_face[try_id[2]] == count_knew_face[try_id[2]]) { // если грань на другом конце определена

					cur_el.push_back(cell);
					outter_part.erase(cell);
					count_knew_face[cell]++;
					continue;
				}
			}
			else if (count_in_face[cell] == count_knew_face[cell]) {  // граница определена
				cur_el.push_back(cell);
				outter_part.erase(cell);
			}
		}
		else if (count_in_face[cell] == count_knew_face[cell]) {
			cur_el.push_back(cell);
		}
	}

	if (cur_el.size() == 0) {

		// плохо, но как есть. Если не смогли найти ни одну ячейку кандидата \
		попробовать пройти отдельно по внутренней границе, 

		std::list<IntId> buf_erase;

		for (auto cell : outter_part) {
			if (count_in_face[cell] == count_knew_face[cell] + 1) {  // граница не определена
				IntId try_id[3] = { -1,-1,-1 };
				FindIdCellInBoundary(direction, inner_part, inner_cells, normals, cell, try_id);
				if (try_id[0] == -1 || try_id[1] == -1 || try_id[2] == -1) continue;
				if (count_in_face[try_id[0]] == count_knew_face[try_id[0]] &&
					count_in_face[try_id[1]] == count_knew_face[try_id[1]] &&
					count_in_face[try_id[2]] == count_knew_face[try_id[2]]) { // если грань на другом конце определена

					cur_el.push_back(cell);
					count_knew_face[cell]++;
					buf_erase.push_back(cell);//outter_part.erase(cell);
					continue;
				}
			}
			else if (count_in_face[cell] == count_knew_face[cell]) {
				buf_erase.push_back(cell); //outter_part.erase(cell);
				cur_el.push_back(cell);
			}
		}

		for (auto el : buf_erase)
			outter_part.erase(el);

		if (cur_el.size() == 0) {
			std::cout << "NextCell is -1\n";
			return -1;
		}
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
			// по всем соседям
		for (size_t j = 0; j < 4; j++)
		{
			int neighbour = all_pairs_id[cell * 4 + j];
			if (neighbour == -1) continue;
			neighbour /= 4;

			if (count_in_face[neighbour] > count_knew_face[neighbour]) {
				next_step_el.emplace(neighbour);  // ячейка была изменена, проверить ее готовность на следующем шаге
				count_knew_face[neighbour]++;  // всегда ли эта грань будет входящей (проверка по нормалям??)
			}
		}

	}

	return 0;
}

#endif // USE_OMP


int DEBUGFindCurCellWithHole(const std::set<IntId>& next_step_el, const std::vector<IntId>& count_in_face, std::vector<IntId>& count_knew_face,
	std::vector<IntId>& cur_el, const std::set<IntId>& inner_part, std::set<IntId>& outter_part,
	const std::map<IntId, FaceCell>& inner_cells, const std::vector<IntId>& all_pairs_id, const Vector3& direction, const std::vector<Normals>& normals) {

	cur_el.clear();

	const int N = next_step_el.size();
	for (auto cell : next_step_el)
	{
		if (outter_part.count(cell) != 0) {
			if (count_in_face[cell] == count_knew_face[cell] + 1) {  // граница не определена
				IntId try_id[3] = { -1,-1,-1 };

				FindIdCellInBoundary(direction, inner_part, inner_cells, normals, cell, try_id);

				if (try_id[0] == -1) continue;

				if (count_in_face[try_id[0]] == count_knew_face[try_id[0]]) { // если грань на другом конце определена

					cur_el.push_back(cell);
					outter_part.erase(cell);
					count_knew_face[cell]++;
					continue;
				}
			}
			else if (count_in_face[cell] == count_knew_face[cell]) {  // граница определена
				cur_el.push_back(cell);
				outter_part.erase(cell);
			}
		}
		else if (count_in_face[cell] == count_knew_face[cell]) {
			cur_el.push_back(cell);
		}
	}

	if (cur_el.size() == 0) {

		// плохо, но как есть. Если не смогли найти ни одну ячейку кандидата \
		попробовать пройти отдельно по внутренней границе, 

		std::list<IntId> buf_erase;

		for (auto cell : outter_part) {
			if (count_in_face[cell] == count_knew_face[cell] + 1) {  // граница не определена
				IntId try_id[3] = { -1,-1,-1 };
				FindIdCellInBoundary(direction, inner_part, inner_cells, normals, cell, try_id);
				if (try_id[0] == -1 ) continue;
				if (count_in_face[try_id[0]] == count_knew_face[try_id[0]]) { // если грань на другом конце определена

					cur_el.push_back(cell);
					count_knew_face[cell]++;
					buf_erase.push_back(cell);//outter_part.erase(cell);
					continue;
				}
			}
			else if (count_in_face[cell] == count_knew_face[cell]) {
				buf_erase.push_back(cell); //outter_part.erase(cell);
				cur_el.push_back(cell);
			}
		}

		for (auto el : buf_erase)
			outter_part.erase(el);

		if (cur_el.size() == 0) {
			std::cout << "NextCell is -1\n";
			return -1;
		}
	}

	return 0;
}