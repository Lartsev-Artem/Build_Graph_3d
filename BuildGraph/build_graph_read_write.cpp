#include "build_graph_read_write.h"

#ifdef OnlyWriteFiles

int ReadFileVtk(const std::string name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {
	vtkSmartPointer<vtkGenericDataObjectReader> reader_vtk =
		vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader_vtk->ReadAllScalarsOn();
	reader_vtk->SetFileName(name_file_vtk.c_str());
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

int ReadInnerBoundaryForWrite(const std::string name_file_boundary_inner, std::set<IntId>& id_inner_boundary_face) {

		std::ifstream ifile;

		ifile.open(name_file_boundary_inner);
		if (!ifile.is_open()) {
			std::cout << "Error read file boundary_inner\n";
			return 1;
		}

		id_inner_boundary_face.clear();

		IntId buf;
		ifile >> buf;
		std::cout << "Inner boundary has " << buf << "faces\n";
		while (ifile >> buf)
		{
			id_inner_boundary_face.emplace(buf);
		}

		ifile.close();
		return 0;
	}

int NormalAndSquareFace(size_t NumberCell, size_t NumberFace, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Eigen::Vector3d& n) {
	vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(NumberCell)->GetFace(NumberFace)->GetPointIds();

	Type P0[3], P1[3], P2[3];
	unstructuredgrid->GetPoint(idp->GetId(0), P0);
	unstructuredgrid->GetPoint(idp->GetId(1), P1);
	unstructuredgrid->GetPoint(idp->GetId(2), P2);

	Type a[3], b[3];
	for (size_t i = 0; i < 3; i++) {
		a[i] = P1[i] - P0[i];
		b[i] = P2[i] - P0[i];
	}
	n[0] = a[1] * b[2] - a[2] * b[1];
	n[1] = -a[0] * b[2] + a[2] * b[0];
	n[2] = a[0] * b[1] - a[1] * b[0];

	n.normalize();

	vtkSmartPointer<vtkIdList> idp2 = unstructuredgrid->GetCell(NumberCell)->GetPointIds();

	size_t id;
	for (size_t i = 0; i < 4; i++) {
		int count = 0;
		for (size_t j = 0; j < 3; j++)
			if (idp2->GetId(i) != idp->GetId(j))
				count++;
		if (count == 3) {
			id = i;
			break;
		}

	}

	Type sum = 0;
	Type P3[3];
	unstructuredgrid->GetPoint(idp2->GetId(id), P3);
	/*for (size_t i = 0; i < 3; i++){
		sum += n[i] * (P3[i] - P0[i]);
	}*/

	sum = P1[0] * (P2[1] - P3[1]) * P0[2] + P0[0] * (P3[1] - P2[1]) * P1[2] +
		P0[0] * (P1[1] - P3[1]) * P2[2] + P2[2] * (P1[0] * P3[1] - P1[0] * P0[1]) +
		P3[0] * (P0[2] * (P1[1] - P2[1]) + P1[2] * (P2[1] - P0[1]) + P2[2] * (P0[1] - P1[1]))
		+ P3[2] * (P1[0] * (P0[1] - P2[1]) + P0[0] * (P2[1] - P1[1])) +
		P2[0] * (P0[2] * (P3[1] - P1[1]) + P1[2] * (P0[1] - P3[1]) + P3[2] * (P1[1] - P0[1]));

	if (sum < 0)
		for (size_t i = 0; i < 3; i++)
			n[i] *= -1;
	return 0;
}

int WriteNormalFile(const std::string name_file_normals, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {
	std::ofstream ofile;

	ofile.open(name_file_normals);
	if (!ofile.is_open()) {
		std::cout << "Error write file normals\n";
		return 1;
	}

	int n = unstructured_grid->GetNumberOfCells();
	Eigen::Vector3d normal;
	ofile << n << '\n';
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			NormalAndSquareFace(i, j, unstructured_grid, normal);
			ofile << setprecision(16) << normal[0] << ' ' << normal[1] << ' ' << normal[2] << '\n';
		}

	}
	ofile.close();
	return 0;
}

int WritePairsId(const std::string name_file_pairs, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {

	auto GetNumberNeighborFace{ [](const int a, const int b, const int c, vtkCell* neighbor_cell)
		{

			vtkIdList* idc;

			int x, y, z;
			for (int i = 0; i < 4; i++)
			{
				idc = neighbor_cell->GetFace(i)->GetPointIds();
				x = idc->GetId(0);
				y = idc->GetId(1);
				z = idc->GetId(2);

				if (a == x && b == y && c == z) return i;
				else if (a == x && b == z && c == y) return i;
				else if (a == y && b == x && c == z) return i;
				else if (a == y && b == z && c == x) return i;
				else if (a == z && b == x && c == y) return i;
				else if (a == z && b == y && c == x) return i;

			}
			return -2;
		} };

	int count_unique_face = 0;
	const int N = unstructured_grid->GetNumberOfCells();
	std::vector<int>all_pairs_face(N * 4);
	for (int i = 0; i < N * 4; i++)
		all_pairs_face[i] = -2;

	vtkSmartPointer<vtkIdList> idp = vtkSmartPointer< vtkIdList>::New();
	vtkSmartPointer<vtkIdList> idc = vtkSmartPointer< vtkIdList>::New();

	int id_a, id_b, id_c;
	for (vtkIdType num_cell = 0; num_cell < N; ++num_cell) {

		for (int num_face = 0; num_face < 4; ++num_face) {
			if (all_pairs_face[num_cell * 4 + num_face] != -2) continue;
			++count_unique_face;

			idp = unstructured_grid->GetCell(num_cell)->GetFace(num_face)->GetPointIds();
			id_a = idp->GetId(0);
			id_b = idp->GetId(1);
			id_c = idp->GetId(2);

			/*ћожет быть проблема с указател€ми на списки!*/
			unstructured_grid->GetCellNeighbors(num_cell, idp, idc);
			int face = num_cell * 4 + num_face;

			if (idc->GetNumberOfIds() == 1) {
				int id_neighbor_cell = idc->GetId(0);
				int id_neighbor_face = GetNumberNeighborFace(id_a, id_b, id_c, unstructured_grid->GetCell(id_neighbor_cell));
				
				all_pairs_face[face] = id_neighbor_cell * 4 + id_neighbor_face;
				all_pairs_face[id_neighbor_cell * 4 + id_neighbor_face] = face;
			}
			else if (idc->GetNumberOfIds() == 0)
				all_pairs_face[face] = -1; // гранична€ €чейка
			else
				std::cout << "More than 1 neighbor????\n";
		}

	}


	std::ofstream ofile;

	ofile.open(name_file_pairs);
	if (!ofile.is_open()) {
		std::cout << "Error write file pairs\n";
		return 1;
	}

	ofile << all_pairs_face.size() << '\n';
	for (size_t i = 0; i < all_pairs_face.size(); i++)
	{
		ofile << all_pairs_face[i] << '\n';
	}
	ofile.close();
	return 0;
}

int WriteInitBoundarySetAndInnerBoundaryFace(const std::string name_file_boundary, const std::string name_file_boundary_inner,
	const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {
	
	std::set<IntId> boundary_cells;
	std::set<IntId> inner_boundary_faces;

	boundary_cells.clear();
	inner_boundary_faces.clear();

	int N = unstructured_grid->GetNumberOfCells();

	vtkSmartPointer<vtkIdList> idc = vtkSmartPointer<vtkIdList>::New();

	for (size_t i = 0; i < N; ++i) {

		for (size_t j = 0; j < 4; ++j) {
			unstructured_grid->GetCellNeighbors(i, unstructured_grid->GetCell(i)->GetFace(j)->GetPointIds(), idc);
			if (idc->GetNumberOfIds() == 0) {

				Vector3 P(unstructured_grid->GetCell(i)->GetFace(j)->GetPoints()->GetPoint(0));
				if ((P-center_point).norm() > R) { // внешн€€ сфера
					boundary_cells.emplace(i);
				}
				else {
					inner_boundary_faces.emplace(i * 4 + j);
					boundary_cells.emplace(i);
				}
				break;
			}
		}
	}


	std::ofstream ofile;

	ofile.open(name_file_boundary);
	if (!ofile.is_open()) {
		std::cout << "Error write file boundary\n";
		return 1;
	}

	ofile << boundary_cells.size() << '\n';
	for (auto el : boundary_cells)
		ofile << el << '\n';
	ofile.close();


	ofile.open(name_file_boundary_inner);
	if (!ofile.is_open()) {
		std::cout << "Error write file boundary_inner\n";
		return 1;
	}

	ofile << inner_boundary_faces.size() << '\n';
	for (auto el : inner_boundary_faces)
		ofile << el << '\n';
	ofile.close();

	return 0;
}

int WriteInnerCellOfSphere(const std::string name_file_inner_sphere, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid,
	const std::string name_file_boundary_inner) {

	

	std::set<IntId> id_inner_boundary_face; // хранит номера граней в формате num_cell*4+num_face

	if (ReadInnerBoundaryForWrite(name_file_boundary_inner, id_inner_boundary_face)) return 1;


	std::ofstream ofile;

	ofile.open(name_file_inner_sphere);
	if (!ofile.is_open()) {
		std::cout << "Error write inner_sphere\n";
		return 1;
	}

	Type A[3];
	Type B[3];
	Type C[3];
	vtkPoints* points_face;
	
	ofile << id_inner_boundary_face.size() << '\n';

	for (auto el : id_inner_boundary_face) {
		points_face = unstructured_grid->GetCell(el / 4)->GetFace(el % 4)->GetPoints();

		points_face->GetPoint(0, A);
		points_face->GetPoint(0, B);
		points_face->GetPoint(0, C);

		ofile << setprecision(16) << A[0] << ' ' << A[1] << ' ' << A[2]
		      << B[0] << ' ' << B[1] << ' ' << B[2]
		      << C[0] << ' ' << C[1] << ' ' << C[2] << '\n';
	}

	ofile.close();
	return 0;
}

int WriteInnerCellAndIdFaces(const std::string name_file_boundary_inner, const std::string name_file_face_and_id,
	const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {

	std::set<IntId> id_inner_boundary_face; // хранит номера граней в формате num_cell*4+num_face
	if (ReadInnerBoundaryForWrite(name_file_boundary_inner, id_inner_boundary_face)) return 1;

	std::ofstream ofile;
	ofile.open(name_file_face_and_id);
	if (!ofile.is_open()) {
		std::cout << "Error write inner_sphere\n";
		return 1;
	}

	Type A[3];
	Type B[3];
	Type C[3];
	vtkPoints* points_face;

	ofile << id_inner_boundary_face.size() << '\n';

	for (auto el : id_inner_boundary_face) {
		points_face = unstructured_grid->GetCell(el / 4)->GetFace(el % 4)->GetPoints();

		points_face->GetPoint(0, A);
		points_face->GetPoint(1, B);
		points_face->GetPoint(2, C);

		ofile<<setprecision(16) << el << ' ' << A[0] << ' ' << A[1] << ' ' << A[2] << ' '
			               << B[0] << ' ' << B[1] << ' ' << B[2] << ' '
		                   << C[0] << ' ' << C[1] << ' ' << C[2] << '\n';
	}

	ofile.close();

	return 0;
}

int BuildSetForClaster(const std::string name_file_vtk, const std::string name_file_pairs,
	const std::string name_file_boundary, const std::string name_file_normals, const std::string name_file_boundary_inner, 
	const std::string name_file_face_and_id) {

	std::cout << "Start build \n";

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	if (ReadFileVtk(name_file_vtk.c_str(), unstructured_grid)) {
		std::cout << "Error reading file vtk\n";
		return 1;
	}
	
	if (WritePairsId(name_file_pairs, unstructured_grid)) return 1;

	if (WriteInitBoundarySetAndInnerBoundaryFace(name_file_boundary, name_file_boundary_inner, unstructured_grid)) return 1;

	if (WriteNormalFile(name_file_normals, unstructured_grid)) return 1;

	//if (WriteInnerCellOfSphere(name_file_inner_sphere, unstructured_grid, name_file_boundary_inner)) return 1;

	if (WriteInnerCellAndIdFaces(name_file_boundary_inner, name_file_face_and_id, unstructured_grid)) return 1;

	std::cout << "Build end\n";
	return 0;
}

#endif

#ifdef OnlyReadFiles

int ReadPairsId(const std::string name_file_pairs, std::vector<IntId>& all_pairs_face) {

	std::ifstream ifile;

	ifile.open(name_file_pairs);
	if (!ifile.is_open()) {
		std::cout << "Error read file pairs\n";
		return 1;
	}

	int N;
	ifile >> N;
	all_pairs_face.resize(N);
	for (int i = 0; i < N; ++i)
	{
		ifile >> all_pairs_face[i];
	}
	ifile.close();
	return 0;
}

int ReadInitBoundarySet(const std::string name_file_boundary, std::set<IntId>& boundary_cells) {

	std::ifstream ifile;

	ifile.open(name_file_boundary);
	if (!ifile.is_open()) {
		std::cout << "Error read file boundary\n";
		return 1;
	}

	IntId N;
	ifile >> N;
	//std::cout << "Boundary has " << N << "cells\n";
	while (ifile >> N) {
		boundary_cells.emplace(N);
	}

	ifile.close();
	return 0;
}

int ReadInnerBoundary(const std::string name_file_boundary_inner, std::set<IntId>& id_inner_boundary_face) {

	std::ifstream ifile;

	ifile.open(name_file_boundary_inner);
	if (!ifile.is_open()) {
		std::cout << "Error read file boundary_inner\n";
		return 1;
	}

	id_inner_boundary_face.clear();

	IntId buf;
	ifile >> buf;
//	std::cout << "Inner boundary has " << buf << "faces\n";
	while (ifile >> buf)
	{
		id_inner_boundary_face.emplace(buf);
	}

	ifile.close();
	return 0;
}

int ReadInnerCellBoundary(const std::string name_file_boundary_inner, std::set<IntId>& id_inner_boundary_cell) {

	std::ifstream ifile;

	ifile.open(name_file_boundary_inner);
	if (!ifile.is_open()) {
		std::cout << "Error read file boundary_inner\n";
		return 1;
	}

	id_inner_boundary_cell.clear();

	IntId buf;
	ifile >> buf;
	//	std::cout << "Inner boundary has " << buf << "faces\n";
	while (ifile >> buf)
	{
		id_inner_boundary_cell.emplace(buf/4);
	}

	ifile.close();
	return 0;
}

int ReadNormalFile(const std::string name_file_normals, std::vector<Normals>& normals) {
	std::ifstream ifile;

	ifile.open(name_file_normals);
	if (!ifile.is_open()) {
		std::cout << "Error read file normals\n";
		return 1;
	}

	int N;
	ifile >> N;
	normals.resize(N);

	Normals norm(4);
	for (size_t i = 0; i < N; i++)
	{
		for (int j = 0; j < 4; j++)
			ifile >> norm.n[j][0] >> norm.n[j][1] >> norm.n[j][2];
		normals[i] = norm;
	}

	ifile.close();
	return 0;
}

int ReadInnerCellOfSphere(const std::string name_file_inner_sphere, std::vector<Face>& inner_faces) {

	std::ifstream ifile;

	ifile.open(name_file_inner_sphere);
	if (!ifile.is_open()) {
		std::cout << "Error read inner_sphere\n";
		return 1;
	}

	int N;
	ifile >> N;

	inner_faces.resize(N);

	Face face;
	for (int i = 0; i < N; ++i){

		ifile >> face.A[0] >> face.A[1] >> face.A[2]
			>> face.B[0] >> face.B[1] >> face.B[2]
			>> face.C[0] >> face.C[1] >> face.C[2];
		inner_faces[i] = face;
	}

	ifile.close();
	return 0;
}

int ReadInnerCellOfSphereAndId(const std::string name_file_face_and_id, std::map<IntId,Face>& inner_faces) {

	std::ifstream ifile;

	ifile.open(name_file_face_and_id);
	if (!ifile.is_open()) {
		std::cout << "Error read inner_sphere\n";
		return 1;
	}

	int N;
	ifile >> N;

	Face face;
	IntId id;
	for (int i = 0; i < N; ++i) {
		ifile >> id;
		ifile >> face.A[0] >> face.A[1] >> face.A[2]
			>> face.B[0] >> face.B[1] >> face.B[2]
			>> face.C[0] >> face.C[1] >> face.C[2];
		inner_faces.emplace(id, face);
	}

	ifile.close();
	return 0;
}

int ReadInnerCellOfSphereAndId(const std::string name_file_face_and_id, std::map<IntId, FaceCell>& inner_cells) {

	std::ifstream ifile;

	ifile.open(name_file_face_and_id);
	if (!ifile.is_open()) {
		std::cout << "Error read inner_sphere\n";
		return 1;
	}

	int N;
	ifile >> N;

	Face face;
	IntId id;
	for (int i = 0; i < N; ++i) {
		ifile >> id;
		ifile >> face.A[0] >> face.A[1] >> face.A[2]
			>> face.B[0] >> face.B[1] >> face.B[2]
			>> face.C[0] >> face.C[1] >> face.C[2];
		inner_cells.emplace(id / 4, FaceCell(id, face));
	}

	ifile.close();
	return 0;
}

int ReadSphereDirectionDecartToSpherical(const std::string name_file_sphere_direction, std::vector<Type>& directions_all) {

	std::ifstream ifile;

	ifile.open(name_file_sphere_direction);
	if (!ifile.is_open()) {
		std::cout << "Error read file sphere direction\n";
		return 1;
	}
	int N = 0;
	ifile >> N;
	directions_all.resize(2 * N);

	Type buf_s;
	int i = 0;
	Type P[3];
	Type theta;
	Type fi;
	for (int i = 0; i < N; i++) {
		ifile >> buf_s;
		ifile >> P[0] >> P[1] >> P[2];
		FromDecartToSphere(P, fi, theta);
		directions_all[i] = theta;
		directions_all[N + i] = fi;
	}
	ifile >> buf_s;
	ifile.close();
	return 0;
}

int ReadSphereDirectionDecart(const std::string name_file_sphere_direction, std::vector<Vector3>& directions_all) {

	std::ifstream ifile;

	ifile.open(name_file_sphere_direction);
	if (!ifile.is_open()) {
		std::cout << "Error read file sphere direction\n";
		return 1;
	}
	int N = 0;
	ifile >> N;
	directions_all.resize( N);

	Type buf_s;
	int i = 0;

	for (int i = 0; i < N; i++) {
		ifile >> buf_s;
		ifile >> directions_all[i][0] >> directions_all[i][1] >> directions_all[i][2];
	}
	ifile >> buf_s;
	ifile.close();
	return 0;
}

int WriteFileGraph(const int i, const std::string& name_file_graph, const std::vector<IntId>& graph) {

#ifdef FastWriteFile


	std::unique_ptr<FILE, int(*)(FILE*)> file_graph(fopen((name_file_graph + std::to_string(i) + ".bin").c_str(), "wb"), fclose);
	if (!file_graph) { printf("file_graph is not opened for writing\n"); return 1; }

	const int n = graph.size();
	fwrite_unlocked(graph.data(), sizeof(IntId), n, file_graph.get());

	fclose(file_graph.get());


	std::unique_ptr<FILE, int(*)(FILE*)> file_id(fopen((std::string(BASE_ADRESS) + "id_defining_faces" + std::to_string(i) + ".bin").c_str(), "wb"), fclose);
	if (!file_id) { printf("file_id is not opened for writing\n"); return 1; }

	int size = id_try_surface.size();
	fwrite_unlocked(&size, sizeof(int), 1, file_id.get());
	fwrite_unlocked(id_try_surface.data(), sizeof(IntId), size, file_id.get());

	fclose(file_id.get());


	std::unique_ptr<FILE, int(*)(FILE*)> file_dist(fopen((std::string(BASE_ADRESS) + "dist_defining_faces" + std::to_string(i) + ".bin").c_str(), "wb"), fclose);
	if (!file_dist) { printf("file_dist is not opened for writing\n"); return 1; }

	size = dist_try_surface.size();
	fwrite_unlocked(&size, sizeof(int), 1, file_dist.get());
	fwrite_unlocked(dist_try_surface.data(), sizeof(Type), size, file_dist.get());

	fclose(file_dist.get());

	std::unique_ptr<FILE, int(*)(FILE*)> file_x(fopen((std::string(BASE_ADRESS) + "x_defining_faces" + std::to_string(i) + ".bin").c_str(), "wb"), fclose);
	if (!file_x) { printf("file_x is not opened for writing\n"); return 1; }

	size = x_try_surface.size();
	fwrite_unlocked(&size, sizeof(int), 1, file_x.get());
	fwrite_unlocked(x_try_surface.data(), sizeof(Vector3), size, file_x.get());

	fclose(file_x.get());

#else

	std::ofstream ofile;
	ofile.open(name_file_graph + std::to_string(i) + ".txt");
	if (!ofile.is_open()) {
		std::cout << "Error open file\n";
		std::cout << "file_graph is not opened for writing\n";
		return 1;
	}

	for (auto el : graph)
		ofile << el << '\n';

	ofile << std::fixed;
	ofile.close();

	ofile.open(std::string(BASE_ADRESS) + "id_defining_faces" + std::to_string(i) + ".txt");
	if (!ofile.is_open()) {
		std::cout << "Error open file\n";
		std::cout << "file id_defining_faces is not opened for writing\n";
		return 1;
	}

	ofile << id_try_surface.size() << '\n';
	for (int i = 0; i < id_try_surface.size(); i += 3)
		ofile << id_try_surface[i] << ' ' << id_try_surface[i + 1] << ' ' << id_try_surface[i + 2] << ' ' << '\n';

	ofile << std::fixed;
	ofile.close();

	ofile.open(std::string(BASE_ADRESS) + "dist_defining_faces" + std::to_string(i) + ".txt");
	if (!ofile.is_open()) {
		std::cout << "Error open file\n";
		std::cout << "file dist_defining_faces is not opened for writing\n";
		return 1;
	}

	ofile << dist_try_surface.size() << '\n';
	for (int i = 0; i < dist_try_surface.size(); i += 3)
		ofile << dist_try_surface[i] << ' ' << dist_try_surface[i + 1] << ' ' << dist_try_surface[i + 2] << ' ' << '\n';

	ofile << std::fixed;
	ofile.close();

#endif // FastWriteFile

	return 0;
}

int WriteFileGraph(std::unique_ptr<FILE, int(*)(FILE*)>& file_graph, std::unique_ptr<FILE, int(*)(FILE*)>& file_id,
	std::unique_ptr<FILE, int(*)(FILE*)>& file_dist, std::unique_ptr<FILE, int(*)(FILE*)>& file_x,
	const int i, const int n, const std::vector<IntId>& graph) {

	fwrite_unlocked(graph.data(), sizeof(IntId), n, file_graph.get());

	id_try_size += id_try_surface.size();	
	fwrite_unlocked(id_try_surface.data(), sizeof(IntId), id_try_surface.size(), file_id.get());


	dist_try_size += dist_try_surface.size();	
	fwrite_unlocked(dist_try_surface.data(), sizeof(Type), dist_try_surface.size(), file_dist.get());


	x_try_size += x_try_surface.size();	
	fwrite_unlocked(x_try_surface.data(), sizeof(Vector3), x_try_surface.size(), file_x.get());

	return 0;
}

#ifdef USE_VTK
int WriteFileBoundary(const std::string name_file_out, const std::string name_file_graph, const std::string name_file_grid) {
	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	if (ReadFileVtk(name_file_grid.c_str(), unstructured_grid))  return 1;

	int n = unstructured_grid->GetNumberOfCells();

	vtkSmartPointer<vtkIntArray> bound_array =
		vtkSmartPointer<vtkIntArray>::New();

	bound_array->SetNumberOfTuples(n);

	std::ifstream ifile;
	ifile.open(name_file_graph);
	if (!ifile.is_open()) {
		std::cout << "Error open file\n";
		std::cout << "file_graph is not opened for reading\n";
		return 1;
	}

	for (int i = 0; i < n; i++)
		bound_array->SetTuple1(i, 0);

	int i = 0, el;
	ifile >> el;
	bound_array->SetTuple1(el, -100);
	while (ifile >> el)
		bound_array->SetTuple1(el, i++);
	ifile.close();

	vtkSmartPointer<vtkUnstructuredGrid> ungrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	ungrid = unstructured_grid;
	ungrid->GetCellData()->SetActiveScalars("energy");
	ungrid->GetCellData()->SetScalars(bound_array);

	vtkSmartPointer<vtkGenericDataObjectWriter> writer =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileTypeToBinary();
	writer->SetFileName(name_file_out.c_str());
	writer->SetInputData(ungrid);
	writer->Write();
	return 0;
}
#endif

#endif