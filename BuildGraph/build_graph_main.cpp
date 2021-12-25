#include "build_graph_main.h"


int FractionInnerBoundary( const Vector3& direction, const std::vector<Normals>& normals, const std::map<IntId, FaceCell>& inner_cells, 
	const std::set<IntId>& full_boundary, std::set<IntId>& inner_part, std::set<IntId>& outter_part) {

	inner_part.clear();
	outter_part.clear();

	std::bitset<4> state;

	int bound_face_id = -1;
	for (auto el: full_boundary)
	{
		FindInAndOutFaces(direction, el, normals, state);

		bound_face_id = inner_cells.find(el)->second.face_id;

		if (state[bound_face_id % 4] != 0)  // выход€ща€ грань
			outter_part.emplace(el);
		else
			inner_part.emplace(el);

		


	}

	return 0;
}

int WriteFileGraph(const int i, const std::string& name_file_graph, const std::vector<IntId>& graph) {
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
	return 0;
}

#ifdef 0

std::vector<IntId> bound_vect;
#endif // 0


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
	writer->SetFileName(name_file_out.c_str());
	writer->SetInputData(ungrid);
	writer->Write();
	return 0;
}
#endif

#ifdef  0
int mainOld(int argc, char* argv[])
{ 
#ifdef  USE_MPI

	int np, myid;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

#endif //  USE_MPI

	std::string name_file_settings;
	if (argc <= 1)
		name_file_settings = "D:\\Desktop\\FilesCourse\\TestClaster\\settings_file.txt";
	else
		name_file_settings = argv[1];

	std::string name_file_vtk;
	std::string name_file_sphere_direction; // = main_file_direction + "SphereDir660.txt";
	std::string name_file_graph;

	std::string name_file_normals;
	std::string name_file_pairs;
	std::string name_file_inner_boundary;
	std::string name_file_init_boundary;
	std::string name_file_face_and_id;

	if (ReadStartSettings(name_file_settings, name_file_vtk, name_file_sphere_direction, name_file_graph, name_file_normals,
		name_file_pairs, name_file_inner_boundary,
		name_file_init_boundary, name_file_face_and_id)) MPI_RETURN(1);

#ifdef OnlyWriteFiles
#ifdef  USE_MPI
	if (myid == 0) 
#endif
	{
		double t = -omp_get_wtime();

		BuildSetForClaster(name_file_vtk, name_file_pairs, name_file_init_boundary, name_file_normals,
			name_file_inner_boundary, name_file_face_and_id);

		t += omp_get_wtime();
		std::cout << "Time for build start set: " << t << "\n";

		MPI_RETURN(0);
}
#endif 


	double tt;
#ifdef USE_MPI

	if (myid == 0) {
		std::cout << "OMP num_threads: " << omp_get_num_threads() << '\n';
		tt = -MPI_Wtime();
	}

#else
	std::cout << "OMP num_threads: " << omp_get_num_threads() << '\n';
	tt = -omp_get_wtime();
#endif
	

	std::vector<IntId> all_pairs_id;
	if (ReadPairsId(name_file_pairs, all_pairs_id)) MPI_RETURN(1);

	std::set<IntId> set_inner_boundary_faces;
	if (ReadInnerBoundary(name_file_inner_boundary, set_inner_boundary_faces)) MPI_RETURN(1);
#ifdef USE_MPI
	if(myid==0) 	
#endif // USE_MPI
		std::cout << "Inner boundary has " << set_inner_boundary_faces.size() << "faces\n";

	std::vector<Normals> normals;
	if (ReadNormalFile(name_file_normals, normals)) MPI_RETURN(1);

	std::map<IntId, Face> inner_faces;
	if (ReadInnerCellOfSphereAndId(name_file_face_and_id, inner_faces)) MPI_RETURN(1);


	std::vector<Type> directions;
	if (ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions)) MPI_RETURN(1);
	
	{
		/*std::ifstream ifile;

		ifile.open(name_file_sphere_direction);
		if (!ifile.is_open()) {
			std::cout << "Error read file sphere direction\n";
			return 1;
        }  
		int N = 0;
		ifile >> N;
		full_dir.resize(N);

		Type buf_s;
		for (int i = 0; i < N; i++) {
			ifile >> buf_s;
			ifile >> full_dir[i][0] >> full_dir[i][1] >> full_dir[i][2];
		}
		ifile.close();*/
	}

#ifdef USE_MPI
	if (myid == 0) {
		tt+= MPI_Wtime();
		std::cout << "Time reading in main procces: " << tt << "\n";
		tt = -MPI_Wtime();
	}
#else
	tt += omp_get_wtime();
	std::cout << "Time reading in main procces: " << tt << "\n";
	tt = -omp_get_wtime();

#endif // USE_MPI

	//int np, myid;
//#pragma omp parallel default(none) shared(bound_vect,all_pairs_id,set_inner_boundary_faces,normals,inner_faces,directions,name_file_graph,name_file_init_boundary, std::cout)  
	{
		std::vector<State> faces_state;
		std::set<IntId> set_boundary_cells;

		std::set<IntId> set_graph;

		Vector3 direction;
		bool flag = true;
		int count = 0;
#ifdef ONLY_ONE_DIRECTION
		for (int cur_direction = 0; cur_direction < 1; cur_direction++)
#else
#ifdef USE_MPI
		for (int cur_direction = myid; cur_direction < (directions.size() / 2); cur_direction += np)
#else
		for (int cur_direction = 0; cur_direction < (directions.size() / 2); cur_direction ++)
#endif
#endif //ONLY_ONE_DIRECTION
		{
			flag = true;
			set_graph.clear();
			ReadInitBoundarySet(name_file_init_boundary, set_boundary_cells);
			bound_vect.assign(set_boundary_cells.begin(), set_boundary_cells.end());
			InitFacesState(all_pairs_id, faces_state, set_inner_boundary_faces);
			FromSphericalToDecart(cur_direction, directions, direction);
			{
				
				//int face_state[4];
				//int i = 0;
				//int count = 0;
				//for (auto el : set_boundary_cells) {
				//	FindInAndOutFaces(direction, el, normals, face_state);
				//	for (int k = 0; k < 4; k++)
				//		if (face_state[k] == 1 && faces_state[el * 4 + k] == 1)
				//		{
				//			std::swap(bound_vect[i], bound_vect.front());//set_boundary_cells.erase(el);
				//			bound_vect.erase(bound_vect.begin() + i);
				//			i--;
				//			count++;
				//			break;
				//		}
				//	i++;
				//}
				////bound_vect.erase(bound_vect.begin(), bound_vect.begin() + count);
				//set_boundary_cells.clear();
				//set_boundary_cells.insert(bound_vect.begin(), bound_vect.end());
			}
			
		
			//direction = Vector3(1, 0, 0);
	
			std::ofstream ofile;

			ofile.open(name_file_graph + std::to_string(cur_direction) + ".txt");
			if (!ofile.is_open()) {
				std::cout << "Error open file\n";
				std::cout << "file_graph is not opened for writing\n";
				flag = false;//break;
			}


			while (set_boundary_cells.size() && flag) {
			//while (bound_vect.size() && flag) {
				IntId id_cell = GetNextCell(direction, set_boundary_cells, normals, faces_state, all_pairs_id, set_inner_boundary_faces, inner_faces);
				//IntId id_cell = GetNextCell2(direction, bound_vect, normals, faces_state, all_pairs_id, set_inner_boundary_faces, inner_faces);
				if (id_cell == -1) {

					ofile.close();
#ifdef USE_MPI
					std::cout << "Error proc: " << myid << '\n';
#endif // USE_MPI
					std::cout << "Error num_direction: " << cur_direction << '\n';
					std::cout << "Error.Complete " << set_graph.size() << " cells\n";
					flag = false; // break;

				}

				set_graph.emplace(id_cell);


				ofile << id_cell << '\n';
				ofile << std::fixed;


				if (flag)
					//ReBuildSetBondary2(id_cell, direction, bound_vect, normals, faces_state, all_pairs_id, set_graph);
					ReBuildSetBondary(id_cell, direction, set_boundary_cells, normals, faces_state, all_pairs_id, set_graph);

			//	if (++count % 5000 == 0 && myid == 0)
				//	printf("Direction %d, count= %d\n", cur_direction, count);

			}

			ofile.close();
#ifdef USE_MPI
			if(myid==0)
			printf("Id_proc: %d. Graph construction in the direction %d is completed\n", myid, cur_direction);
#else
			printf("Graph construction in the direction %d is completed\n", cur_direction);
#endif // USE_MPI

		}
	}
#ifdef USE_MPI
	if (myid == 0) {
		tt += MPI_Wtime();
		std::cout << "Full time: " << tt << '\n';
	}
#else
	tt += omp_get_wtime();
	std::cout << "Full time: " << tt << '\n';
#endif // USE_MPI	
	
#ifdef USE_VTK
#ifdef USE_MPI
	if(myid==0)
#endif
		WriteFileBoundary("D:\\Desktop\\FilesCourse\\MySphereGraph.vtk", name_file_graph + "0.txt", name_file_vtk);
#endif
	MPI_RETURN(0);
}
#endif //  0

int main(int argc, char** argv)
{
#ifdef  USE_MPI

	int np, myid;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

#endif //  USE_MPI

	std::string name_file_settings;
	if (argc <= 1)
		name_file_settings = "D:\\Desktop\\FilesCourse\\TestClaster\\settings_file.txt";
	else
		name_file_settings = argv[1];

	std::string name_file_vtk;
	std::string name_file_sphere_direction; // = main_file_direction + "SphereDir660.txt";
	std::string name_file_graph;

	std::string name_file_normals;
	std::string name_file_pairs;
	std::string name_file_inner_boundary;
	std::string name_file_init_boundary;
	std::string name_file_face_and_id;

	if (ReadStartSettings(name_file_settings, name_file_vtk, name_file_sphere_direction, name_file_graph, name_file_normals,
		name_file_pairs, name_file_inner_boundary,
		name_file_init_boundary, name_file_face_and_id)) MPI_RETURN(1);

#ifdef OnlyWriteFiles
#ifdef  USE_MPI
	if (myid == 0)
#endif
	{
		double t = -omp_get_wtime();

		BuildSetForClaster(name_file_vtk, name_file_pairs, name_file_init_boundary, name_file_normals,
			name_file_inner_boundary, name_file_face_and_id);

		t += omp_get_wtime();
		std::cout << "Time for build start set: " << t << "\n";

		MPI_RETURN(0);
	}
#endif 


	double tt;
#ifdef USE_MPI

	if (myid == 0) {
		std::cout << "OMP num_threads: " << omp_get_num_threads() << '\n';
		tt = -MPI_Wtime();
	}

#else
	std::cout << "OMP num_threads: " << omp_get_num_threads() << '\n';
	tt = -omp_get_wtime();
#endif

	std::vector<IntId> all_pairs_id;
	if (ReadPairsId(name_file_pairs, all_pairs_id)) MPI_RETURN(1);

	std::set<IntId> set_inner_boundary_cells;
	if (ReadInnerCellBoundary(name_file_inner_boundary, set_inner_boundary_cells)) MPI_RETURN(1);
	std::cout << "Inner boundary has " << set_inner_boundary_cells.size() << "faces\n";

	std::vector<Normals> normals;
	if (ReadNormalFile(name_file_normals, normals)) MPI_RETURN(1);

	std::map<IntId, FaceCell> inner_faces;
	if (ReadInnerCellOfSphereAndId(name_file_face_and_id, inner_faces)) MPI_RETURN(1);

	std::vector<Vector3> directions;
	if (ReadSphereDirectionDecart(name_file_sphere_direction, directions)) MPI_RETURN(1);

#ifdef USE_MPI
	if (myid == 0) {
		tt += MPI_Wtime();
		std::cout << "Time reading in main procces: " << tt << "\n";
		tt = -MPI_Wtime();
	}
#else
	tt += omp_get_wtime();
	std::cout << "Time reading in main procces: " << tt << "\n";
	tt = -omp_get_wtime();

#endif // USE_MPI


	{
		const int num_cells = all_pairs_id.size() / 4;

		std::vector<IntId> count_in_face(num_cells, 0);
		std::vector<IntId> count_knew_face(num_cells, 0);
		std::vector<IntId> graph(num_cells, 0);

		std::set<IntId> inner_part;
		std::set<IntId> outter_part;

		std::vector<State> faces_state;
		Vector3 direction;

		bool flag = true;
		int count = 0;

#ifdef USE_OMP
		std::vector<IntId> next_step_el_OMP;
		std::list<IntId> cur_el_OMP;
#else
		std::vector<IntId> cur_el;
		std::set<IntId> next_step_el;
#endif // USE_OMP


		const int count_cirection = directions.size();

#ifdef ONLY_ONE_DIRECTION
		for (int cur_direction = 0; cur_direction < 1; cur_direction++)
#else
#ifdef USE_MPI
		for (int cur_direction = myid; cur_direction < count_cirection; cur_direction += np)
#else
		for (int cur_direction = 0; cur_direction < count_cirection; ++cur_direction)
#endif
#endif //ONLY_ONE_DIRECTION

		{
			flag = true;
			InitFacesState(all_pairs_id, faces_state, inner_faces);
			direction = directions[cur_direction];

			FractionInnerBoundary(direction, normals, inner_faces, set_inner_boundary_cells, inner_part, outter_part);

			//-------------------------------------
#ifdef USE_OMP
			FindNumberOfAllInnerFaceAndKnew(direction, normals, faces_state, count_in_face, count_knew_face, next_step_el_OMP); 
#else
			FindNumberOfAllInnerFaceAndKnew(direction, normals, faces_state, count_in_face, count_knew_face, next_step_el);
#endif // USE_OMP
			//-------------------------------------

			int count_graph = 0; // число €чеек вошедших в граф 

#ifdef USE_OMP
			while (next_step_el_OMP.size() && flag) 
#else
			while (next_step_el.size() && flag) 
#endif // USE_OMP
			{ 

#ifdef USE_OMP

#ifdef GRID_WITH_INNER_BOUNDARY	
				IntId id_cell = FindCurCellWithHole(next_step_el_OMP, count_in_face, count_knew_face, cur_el_OMP, set_inner_boundary_cells, \
					inner_faces, all_pairs_id, direction, normals);
#else
				IntId id_cell = FindCurCell(next_step_el_OMP, count_in_face, count_knew_face, cur_el_OMP);
#endif //GRID_WITH_INNER_BOUNDARY

#else  //no use omp

#ifdef GRID_WITH_INNER_BOUNDARY	
				IntId id_cell = FindCurCellWithHole(next_step_el, count_in_face, count_knew_face, cur_el, inner_part, outter_part, \
					inner_faces, all_pairs_id, direction, normals);
#else
				IntId id_cell = FindCurCell(next_step_el, count_in_face, count_knew_face, cur_el);
#endif //GRID_WITH_INNER_BOUNDARY
#endif // USE_OMP

				if (id_cell == -1) {
#ifdef USE_MPI
					std::cout << "Error proc: " << myid << '\n';
#endif // USE_MPI
					WriteFileGraph(cur_direction, name_file_graph, graph);
					std::cout << "Error num_direction: " << cur_direction << '\n';
					std::cout << "Error.Complete " << count_graph << " cells\n";
					flag = false; // break;
				}

#ifdef USE_OMP
								
			   NewStep(all_pairs_id, count_in_face, count_knew_face, cur_el_OMP, next_step_el_OMP);
				
				for (auto el : cur_el_OMP) {
					graph[count_graph] = el;
					count_graph++;
				}
#else
			
				NewStep(all_pairs_id, count_in_face, count_knew_face, cur_el, next_step_el);
				for (auto el : cur_el) {
					graph[count_graph] = el;
					count_graph++;
				}

#endif // USE_OMP
				
				//	if (++count % 5000 == 0 && myid == 0)
				//	printf("Direction %d, count= %d\n", cur_direction, count);
			
			}//while
			
			if (WriteFileGraph(cur_direction, name_file_graph, graph))
				printf("Error writing graph file numbeb %d\n", cur_direction);

#ifdef USE_MPI
			if (myid == 0)
				printf("Id_proc: %d. Graph construction in the direction %d is completed\n", myid, cur_direction);
#else
			printf("Graph construction in the direction %d is completed\n", cur_direction);
#endif // USE_MPI

		}
	}

#ifdef USE_MPI
	if (myid == 0) {
		tt += MPI_Wtime();
		std::cout << "Full time: " << tt << '\n';
	}
#else
	tt += omp_get_wtime();
	std::cout << "Full time: " << tt << '\n';
#endif // USE_MPI	

#ifdef USE_VTK
#ifdef USE_MPI
	if (myid == 0)
#endif
		WriteFileBoundary("D:\\Desktop\\FilesCourse\\MySphereGraph.vtk", name_file_graph + "0.txt", name_file_vtk);
#endif
	MPI_RETURN(0);
}

