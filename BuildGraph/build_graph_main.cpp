#include "build_graph_main.h"

std::string BASE_ADRESS;

std::vector<int> id_try_surface;
std::vector<Type> dist_try_surface;
std::vector<Vector3> x_try_surface;
TrySolve buf_try;

int id_try_size;
int dist_try_size;
int x_try_size;

#ifdef USE_VTK
int WriteFile(const std::string name_file_out, const std::string name_file_graph, const std::string name_file_grid) {
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
		bound_array->SetTuple1(el/4, i++);
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
int WriteFile(const std::string name_file_out, const std::set<int> bound, const std::string name_file_grid) {
	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	if (ReadFileVtk(name_file_grid.c_str(), unstructured_grid))  return 1;

	int n = unstructured_grid->GetNumberOfCells();

	vtkSmartPointer<vtkIntArray> bound_array =
		vtkSmartPointer<vtkIntArray>::New();

	bound_array->SetNumberOfTuples(n);


	for (int i = 0; i < n; i++)
		bound_array->SetTuple1(i, 0);

	for (auto el: bound)
	{
		bound_array->SetTuple1(el, 1);
	}

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
#endif // USE_VTK

int main(int argc, char** argv)
{
	id_try_size = 0;
	dist_try_size = 0;
	x_try_size = 0;

#ifdef  USE_MPI

	int np, myid;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

#endif //  USE_MPI

	std::string name_file_settings;
	if (argc <= 1)
		name_file_settings = "D:\\Desktop\\FilesCourse\\settings_file_build.txt";
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

	if (ReadStartSettings(name_file_settings, name_file_vtk, name_file_sphere_direction, BASE_ADRESS, name_file_normals,
		name_file_pairs, name_file_inner_boundary,
		name_file_init_boundary, name_file_face_and_id)) MPI_RETURN(1);

	// ????? ?????? ???? ifndef
#ifdef OnlyWriteFiles
#ifdef  USE_MPI
	if (myid == 0)
#endif
	{
		double t = -omp_get_wtime();

		BuildSetForClaster(name_file_vtk, name_file_pairs, name_file_init_boundary, name_file_normals,
			name_file_inner_boundary, name_file_face_and_id, BASE_ADRESS+"squares.bin", BASE_ADRESS+"volume.bin", BASE_ADRESS + "centers.bin");

		t += omp_get_wtime();
		std::cout << "Time for build start set: " << t << "\n";

#ifndef OnlyReadFiles
		MPI_RETURN(0);
#endif
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

		std::ofstream ofile(std::string(BASE_ADRESS) + "Logs.txt");

		std::unique_ptr<FILE, int(*)(FILE*)> file_graph(fopen((std::string(BASE_ADRESS)+"graph" + ".bin").c_str(), "wb"), fclose);
		if (!file_graph) { printf("file_graph is not opened for writing\n"); return 1; }

		std::unique_ptr<FILE, int(*)(FILE*)> file_id(fopen((std::string(BASE_ADRESS) + "id_defining_faces" + ".bin").c_str(), "wb"), fclose);
		if (!file_id) { printf("file_id is not opened for writing\n"); return 1; }

		std::unique_ptr<FILE, int(*)(FILE*)> file_dist(fopen((std::string(BASE_ADRESS) + "dist_defining_faces" + ".bin").c_str(), "wb"), fclose);
		if (!file_dist) { printf("file_dist is not opened for writing\n"); return 1; }


		std::unique_ptr<FILE, int(*)(FILE*)> file_x(fopen((std::string(BASE_ADRESS) + "x_defining_faces" +  ".bin").c_str(), "wb"), fclose);
		if (!file_x) { printf("file_x is not opened for writing\n"); return 1; }

		std::unique_ptr<FILE, int(*)(FILE*)> file_shift_try(fopen((std::string(BASE_ADRESS) + "ShiftTry" + ".bin").c_str(), "wb"), fclose);
		if (!file_x) { printf("file_x is not opened for writing\n"); return 1; }
		

		const int count_cirection = directions.size();
		std::vector<int> ArrayData(num_cells, 0);

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
			direction = directions[cur_direction];  // 13 -- error direction for test.vtk grid
			//direction = Vector3(1, 0, 0);

			//direction = directions[3];// Vector3(1, 0, 0);

			FractionInnerBoundary(direction, normals, inner_faces, set_inner_boundary_cells, inner_part, outter_part);

			id_try_surface.clear();
			dist_try_surface.clear();
			x_try_surface.clear();

			x_try_surface.reserve(outter_part.size() * 3);
			id_try_surface.reserve(outter_part.size() * 3);
			dist_try_surface.reserve(outter_part.size() * 3);

			//-------------------------------------
#ifdef USE_OMP
			FindNumberOfAllInnerFaceAndKnew(direction, normals, faces_state, count_in_face, count_knew_face, next_step_el_OMP); 
#else
			FindNumberOfAllInnerFaceAndKnew(direction, normals, faces_state, count_in_face, count_knew_face, next_step_el);
#endif // USE_OMP
			//-------------------------------------

			int count_graph = 0; // ????? ????? ???????? ? ???? 
			graph.assign(num_cells, -1); // ??? ???????????? ????????? ???????????
			bool try_restart = true;

#ifdef USE_OMP
			while (next_step_el_OMP.size() && flag) 
#else
			while (next_step_el.size() && flag) 
#endif // USE_OMP
			{ 

#ifdef USE_OMP

#ifdef GRID_WITH_INNER_BOUNDARY	
				IntId id_cell = FindCurCellWithHole(next_step_el_OMP, count_in_face, count_knew_face, cur_el_OMP, inner_part, outter_part, \
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
					//WriteFileGraph(cur_direction, name_file_graph, graph);
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
			
				static int ccoun = 1;
				{
					for (auto el : cur_el)
						ArrayData[el] = ccoun;
				}
				ccoun++;


				NewStep(all_pairs_id, count_in_face, count_knew_face, cur_el, next_step_el);
				for (auto el : cur_el) {
					graph[count_graph] = el;
					count_graph++;
				}


#endif // USE_OMP
				
				//	if (++count % 5000 == 0 && myid == 0)
				//	printf("Direction %d, count= %d\n", cur_direction, count);

				
				if (count_graph < graph.size() && next_step_el.size() == 0 && try_restart) {
					try_restart = !try_restart;

					flag = true;
					//printf("Size outter boundary %d\n", outter_part.size());
#ifdef GRID_WITH_INNER_BOUNDARY
					if (outter_part.size() != 0) 
					{
						cur_el.clear();
						next_step_el.clear();
						next_step_el.insert(outter_part.begin(), outter_part.end());
						printf("Error. try_short_restart %d\n", cur_direction);
						ofile << "Error. try_short_restart" << cur_direction << "\n";

						//???????? ?????????? ?????????? ??????? ?? ??????????
						static int cc = 0;
						if (cc++ < 3)
							try_restart = !try_restart;

						continue;
					}
#endif //GRID_WITH_INNER_BOUNDARY

					printf("\n        Error. try_restart %d            \n", cur_direction);
					ofile << "----------   Error. try_restart   ----------" << cur_direction << "\n";
					
					cur_el.clear();
					next_step_el.clear();
					
					cout << outter_part.size() << '\n';
					for (size_t i = 0; i < num_cells; i++)
					{
						if (count_in_face[i] > count_knew_face[i]) {
							next_step_el.emplace(i);  // ?????? ???? ????????, ????????? ?? ?????????? ?? ????????? ????
							
						}

					}

					//NewStep(all_pairs_id, count_in_face, count_knew_face, cur_el, next_step_el);
				}
			
			}//while

			try_restart = !try_restart;
			if (count_graph < graph.size()) {
				printf("-------------------------Error size graph-------------------------------\n");
				ofile << "----------  Error size graph   ----------(" << cur_direction <<") Size: "<<count_graph << "\n";
			}

			fwrite_unlocked(&id_try_size, sizeof(int), 1, file_shift_try.get());

			if (WriteFileGraph(file_graph, file_id, file_dist, file_x, cur_direction, num_cells, graph))
				printf("Error writing graph file numbeb %d\n", cur_direction);
			
			

#ifdef USE_MPI
			if (myid == 0)
				printf("Id_proc: %d. Graph construction in the direction %d is completed\n", myid, cur_direction);
#else
			printf("Graph construction in the direction %d is completed\n", cur_direction);
#endif // USE_MPI

		}

		std::ofstream steps(BASE_ADRESS + "steps.txt");
		for (size_t i = 0; i < ArrayData.size(); i++)
		{
			steps << ArrayData[i] << '\n';
		}

		steps.close();

		ofile.close();
		fclose(file_graph.get());
		fclose(file_id.get());
		fclose(file_dist.get());
		fclose(file_x.get());
		fclose(file_shift_try.get());
		

		std::ofstream ofile2(std::string(BASE_ADRESS) + "Size.txt");
		ofile2 << id_try_size << '\n'; // "id_try_size: " << id_try_size;
		/*ofile2 << "\ndist_try_size: " << dist_try_size;
		ofile2 << "\nx_try_size: " << x_try_size << '\n';*/
		ofile2.close();

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
	{
		WriteFileBoundary("D:\\Desktop\\FilesCourse\\Steps.vtk", BASE_ADRESS + "steps.txt", name_file_vtk);
		//WriteFileBoundary("D:\\Desktop\\FilesCourse\\MySphereGraph.vtk", name_file_graph + "0.txt", name_file_vtk);
	}
#endif
	MPI_RETURN(0);
}

