#ifndef BUILD_GRAPH_PRJ_CONFIG
#define BUILD_GRAPH_PRJ_CONFIG

// ����� ���� ������� ���. ����� ������� ���� ������ ����� ������
#define OnlyWriteFiles  // ��������� ������ vtk � �������� ������ ����������� ��� ��������� ������
#define OnlyReadFiles     // ������ �������������� ������ ��� ���������� ������ � �� ���������� ������

//#define USE_VTK          // ������������� vtk ��� ������ ����������� � ���� �����
//#define USE_MPI        // ����������� ���������� mpi
//#define USE_OMP        // ����������� ���������� omp (�������������� / ������ mpi) 

//#define GRID_WITH_INNER_BOUNDARY  // ���� ��� ����� � ���������� �������� (��� ����������� ������� ���������, ��� ������������� �������� �����)

// DEBUG
//#define ONLY_ONE_DIRECTION  // ������  � ���������������� ������ 

#include <iomanip>
#include <eigen3/Eigen/Dense>       // ��� ����������� ���������� ����� � ��������

typedef Eigen::Vector3d Vector3;
typedef double Type;
typedef int IntId;
typedef uint8_t State;

const Vector3 center_point(0, 0, 0);//center_point(1, 0, 0);//center_point(0, 0, 0);
const Type R = 0.51;//R = 0.11; //0.51; // ������ ���������� ����� (� �������)

//const Vector3 center_point(1, 0, 0);
//const Type R = 0.11;  // ������ ���������� ����� (� �������)

#endif