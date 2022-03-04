#ifndef BUILD_GRAPH_PRJ_CONFIG
#define BUILD_GRAPH_PRJ_CONFIG

// могут быть активны оба. Тогда сначала идет запись потои чтение
#define OnlyWriteFiles  // включение файлов vtk и создание файлов необходимых для постоения графов
#define OnlyReadFiles     // чтение сформированных файлов для построения графов и из дальнейший расчёт

#define FastWriteFile

#define USE_VTK          // использование vtk для вывода результатов в виде сетки
//#define USE_MPI        // подключение технологии mpi
//#define USE_OMP        // подключение технологии omp (самостоятельно / вместе mpi) 

#define GRID_WITH_INNER_BOUNDARY  // граф для сетки с внутренней границей (для оптимизации следует отключить, при использование сплошной сетки)

// DEBUG
//#define ONLY_ONE_DIRECTION  // только  в последовательном режиме 

#include <iomanip>
#include <eigen3/Eigen/Dense>       // для определения внутренней сферы и нормалей

typedef Eigen::Vector3d Vector3;
typedef double Type;
typedef int IntId;
typedef uint8_t State;

extern std::string BASE_ADRESS;// "D:\\Desktop\\FilesCourse\\IllumGrid\\"

//const Vector3 center_point(0, 0, 0);
//const Type R = 0.51; // радиус внутренней сферы (с запасом)

const Vector3 center_point(1, 0, 0);
const Type R = 0.11;  // радиус внутренней сферы (с запасом)

#endif