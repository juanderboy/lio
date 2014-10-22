/* includes */
#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

#include "common.h"
#include "init.h"
#include "partition.h"

using namespace std;
using namespace G2G;

/************************************************************
 * Construct partition
 ************************************************************/

//Sorting the cubes in increasing order of size in bytes in GPU.
template <typename T>
bool comparison_by_size(const T & a, const T & b) {
    return a.size_in_gpu() < b.size_in_gpu();
}
template <typename T>
void sortBySize(std::vector<T> input) {
    sort(input.begin(), input.end(), comparison_by_size<T>);
}

/* methods */
namespace G2G {
  template <class T> unsigned long do_benchmark();
}

void Partition::generate_gpu_profile() {
    int total_threads = 1;
    int gpu_count = 0;
    int prev_gpu = 0;
    #if CPU_KERNELS
    return;
    #endif

    #ifndef _OPENMP
    total_threads = 1;
    #endif
    #if FULL_DOUBLE
    typedef double scalar_type;
    #else
    typedef float scalar_type;
    #endif
    cudaGetDevice(&prev_gpu);
    cudaGetDeviceCount(&gpu_count);
    total_threads = gpu_count;

    correction = std::vector<std::vector<double> >(gpu_count, std::vector<double>(gpu_count, 0));

    std::vector<double> runtimes(gpu_count, 0.0);
    for(int i = 0; i < gpu_count ; i++) {
      cudaSetDevice(i);
      double slow_coef = (i==1)? (1) : 1;
      double result = slow_coef * (double) G2G::do_benchmark<scalar_type>();
      runtimes[i] = result;
    }
    std::cout << "Matriz de correccion de tiempos de GPU: " << std::endl;
    for(int i = 0; i < gpu_count ; i++) {
      for(int j = 0; j < gpu_count ; j++) {
        double corr = runtimes[i]/ runtimes[j];
        correction[j][i] = corr;
        std::cout << corr << "  ";
      }
      std::cout << std::endl;
    }
    cudaSetDevice(prev_gpu);
  }

void Partition::compute_work_partition()
{
  int total_threads = 1;
  #if !CPU_KERNELS
  int gpu_count; cudaGetDeviceCount(&gpu_count);
  total_threads = gpu_count;
  #endif
  #ifndef _OPENMP
  total_threads = 1;
  #endif

  std::vector<long long> thread_estimate(total_threads,0);
  std::vector<std::vector<long long> > work_estimate(total_threads);
  work = std::vector<std::vector<int> >(total_threads);
  for(int i = 0; i < cubes.size()+spheres.size(); i++) {
    int thread_assigned = 0;
    thread_assigned = i % total_threads;
    work[thread_assigned].push_back(i);
    size_t size = (i<cubes.size())? cubes[i].size_in_gpu() : spheres[i-cubes.size()].size_in_gpu();
    //long long size = (i<cubes.size())? cubes[i].cost() : spheres[i-cubes.size()].cost();
    thread_estimate[thread_assigned] += size;
    work_estimate[thread_assigned].push_back(size);
  }
  balance_load(thread_estimate, work_estimate);
}

void Partition::regenerate(void)
{
//	cout << "<============ G2G Partition (" << fortran_vars.grid_type << ")============>" << endl;

    // Determina el exponente minimo para cada tipo de atomo.
    // uno por elemento de la tabla periodica.
    vector<double> min_exps(120, numeric_limits<double>::max());
    for (uint i = 0; i < fortran_vars.m; i++)
    {
        uint contractions = fortran_vars.contractions(i);
        uint nuc = fortran_vars.nucleii(i) - 1;
        uint nuc_type = fortran_vars.atom_types(nuc);
        for (uint j = 0; j < contractions; j++)
        {
            min_exps[nuc_type] = min(min_exps[nuc_type], fortran_vars.a_values(i, j));
        }
    }

    // Un exponente y un coeficiente por funcion.
    vector<double> min_exps_func(fortran_vars.m, numeric_limits<double>::max());
    vector<double> min_coeff_func(fortran_vars.m);
    for (uint i = 0; i < fortran_vars.m; i++)
    {
        uint contractions = fortran_vars.contractions(i);
        for (uint j = 0; j < contractions; j++)
        {
            if (fortran_vars.a_values(i, j) < min_exps_func[i])
            {
                min_exps_func[i] = fortran_vars.a_values(i, j);
                min_coeff_func[i] = fortran_vars.c_values(i, j);
            }
        }
    }

    // Encontrando el prisma conteniendo el sistema.
    double3 x0 = make_double3(0,0,0);
    double3 x1 = make_double3(0,0,0);
    for (uint atom = 0; atom < fortran_vars.atoms; atom++)
    {
        double3 atom_position(fortran_vars.atom_positions(atom));
        uint atom_type = fortran_vars.atom_types(atom);
        // TODO: max_radios esta al doble de lo que deberia porque el criterio no esta bien aplicado.
        double max_radius = 2 * sqrt(max_function_exponent / min_exps[atom_type]);
        _DBG(cout << "tipo: " << atom_type << " " << min_exps[atom_type] << " radio: " << max_radius << endl);
        double3 tuple_max_radius = make_double3(max_radius, max_radius, max_radius);
        if (atom == 0)
        {
            x0 = atom_position - tuple_max_radius;
            x1 = atom_position + tuple_max_radius;
        }
        else
        {
            x0.x = min(x0.x, atom_position.x - max_radius);
            x0.y = min(x0.y, atom_position.y - max_radius);
            x0.z = min(x0.z, atom_position.z - max_radius);

            x1.x = max(x1.x, atom_position.x + max_radius);
            x1.y = max(x1.y, atom_position.y + max_radius);
            x1.z = max(x1.z, atom_position.z + max_radius);
        }
    }

    // El prisma tiene vertices (x,y), con x0 el vertice inferior, izquierdo y mas lejano
    // y x1 el vertice superior, derecho y mas cercano.

    // Generamos la particion en cubos.
    uint3 prism_size = ceil_uint3((x1 - x0) / little_cube_size);

    vector<vector<vector<Cube> > >
      prism(prism_size.x, vector<vector<Cube> >(prism_size.y, vector<Cube>(prism_size.z)));

    // Inicializamos las esferas.
    vector<Sphere> sphere_array;
    if (sphere_radius > 0)
    {
        sphere_array.resize(fortran_vars.atoms);
        for (uint atom = 0; atom < fortran_vars.atoms; atom++)
        {
            uint atom_shells = fortran_vars.shells(atom);
            uint included_shells = (uint)ceil(sphere_radius * atom_shells);
            double radius;
            if (included_shells == 0)
            {
                radius = 0;
            }
            else
            {
                double x = cos((M_PI / (atom_shells + 1)) * (atom_shells - included_shells + 1));
                double rm = fortran_vars.rm(atom);
                radius = rm * (1.0 + x) / (1.0 - x);
            }
            _DBG(cout << "esfera incluye " << included_shells << " capas de " << atom_shells << " (radio: " << radius << ")" << endl);
            sphere_array[atom] = Sphere(atom, radius);
        }
    }

    // Precomputamos las distancias entre atomos.
    for (uint i = 0; i < fortran_vars.atoms; i++)
    {
        const double3& atom_i_position(fortran_vars.atom_positions(i));
        double nearest_neighbor_dist = numeric_limits<double>::max();

        double sphere_i_radius = (sphere_radius > 0 ? sphere_array[i].radius : 0);

        for (uint j = 0; j < fortran_vars.atoms; j++)
        {
            const double3& atom_j_position(fortran_vars.atom_positions(j));
            double dist = length(atom_i_position - atom_j_position);
            fortran_vars.atom_atom_dists(i, j) = dist;
            if (i != j)
                nearest_neighbor_dist = min(nearest_neighbor_dist, dist);
        }
        fortran_vars.nearest_neighbor_dists(i) = nearest_neighbor_dist;
    }

    // Computamos los puntos y los asignamos a los cubos y esferas.
    uint puntos_totales = 0;
    uint puntos_finales = 0;
    uint funciones_finales = 0;
    uint costo = 0;

    // Limpiamos las colecciones de las esferas y cubos que tengamos guardadas.
    this->clear();

    // Computa las posiciones de los puntos (y los guarda).
    for (uint atom = 0; atom < fortran_vars.atoms; atom++)
    {
        uint atom_shells = fortran_vars.shells(atom);
        const double3& atom_position(fortran_vars.atom_positions(atom));

        double t0 = M_PI / (atom_shells + 1);
        double rm = fortran_vars.rm(atom);

        puntos_totales += (uint) fortran_vars.grid_size * atom_shells;
        for (uint shell = 0; shell < atom_shells; shell++)
        {
            double t1 = t0 * (shell + 1);
            double x = cos(t1);
            double w = t0 * abs(sin(t1));
            double r1 = rm * (1.0 + x) / (1.0 - x);
            double wrad = w * (r1 * r1) * rm * 2.0 / ((1.0 - x) * (1.0 - x));

            for (uint point = 0; point < (uint)fortran_vars.grid_size; point++)
            {
                double3 rel_point_position = make_double3(fortran_vars.e(point,0), fortran_vars.e(point,1), fortran_vars.e(point,2));
                double3 point_position = atom_position + rel_point_position * r1;
                bool inside_prism = ((x0.x <= point_position.x && point_position.x <= x1.x) &&
                                     (x0.y <= point_position.y && point_position.y <= x1.y) &&
                                     (x0.z <= point_position.z && point_position.z <= x1.z));
                if (inside_prism)
                {
                    double point_weight = wrad * fortran_vars.wang(point); // integration weight
                    Point point_object(atom, shell, point, point_position, point_weight);
                    uint included_shells = (uint)ceil(sphere_radius * atom_shells);

                    // Si esta capa esta muy lejos del nucleo, la modelamos como esfera, sino como cubo.
                    if (shell >= (atom_shells - included_shells))
                    {
                        // Asignamos este punto a la esfera de este atomo.
                        Sphere& sphere = sphere_array[atom];
                        sphere.add_point(point_object);
                    }
                    else
                    {
                        // Insertamos este punto en el cubo correspondiente.
                        uint3 cube_coord = floor_uint3((point_position - x0) / little_cube_size);
                        if (cube_coord.x >= prism_size.x || cube_coord.y >= prism_size.y || cube_coord.z >= prism_size.z)
                            throw std::runtime_error("Se accedio a un cubo invalido");
                        prism[cube_coord.x][cube_coord.y][cube_coord.z].add_point(point_object);
                    }
                }
            }
        }
    }

    // La grilla computada ahora tiene |puntos_totales| puntos, y |fortran_vars.m| funciones.
    uint nco_m = 0;
    uint m_m = 0;

    puntos_finales = 0;
    // Completamos los parametros de los cubos y los agregamos a la particion.
    for (uint i = 0; i < prism_size.x; i++)
    {
        for (uint j = 0; j < prism_size.y; j++)
        {
            for (uint k = 0; k < prism_size.z; k++)
            {
                Cube& cube_ijk = prism[i][j][k];

                double3 cube_coord_abs = x0 + make_uint3(i,j,k) * little_cube_size;

                cube_ijk.assign_significative_functions(cube_coord_abs, min_exps_func, min_coeff_func);
                if (cube_ijk.total_functions_simple() == 0) // Este cubo no tiene funciones.
                    continue;
                if (cube_ijk.number_of_points < min_points_per_cube) // Este cubo no tiene suficientes puntos.
                    continue;

                Cube cube(cube_ijk);
                assert(cube.number_of_points != 0);
                cube.compute_weights();

                if (cube.number_of_points < min_points_per_cube)
                {
                    cout << "not enough points" << endl;
                    continue;
                }
                cubes.push_back(cube);

                // para hacer histogramas
//#ifdef HISTOGRAM
                //cout << "[" << fortran_vars.grid_type << "] cubo: (" << i << "," << j << "," << k << "): " << cube.number_of_points << " puntos; " <<
                     //cube.total_functions() << " funciones, vecinos: " << cube.total_nucleii() << endl;
//#endif

                puntos_finales += cube.number_of_points;
                funciones_finales += cube.number_of_points * cube.total_functions();
                costo += cube.number_of_points * (cube.total_functions() * cube.total_functions());
                nco_m += cube.total_functions() * fortran_vars.nco;
                m_m += cube.total_functions() * cube.total_functions();
            }
        }
    }
    sortBySize<Cube>(cubes);

    // Si esta habilitada la particion en esferas, entonces clasificamos y las agregamos a la particion tambien.
    if (sphere_radius > 0)
    {
        for (uint i = 0; i < fortran_vars.atoms; i++)
        {
            Sphere& sphere_i = sphere_array[i];

            assert(sphere_i.number_of_points != 0);

            sphere_i.assign_significative_functions(min_exps_func, min_coeff_func);
            assert(sphere_i.total_functions_simple() != 0);
            if (sphere_i.number_of_points < min_points_per_cube)
            {
                cout << "not enough points" << endl;
                continue;
            }

            Sphere sphere(sphere_i);
            assert(sphere.number_of_points != 0);
            sphere.compute_weights();
            if (sphere.number_of_points < min_points_per_cube)
            {
                cout << "not enough points" << endl;
                continue;
            }
            assert(sphere.number_of_points != 0);
            spheres.push_back(sphere);

//#ifdef HISTOGRAM
            //cout << "sphere: " << sphere.number_of_points << " puntos, " << sphere.total_functions() <<
                // " funciones | funcion x punto: " << sphere.total_functions() / (double)sphere.number_of_points <<
                // " vecinos: " << sphere.total_nucleii() << endl;
//#endif

            puntos_finales += sphere.number_of_points;
            funciones_finales += sphere.number_of_points * sphere.total_functions();
            costo += sphere.number_of_points * (sphere.total_functions() * sphere.total_functions());
            nco_m += sphere.total_functions() * fortran_vars.nco;
            m_m += sphere.total_functions() * sphere.total_functions();
        }
    }
    //Sorting the spheres in increasing order
    sortBySize<Sphere>(spheres);

    //Initialize the global memory pool for CUDA, with the default safety factor
    //If it is CPU, then this doesn't matter
    GlobalMemoryPool::init(G2G::free_global_memory);
    generate_gpu_profile();
    compute_work_partition();
    //cout << "Grilla final: " << puntos_finales << " puntos (recordar que los de peso 0 se tiran), " << funciones_finales << " funciones" << endl ;
    //cout << "Costo: " << costo << endl;
    //cout << "NCOxM: " << nco_m << " MxM: " << m_m << endl;
    //cout << "Particion final: " << spheres.size() << " esferas y " << cubes.size() << " cubos" << endl;
}
