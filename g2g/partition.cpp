#include <iostream>
#include <limits>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "common.h"
#include "init.h"
#include "matrix.h"
#include "partition.h"
#include "timer.h"
using namespace std;

namespace G2G {
Partition partition;

ostream& operator<<(ostream& io, const Timers& t) {
#ifdef TIMINGS
  cout << "iteration: " << t.total << endl;
  cout << "rmm: " << t.rmm << " density: " << t.density << " pot: " << t.pot << " forces: " << t.forces << " resto: " << t.resto << " functions: " << t.functions << endl;
#endif
  return io;
}

/********************
 * PointGroup
 ********************/

template<class scalar_type>
void PointGroup<scalar_type>::get_rmm_input(HostMatrix<scalar_type>& rmm_input,
    FortranMatrix<double>& source) const {
  rmm_input.zero();
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);

          rmm_input(ii, jj) = (scalar_type)source.data[big_index];

          rmm_input(jj, ii) = rmm_input(ii, jj);
        }
      }
    }
  }
}

template<class scalar_type>
void PointGroup<scalar_type>::get_rmm_input(HostMatrix<scalar_type>& rmm_input) const {
  get_rmm_input(rmm_input, fortran_vars.rmm_input_ndens1);
}

template<class scalar_type>
void PointGroup<scalar_type>::get_rmm_input(HostMatrix<scalar_type>& rmm_input_a, HostMatrix<scalar_type>& rmm_input_b) const {
  get_rmm_input(rmm_input_a, fortran_vars.rmm_dens_a);
  get_rmm_input(rmm_input_b, fortran_vars.rmm_dens_b);
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_output(const HostMatrix<scalar_type>& rmm_output,
    FortranMatrix<double>& target ) const {
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
          target(big_index) += (double)rmm_output(ii, jj);
        }
      }
    }
  }
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_output(const HostMatrix<scalar_type>& rmm_output) const {
  add_rmm_output(rmm_output, fortran_vars.rmm_output);
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_output(
    const HostMatrix<scalar_type>& rmm_output, HostMatrix<double>& target) const {
  for (uint i = 0, ii = 0; i < total_functions_simple(); i++) {
    uint inc_i = small_function_type(i);

    for (uint k = 0; k < inc_i; k++, ii++) {
      uint big_i = local2global_func[i] + k;
      for (uint j = 0, jj = 0; j < total_functions_simple(); j++) {
        uint inc_j = small_function_type(j);

        for (uint l = 0; l < inc_j; l++, jj++) {
          uint big_j = local2global_func[j] + l;
          if (big_i > big_j) continue;
          uint big_index = (big_i * fortran_vars.m - (big_i * (big_i - 1)) / 2) + (big_j - big_i);
          target(big_index) += (double)rmm_output(ii, jj);
        }
      }
    }
  }

}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_output_a(const HostMatrix<scalar_type>& rmm_output) const {
  add_rmm_output(rmm_output, fortran_vars.rmm_output_a);
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_output_b(const HostMatrix<scalar_type>& rmm_output) const {
  add_rmm_output(rmm_output, fortran_vars.rmm_output_b);
}

template<class scalar_type>
void PointGroup<scalar_type>::add_rmm_open_output(const HostMatrix<scalar_type>& rmm_output_a,
    const HostMatrix<scalar_type>& rmm_output_b) const {
  add_rmm_output(rmm_output_a, fortran_vars.rmm_output_a);
  add_rmm_output(rmm_output_b, fortran_vars.rmm_output_b);
}

template<class scalar_type>
void PointGroup<scalar_type>::compute_nucleii_maps(void)
{
  if (total_functions_simple() != 0) {
    func2global_nuc.resize(total_functions_simple());
    for (uint i = 0; i < total_functions_simple(); i++) {
      func2global_nuc(i) = fortran_vars.nucleii(local2global_func[i]) - 1;
    }

    func2local_nuc.resize(total_functions());
    uint ii = 0;
    for (uint i = 0; i < total_functions_simple(); i++) {
      uint global_atom = func2global_nuc(i);
      uint local_atom = std::distance(local2global_nuc.begin(),
          std::find(local2global_nuc.begin(), local2global_nuc.end(), global_atom));
      uint inc = small_function_type(i);
      for (uint k = 0; k < inc; k++, ii++) func2local_nuc(ii) = local_atom;
    }
  }
}

template<class scalar_type>
void PointGroup<scalar_type>::add_point(const Point& p) {
  points.push_back(p);
  number_of_points++;
}

#define EXP_PREFACTOR 1.01057089636005 // (2 * pow(4, 1/3.0)) / M_PI

template<class scalar_type>
bool PointGroup<scalar_type>::is_significative(FunctionType type, double exponent, double coeff, double d2) {
  switch(type) {
    case FUNCTION_S:
      return (exponent * d2 < max_function_exponent-log(pow((2.*exponent/M_PI),3))/4);
    break;
    default:
    {
      double x = 1;
      double delta;
      double e = 0.1;
      double factor = pow((2.0*exponent/M_PI),3);
      factor = sqrt(factor*4.0*exponent) ;
      double norm = (type == FUNCTION_P ? sqrt(factor) : abs(factor)) ;
      do {
        double div = (type == FUNCTION_P ? log(x) : 2 * log(x));
        double x1 = sqrt((max_function_exponent - log(norm) + div) / exponent);
        delta = abs(x-x1);
        x = x1;
      } while (delta > e);
      return (sqrt(d2) < x);
    }
    break;
  }
}
template<class scalar_type>
bool PointGroup<scalar_type>::operator<(const PointGroup<scalar_type>& T) const{
    int my_cost = number_of_points * total_functions();
    int T_cost = T.number_of_points * T.total_functions();
    return my_cost < T_cost;
}
template<class scalar_type>
size_t PointGroup<scalar_type>::size_in_gpu() const
{
    uint total_cost=0;
    uint single_matrix_cost = COALESCED_DIMENSION(number_of_points) * total_functions();

    total_cost += 2*single_matrix_cost;       //1 scalar_type functions * 2 (matrix and its transposed)
    if (fortran_vars.do_forces || fortran_vars.gga)
      total_cost += (single_matrix_cost*4) * 2; //4 vec_type gradient and its transposed
    if (fortran_vars.gga)
      total_cost+= (single_matrix_cost*8);  //2*4 vec_type hessian
    return total_cost*sizeof(scalar_type);  // size in bytes according to precision
}
template<class scalar_type>
PointGroup<scalar_type>::~PointGroup<scalar_type>() {
#if !CPU_KERNELS
  if(inGlobal) {
    GlobalMemoryPool::dealloc(size_in_gpu());
    function_values.deallocate();
    function_values_transposed.deallocate();
    gradient_values.deallocate();
    gradient_values_transposed.deallocate();
    hessian_values.deallocate();
    hessian_values_transposed.deallocate();
  }
#endif
}

void Partition::solve(Timers& timers, bool compute_rmm,bool lda,bool compute_forces,
                      bool compute_energy, double* fort_energy_ptr, double* fort_forces_ptr, bool OPEN) {
  double cubes_energy = 0, spheres_energy = 0;
  double cubes_energy_i = 0, spheres_energy_i = 0;
  double cubes_energy_c = 0, spheres_energy_c = 0;
  double cubes_energy_c1 = 0, spheres_energy_c1 = 0;
  double cubes_energy_c2 = 0, spheres_energy_c2 = 0;

  int gpu_count;
  cudaGetDeviceCount(&gpu_count);
  if (gpu_count == 0) {
    std::cout << "Error: No GPU found" << std::endl;
    exit(1);
  }
  int total_threads = gpu_count;
  #ifndef _OPENMP
  total_threads = 1;
  gpu_count = 1;
  #else
  omp_set_num_threads(total_threads);
  #endif
  double energy_cubes[total_threads];
  double energy_spheres[total_threads];

  HostMatrix<double> fort_forces_ms[total_threads];

  if (compute_forces) {
      for(int i = 0; i < total_threads; i++) {
          fort_forces_ms[i].resize(fortran_vars.max_atoms, 3);
          fort_forces_ms[i].zero();
      }
  }

  HostMatrix<double> rmm_outputs[total_threads];
  if (compute_rmm) {
      for(int i = 0; i < total_threads; i++) {
          rmm_outputs[i].resize(fortran_vars.rmm_output.width, fortran_vars.rmm_output.height);
          rmm_outputs[i].zero();
      }
  }
#pragma omp parallel shared(energy_spheres, energy_cubes)
  {
    int my_thread = 0;
    #ifdef _OPENMP
    my_thread = omp_get_thread_num();
    #endif
    energy_spheres[my_thread] = 0.0f;
    energy_cubes[my_thread] = 0.0f;
    cudaSetDevice(my_thread % gpu_count);
    Timer t0;
    t0.start_and_sync();
    for(int i = my_thread; i < cubes.size() + spheres.size(); i+= total_threads) {
      if(i < cubes.size()) {
        cubes[i].solve(
            timers, compute_rmm,lda,compute_forces, compute_energy, energy_cubes[my_thread], cubes_energy_i,
            cubes_energy_c, cubes_energy_c1, cubes_energy_c2, fort_forces_ms[my_thread], rmm_outputs[my_thread], OPEN);
      }
      else {
        spheres[i - cubes.size()].solve(
            timers, compute_rmm,lda,compute_forces, compute_energy, energy_spheres[my_thread],
            spheres_energy_i, spheres_energy_c, spheres_energy_c1, spheres_energy_c2, fort_forces_ms[my_thread],
            rmm_outputs[my_thread], OPEN);
      }
    }
    t0.stop_and_sync();
    std::cout << "Workload " << my_thread << " " << t0 << std::endl;
  }

  if (compute_forces) {
      FortranMatrix<double> fort_forces_out(fort_forces_ptr, fortran_vars.atoms, 3, fortran_vars.max_atoms);
      for(int k = 0; k < total_threads; k++) {
          for(int i = 0; i < fortran_vars.atoms; i++) {
              for(int j = 0; j < 3; j++) {
                  fort_forces_out(i,j) += fort_forces_ms[k](i,j);
              }
          }
      }
  }

  if(compute_rmm) {
    for(int k = 0; k < total_threads; k++) {
      for(int i = 0; i < rmm_outputs[k].width; i++) {
        for(int j = 0; j < rmm_outputs[k].height; j++) {
          fortran_vars.rmm_output(i,j) += rmm_outputs[k](i,j);
        }
      }
    }
  }

  for(int i = 0; i< total_threads; i++) {
    cubes_energy += energy_cubes[i];
    spheres_energy += energy_spheres[i];
  }

  *fort_energy_ptr = cubes_energy + spheres_energy;
  if(*fort_energy_ptr != *fort_energy_ptr) {
      std::cout << "I see dead peaple " << std::endl;
#ifndef CPU_KERNELS
      cudaDeviceReset();
#endif
      exit(1);
   }

}



/**********************
 * Sphere
 **********************/
Sphere::Sphere(void) : atom(0), radius(0) { }
Sphere::Sphere(uint _atom, double _radius) : atom(_atom), radius(_radius) { }

/**********************
 * Cube
 **********************/

template class PointGroup<double>;
template class PointGroup<float>;
}
