// *************************************************************************
// Copyright (C) 2014 by Arash Bakhtiari
// You may not use this file except in compliance with the License.
// You obtain a copy of the License in the LICENSE file.

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// *************************************************************************

#ifndef SRC_UTILS_COMMON_H_
#define SRC_UTILS_COMMON_H_

// SYSTEM HEADERS
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// EXTERNAL HEADERS
#include <mpi.h>

// PROJECT HEADERS
#include "utils/sim_config.h"
#include "utils/singleton.h"
#include "utils/types.h"

namespace tbslas {

// ************************************************************
// GLOBAL DECLARATION
// ************************************************************
// simulation configuration singleton
typedef Singleton<SimConfig> SimConfigSingleton;

// ************************************************************
// COMMON FUNCTIONS
// ************************************************************

template <typename T>
std::string ToString(T val) {
  std::stringstream stream;
  stream << val;
  return stream.str();
}

template <typename PointerType>
void swap_pointers(PointerType **ta, PointerType **tb) {
  PointerType *tmp = *ta;
  *ta = *tb;
  *tb = tmp;
}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string get_current_datetime() {
  time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
  return buf;
}

const std::string get_result_dir() {
  char res_dir_buffer[200];
  snprintf(res_dir_buffer, sizeof(res_dir_buffer), "%s/",
           getenv("TBSLAS_RESULT_DIR"));

  return res_dir_buffer;
}

std::string GetVTKFileName(int timestep, std::string variable_id) {
  tbslas::SimConfig *sim_config = tbslas::SimConfigSingleton::Instance();
  char buffer[500];
  snprintf(buffer, sizeof(buffer), sim_config->vtk_filename_format.c_str(),
           tbslas::get_result_dir().c_str(),
           // sim_config->vtk_filename_prefix.c_str(),
           variable_id.c_str(), timestep);
  std::ostringstream stringStream;
  stringStream << buffer;
  return stringStream.str();
}

typedef void (*dummy_fn_ptr)(const double *coord, int n, double *out);

template <typename real_t>
void dummy_fn(const real_t *points_pos, int num_points, real_t *values) {}

inline bool is_little_endian() {
  int n = 1;
  return (*reinterpret_cast<char *>(&n) == 1);
}

template <class Real_t>
std::vector<Real_t> point_distrib(DistribType dist_type, size_t N,
                                  MPI_Comm comm) {
  int np, myrank;
  MPI_Comm_size(comm, &np);
  MPI_Comm_rank(comm, &myrank);
  std::vector<Real_t> coord;
  srand48(myrank + 1);

  switch (dist_type) {
    case UnifGrid: {
      size_t NN = (size_t)round(pow((double)N, 1.0 / 3.0));
      size_t N_total = NN * NN * NN;
      size_t start = myrank * N_total / np;
      size_t end = (myrank + 1) * N_total / np;
      for (size_t i = start; i < end; i++) {
        coord.push_back(((Real_t)((i / 1) % NN) + 0.5) / NN);
        coord.push_back(((Real_t)((i / NN) % NN) + 0.5) / NN);
        coord.push_back(((Real_t)((i / (NN * NN)) % NN) + 0.5) / NN);
      }
    } break;
    case RandUnif: {
      size_t N_total = N;
      size_t start = myrank * N_total / np;
      size_t end = (myrank + 1) * N_total / np;
      size_t N_local = end - start;
      coord.resize(N_local * 3);

      for (size_t i = 0; i < N_local * 3; i++) coord[i] = ((Real_t)drand48());
    } break;
    case RandGaus: {
      Real_t e = 2.7182818284590452;
      Real_t log_e = log(e);
      size_t N_total = N;
      size_t start = myrank * N_total / np;
      size_t end = (myrank + 1) * N_total / np;

      for (size_t i = start * 3; i < end * 3; i++) {
        Real_t y = -1;
        while (y <= 0 || y >= 1) {
          Real_t r1 =
              sqrt(-2 * log(drand48()) / log_e) * cos(2 * M_PI * drand48());
          Real_t r2 = pow(0.5, i * 10 / N_total);
          y = 0.5 + r1 * r2;
        }
        coord.push_back(y);
      }
    } break;
    case RandElps: {
      size_t N_total = N;
      size_t start = myrank * N_total / np;
      size_t end = (myrank + 1) * N_total / np;
      size_t N_local = end - start;
      coord.resize(N_local * 3);

      const Real_t r = 0.25;
      const Real_t center[3] = {0.5, 0.5, 0.5};
      for (size_t i = 0; i < N_local; i++) {
        Real_t *y = &coord[i * 3];
        Real_t phi = 2 * M_PI * drand48();
        Real_t theta = M_PI * drand48();
        y[0] = center[0] + 0.25 * r * sin(theta) * cos(phi);
        y[1] = center[1] + 0.25 * r * sin(theta) * sin(phi);
        y[2] = center[2] + r * cos(theta);
      }
    } break;
    case RandSphr: {
      size_t N_total = N;
      size_t start = myrank * N_total / np;
      size_t end = (myrank + 1) * N_total / np;
      size_t N_local = end - start;
      coord.resize(N_local * 3);

      const Real_t center[3] = {0.5, 0.5, 0.5};
      for (size_t i = 0; i < N_local; i++) {
        Real_t *y = &coord[i * 3];
        Real_t r = 1;
        while (r > 0.5 || r == 0) {
          y[0] = drand48();
          y[1] = drand48();
          y[2] = drand48();
          r = sqrt((y[0] - center[0]) * (y[0] - center[0]) +
                   (y[1] - center[1]) * (y[1] - center[1]) +
                   (y[2] - center[2]) * (y[2] - center[2]));
          y[0] = center[0] + 0.1 * (y[0] - center[0]) / r;
          y[1] = center[1] + 0.1 * (y[1] - center[1]) / r;
          y[2] = center[2] + 0.1 * (y[2] - center[2]) / r;
        }
      }
    } break;
    default:
      break;
  }
  return coord;
}

template <typename real_t, int sdim>
void get_reg_grid_points(const size_t N, real_t *out) {
  size_t tN = std::pow(static_cast<double>(N), sdim);
  real_t dx = 1.0 / (N - 1);
  // regular grid points positions
  for (size_t i = 0; i < tN; i++)
    for (int idim = 0; idim < sdim; idim++) {
      size_t denom_shift = (size_t)std::pow(static_cast<double>(N), idim);
      out[i * sdim + idim] = (((real_t)((i / denom_shift) % N)) * dx);
    }
}

template <typename real_t, int dim>
std::vector<real_t> generate_reg_grid_points(const size_t N) {
  size_t tN = std::pow(static_cast<double>(N), dim);
  real_t dx = 1.0 / (N - 1);
  // regular grid points positions
  std::vector<real_t> reg_grid_points_pos(tN * dim);
  for (int idim = 0; idim < dim; idim++) {
    // idim's cooridnate
    int index_shift = idim * tN;
    size_t denom_shift = (size_t)std::pow(N, idim);
    for (size_t i = 0; i < tN; i++)
      reg_grid_points_pos[i + index_shift] =
          (((real_t)((i / denom_shift) % N)) * dx);
  }
  return reg_grid_points_pos;
}

template <typename real_t>
inline int find_grid_index_1d(const std::vector<real_t> &grid,
                              const real_t query) {
  int num_el = grid.size();
  assert(num_el > 2);
  // assume grid is sorted
  real_t init = grid[0];
  real_t final = grid[num_el - 1];
  // assume grid points are equidistant
  real_t spacing = grid[1] - grid[0];
  assert(query >= init);
  assert(query <= final);
  int index = (int)std::floor((query - init) / spacing);
  assert(index >= 0);
  assert(index < num_el);
  return index;
}

/* cheb_cfl = dt/(2^-L*q^-2) */
template <typename real_t>
inline real_t compute_cheb_CFL(real_t dt, int max_depth, int cheb_degree) {
  tbslas::SimConfig *sim_config = tbslas::SimConfigSingleton::Instance();
  int myrank;
  MPI_Comm_rank(sim_config->comm, &myrank);
  real_t cfl = dt * std::pow(2.0, max_depth) * cheb_degree * cheb_degree;
  if (!myrank) std::cout << "# CFL: " << cfl << std::endl;
  return cfl;
}

}  // namespace tbslas

#endif  // SRC_UTILS_COMMON_H_
