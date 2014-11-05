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
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <string>

// EXTERNAL HEADERS
#include <mpi.h>

namespace tbslas {

template<typename PointerType>
void swap_pointers(PointerType** ta,
                   PointerType** tb) {
  PointerType* tmp = *ta;
  *ta = *tb;
  *tb = tmp;
}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string get_current_datetime() {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
  return buf;
}

const std::string get_result_dir() {
  char res_dir_buffer[200];
  snprintf(res_dir_buffer, sizeof(res_dir_buffer),
           "%s/", getenv("TBSLAS_RESULT_DIR"));

  return res_dir_buffer;
}

template<typename real_t>
void dummy_fn(const real_t* points_pos,
              int num_points,
              real_t* values) {
}

enum DistribType{
  UnifGrid,
  RandUnif,
  RandGaus,
  RandElps,
  RandSphr
};


template <class Real_t>
std::vector<Real_t> point_distrib(DistribType, size_t N, MPI_Comm comm);

template <typename real_t, int dim>
void
get_reg_grid_points(const size_t N, real_t* out);

template<typename real_t, int sdim>
void
get_gaussian_field(const real_t* pnts_pos,
                   int num_points,
                   real_t* out);

template<typename real_t, int sdim>
void
get_gaussian_field_yext(const real_t* points_pos,
            int num_points,
            real_t* out);

template<typename real_t, int sdim>
void
get_gaussian_field_1d(const real_t* points_pos,
                      int num_points,
                      real_t* out);

template<typename real_t, int sdim>
void
get_gaussian_field_3d(const real_t* points_pos,
                      int num_points,
                      real_t* out);

template<typename real_t, int sdim>
void
get_vorticity_field(const real_t* points_pos,
                    int num_points,
                    real_t* out);

template<typename real_t, int sdim>
void
get_vel_field_hom(const real_t* points_pos,
                  int num_points,
                  real_t* points_values);

template <typename real_t, int dim>
std::vector<real_t>
generate_reg_grid_points(const size_t N);

template<typename real_t>
std::vector<real_t>
generate_vorticity_field(const std::vector<real_t>& pnts_pos,
                         const real_t time  = 0,
                         const real_t omega = 1);

template<typename real_t>
std::vector<real_t>
generate_gaussian_field(const std::vector<real_t>& pnts_pos,
                        const real_t xcenter = 0.6,
                        const real_t ycenter = 0.6,
                        const real_t theta   = 0.0,
                        const real_t sigma_x = 0.06,
                        const real_t sigma_y = 0.06);

template<typename real_t>
inline int
find_grid_index_1d(const std::vector<real_t>& grid,
                   const real_t query);

inline bool
is_little_endian() {
  int n = 1;
  return (*reinterpret_cast<char*>(&n) == 1);
}

}  // namespace tbslas

#include "common.inc"

#endif  // SRC_UTILS_COMMON_H_
