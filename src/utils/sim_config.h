// *************************************************************************
// Copyright (C) 2015 by Arash Bakhtiari
// You may not use this file except in compliance with the License.
// You obtain a copy of the License in the LICENSE file.

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// *************************************************************************

#ifndef SRC_UTILS_SIM_CONFIG_H_
#define SRC_UTILS_SIM_CONFIG_H_

#include <mpi.h>
#include <string>

namespace tbslas {

class SimConfig {

 public:
  SimConfig():
      total_num_timestep(1),
      dt(0.1),
      num_rk_step(1),
      vtk_order(14),
      vtk_filename_format("%s/%s-VAR_%s-TS_%04d-RNK"),
      vtk_filename_prefix("output"),
      vtk_filename_variable("vals"),
      use_cubic(false),
      cubic_upsampling_factor(4),
      cubic_use_analytical(false),
      comm(MPI_COMM_WORLD),
      num_omp_threads(1),
      tree_num_point_sources(1),
      tree_num_points_per_octanct(1),
      tree_chebyshev_order(14),
      tree_max_depth(15),
      tree_adap(false),
      profile(false) {
  }
  // *************************************************************************
  // time stepping
  // *************************************************************************
  int total_num_timestep;
  double dt;
  // *************************************************************************
  // semi-lagrangian solver
  // *************************************************************************
  int num_rk_step;
  // *************************************************************************
  // output
  // *************************************************************************
  int vtk_order;
  std::string vtk_filename_format;
  std::string vtk_filename_prefix;
  std::string vtk_filename_variable;
  // *************************************************************************
  // cubic interpolation
  // *************************************************************************
  bool use_cubic;
  int cubic_upsampling_factor;
  bool cubic_use_analytical;
  // *************************************************************************
  // parallelization
  // *************************************************************************
  MPI_Comm comm;
  int num_omp_threads;
  // *************************************************************************
  // chebyshev tree
  // *************************************************************************
  size_t tree_num_point_sources;
  size_t tree_num_points_per_octanct;
  int tree_chebyshev_order;
  int tree_max_depth;
  double tree_tolerance;
  bool tree_adap;
  // *************************************************************************
  // MISC
  // *************************************************************************
  bool profile;
};

}  // namespace tbslas

#endif  // SRC_UTILS_SIM_CONFIG_H_
