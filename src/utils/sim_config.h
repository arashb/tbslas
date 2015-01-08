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
      cubic(false),
      upsampling_factor(4) {
  }

  // time stepping
  int total_num_timestep;
  double dt;

  // semi-lagrangian solver
  int num_rk_step;

  // output
  int vtk_order;
  std::string vtk_filename_format;
  std::string vtk_filename_prefix;
  std::string vtk_filename_variable;

  // cubic interpolation
  bool cubic;
  int upsampling_factor;
};

}  // namespace tbslas

#endif  // SRC_UTILS_SIM_CONFIG_H_
