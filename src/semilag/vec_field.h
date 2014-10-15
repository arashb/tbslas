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

/**
 * @file   vec_field.h
 * @author Arash Bakhtiari <bakhtiar@in.tum.de>
 * @date   Thu Sep  4 11:00:48 2014
 *
 * @brief  This class represent the vector field concept.
 *
 *
 */

#ifndef SRC_SEMILAG_VEC_FIELD_H_
#define SRC_SEMILAG_VEC_FIELD_H_

#include <cstdlib>
#include <vector>

namespace tbslas {

template<typename real_t, int sdim = 3, int vdim = 3>
class VecField {
 public:
  VecField();
  VecField(const std::vector<real_t>& field_points,
           const std::vector<real_t>& field_values,
           const real_t time_init = 0);

  virtual ~VecField();

 public:
  void
  init(const std::vector<real_t>& field_points,
       const std::vector<real_t>& field_values,
       const real_t time_init = 0);


  void
  push_back_values(const std::vector<real_t>& field_values,
                   real_t time = 0);

  template<typename InterpPolicy>
  void
  interp(const std::vector<real_t>& query_points,
         const InterpPolicy& interpolant,
         std::vector<real_t>& query_values,
         const int timestep = 0) const;

  template<typename InterpPolicy>
  void
  interp(const std::vector<real_t>& query_points,
         const real_t query_time,
         const InterpPolicy& interpolant,
         std::vector<real_t>& query_values
         ) const;

  void write2file(const char* file_name, int timestep = 0, int xyz_mem_layout = 1);

  void save(const char* file_name, int xyz_mem_layout = 1);

  const std::vector<real_t>&
  get_points() const {return field_points_;}

  size_t
  get_num_points() const {return num_field_points_;}

 private:
  std::vector<real_t> field_points_;  // vector field's space subset
  size_t num_field_points_;           // number of field points
  // values of the vector field's
  std::vector< std::vector<real_t> > field_values_;
  std::vector<real_t> time_;
};  // class VecField

}  // namespace tbslas

#include "vec_field.inc"

#endif  // SRC_SEMILAG_VEC_FIELD_H_
