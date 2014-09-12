// Copyright (C) 2014 by Arash Bakhtiari

// *************************************************************************
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

#ifndef INCLUDE_VEC_FIELD_H_
#define INCLUDE_VEC_FIELD_H_

#include <cstdlib>
#include <vector>

namespace tbslas {

template<typename real_t, int sdim = 3, int vdim = 3>
class VecField {
 public:
  VecField();
  VecField(std::vector<real_t> field_points,
           std::vector<real_t> field_values);
  virtual ~VecField();

 public:
  void
  init(std::vector<real_t> field_points,
       std::vector<real_t> field_values);

  template<typename InterpPolicy>
  void
  interp(const std::vector<real_t>& query_points,
         const InterpPolicy& interpolant,
         std::vector<real_t>& query_values) const;

  void write2file(const char* file_name);

 private:
  std::vector<real_t> field_points_;  // vector field's space subset
  std::vector<real_t> field_values_;  // values of the vector field's
};  // class VecField

}  // namespace tbslas

#include "vec_field.inc"

#endif  // INCLUDE_VEC_FIELD_H_
