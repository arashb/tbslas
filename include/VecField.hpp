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
 * @file   VecField.hpp
 * @author Arash Bakhtiari <bakhtiar@in.tum.de>
 * @date   Thu Sep  4 11:00:48 2014
 *
 * @brief  This class represent the vector field concept.
 *
 *
 */

#ifndef VECFIELD_H
#define VECFIELD_H

#include <cstdlib>
#include <vector>

namespace tbslas {

template<typename real_t, int sdim = 3, int vdim = 3>
class VecField {
 public:
  VecField();
  VecField(std::vector<real_t> pnts,
           std::vector<real_t> vls);
  virtual ~VecField();

 public:
  void
  init(std::vector<real_t> pnts,
       std::vector<real_t> vls);

  template<typename C>
  void
  interp(std::vector<real_t>& qry_pnts,
         C& interpolant,
         std::vector<real_t>& qry_vls);

  void write2file(const char* fname);

 private:
  std::vector<real_t> _pnts;   // vector field's space subset
  std::vector<real_t> _vls;    // values of the vector field's
};  // class VecField

}  // namespace tbslas

#include "VecField.ipp"

#endif /* VECFIELD_H */
