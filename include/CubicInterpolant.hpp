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

#ifndef CUBICINTERPOLANT_H
#define CUBICINTERPOLANT_H

#include <TriCubicInterpolator.h>
#include <vector>

namespace tbslas {

template<typename real_t>
class CubicInterpolant {
 public:
  CubicInterpolant();
  virtual ~CubicInterpolant();

 public:
  std::vector<real_t> interp(std::vector<real_t>& grid_pnts, int sdim,
                             std::vector<real_t>& grid_vals, int vdim,
                             std::vector<real_t>& query_pnts);

 private:
  //  likely::TriCubicInterpolator tri_cubic_inerpolator;
};

}  // namespace tbslas

#include "CubicInterpolant.ipp"

#endif /* CUBICINTERPOLANT_H */
