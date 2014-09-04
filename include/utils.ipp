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

#include <iostream>
#include <vector>

namespace semilag {

  template <typename real_t, int dim>
  std::vector<real_t> gen_reg_grid_points(size_t N) {
    std::vector<real_t> rgp_pos;    // regular grid points positions
    size_t total = N*N*N;
    real_t dx = 1.0/(N-1);
    real_t xcor;
    real_t ycor;
    real_t zcor;
  for (size_t i = 0; i < total; i++) {
       xcor = ((real_t)((i / 1)    %N))*dx;
       ycor = ((real_t)((i / N)    %N))*dx;
       zcor = ((real_t)((i /(N*N)) %N))*dx;
      // std::cout << xcor << " " << ycor << " " << zcor << std::endl;
      rgp_pos.push_back(xcor);  // x cooridnate
      rgp_pos.push_back(ycor);  // y coordinate
      rgp_pos.push_back(zcor);  // z coordinate
    }
    return rgp_pos;
  }

}  // namespace semilag
