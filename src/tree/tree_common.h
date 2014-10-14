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

#ifndef SRC_TREE_TREE_COMMON_H_
#define SRC_TREE_TREE_COMMON_H_

#include <pvfmm_common.hpp>
#include <mpi_tree.hpp>
#include <cheb_node.hpp>

namespace tbslas {

// only compiles with c++11 flag
template <typename real_t>
using Node_t = pvfmm::Cheb_Node<real_t>;

template <typename real_t>
using Tree_t = pvfmm::MPI_Tree< Node_t<real_t> >;

}  // namespace tbslas

#endif  // SRC_TREE_TREE_COMMON_H_
