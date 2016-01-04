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

#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <string>

#include <pvfmm_common.hpp>
#include <mpi_tree.hpp>
#include <cheb_node.hpp>
#include <utils.hpp>
#include <vector.hpp>
#include <cheb_utils.hpp>
#include <profile.hpp>

#include <utils/common.h>
#include <utils/metadata.h>
#include <utils/reporter.h>
#include <utils/fields.h>

#include <tree/semilag_tree.h>
#include <tree/utils_tree.h>

typedef pvfmm::Cheb_Node<double> Node_t;
typedef pvfmm::MPI_Tree<Node_t> Tree_t;

typedef tbslas::MetaData<std::string,
                         std::string,
                         std::string> MetaData_t;
double tcurr = 0;

template<typename real_t, int sdim>
void
get_gaussian_field_cylinder_atT(const real_t* points_pos,
                                int num_points,
                                real_t* out) {
  real_t xc      = 0.6;
  real_t yc      = 0.5;
  real_t r = sqrt((xc-0.5)*(xc-0.5) + (yc-0.5)*(yc-0.5));
  xc = 0.5+r*cos(tcurr);
  yc = 0.5+r*sin(tcurr);
  const real_t theta   = 0.0;
  const real_t sigma_x = 0.06;
  const real_t sigma_y = 0.06;
  const real_t A       = 1.0;
  tbslas::get_gaussian_field_cylinder<real_t, sdim>(points_pos,
                                                    num_points,
                                                    out,
                                                    xc,
                                                    yc,
                                                    theta,
                                                    sigma_x,
                                                    sigma_y,
                                                    A);
}

int main (int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm=MPI_COMM_WORLD;

  // commandline_option_end(argc, argv);
  parse_command_line_options(argc, argv);

  {
    int myrank;
    MPI_Comm_rank(comm, &myrank);
    // =========================================================================
    // SIMULATION PARAMETERS
    // =========================================================================
    tbslas::SimConfig* sim_config       = tbslas::SimConfigSingleton::Instance();
    sim_config->vtk_filename_variable   = "conc";
    // sim_config->use_cubic = true;
    // sim_config->cubic_use_analytical = true;
    pvfmm::Profile::Enable(sim_config->profile);
    // =========================================================================
    // PRINT METADATA
    // =========================================================================
    if (!myrank) {
      MetaData_t::Print();
    }
    // =========================================================================
    // CONSTRUCT TREE
    // =========================================================================
    tcurr = 0;
    Tree_t tree(comm);
    tbslas::ConstructTree<Tree_t>(sim_config->tree_num_point_sources,
                                  sim_config->tree_num_points_per_octanct,
                                  sim_config->tree_chebyshev_order,
                                  sim_config->tree_max_depth,
                                  sim_config->tree_adap,
                                  sim_config->tree_tolerance,
                                  comm,
                                  get_gaussian_field_cylinder_atT<double,3>,
                                  1,
                                  tree);
    // =========================================================================
    // COMPUTE ERROR
    // =========================================================================
    double rli, rl2;
    double ali, al2;
    CheckChebOutput<Tree_t>(&tree,
                            get_gaussian_field_cylinder_atT<double,3>,
                            1,
                            al2,rl2,ali,rli,
                            std::string("Output"));

    std::vector<double> grid_points;
    std::vector<double> node_pos = pvfmm::cheb_nodes<double>(sim_config->tree_chebyshev_order, 3);
    tbslas::CollectTreeGridPoints(tree,
                                  node_pos,
                                  grid_points);

    // tbslas::ComputeTreeError(tree,
    //                          get_gaussian_field_cylinder_atT<double,3>,
    //                          grid_points,
    //                          ali,
    //                          al2);
    int num_leaves = tbslas::CountNumLeafNodes(tree);
    // =========================================================================
    // REPORT RESULTS
    // =========================================================================
  typedef tbslas::Reporter<double> Rep;
  if(!myrank) {
    Rep::AddData("TOL", sim_config->tree_tolerance);
    Rep::AddData("NOCT", num_leaves, tbslas::REP_INT);

    Rep::AddData("CUBIC", sim_config->use_cubic, tbslas::REP_INT);
    Rep::AddData("CUF", sim_config->cubic_upsampling_factor, tbslas::REP_INT);
    // Rep::AddData("ANAL", sim_config->cubic_use_analytical);

    Rep::AddData("OutAL2", al2);
    Rep::AddData("OutALINF", ali);

    Rep::Report();
  }
    //Output Profiling results.
    pvfmm::Profile::print(&comm);
  }

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
