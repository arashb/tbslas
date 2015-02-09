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
#include <tree/semilag_tree.h>
#include <tree/utils_tree.h>

typedef pvfmm::Cheb_Node<double> Node_t;
typedef pvfmm::MPI_Tree<Node_t> Tree_t;

typedef tbslas::MetaData<std::string,
                         std::string,
                         std::string> MetaData_t;
double tcurr = 0;

int main (int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm=MPI_COMM_WORLD;
  int myrank;
  MPI_Comm_rank(comm, &myrank);

  parse_command_line_options(argc, argv);

    int   test = strtoul(commandline_option(argc, argv, "-test",     "1", false,
                                          "-test <int> = (1)    : 1) Gaussian profile 2) Zalesak disk"),NULL,10);

  {
    tbslas::SimConfig* sim_config       = tbslas::SimConfigSingleton::Instance();
    pvfmm::Profile::Enable(sim_config->profile);
    // =========================================================================
    // SIMULATION PARAMETERS
    // =========================================================================
    sim_config->vtk_filename_prefix     = "advection";
    sim_config->vtk_filename_variable   = "conc";
    // =========================================================================
    // PRINT METADATA
    // =========================================================================
    if (!myrank) {
      MetaData_t::Print();
    }
    // =========================================================================
    // INIT THE VELOCITY TREE
    // =========================================================================    
    Tree_t tvel_curr(comm);
    tbslas::ConstructTree<Tree_t>(sim_config->tree_num_point_sources,
                                  sim_config->tree_num_points_per_octanct,
                                  sim_config->tree_chebyshev_order,
                                  sim_config->tree_max_depth,
                                  sim_config->tree_adap,
                                  sim_config->tree_tolerance,
                                  comm,
                                  tbslas::get_vorticity_field<double,3>,
                                  3,
                                  tvel_curr);
    char out_name_buffer[300];
    snprintf(out_name_buffer,
         sizeof(out_name_buffer),
         sim_config->vtk_filename_format.c_str(),
         tbslas::get_result_dir().c_str(),
         sim_config->vtk_filename_prefix.c_str(),
         "vel",
         0);
    tvel_curr.Write2File(out_name_buffer, sim_config->vtk_order);
    // =========================================================================
    // INIT THE CONCENTRATION TREE
    // =========================================================================
    tcurr = 0;
    Tree_t tconc_curr(comm);
    tbslas::ConstructTree<Tree_t>(sim_config->tree_num_point_sources,
                                  sim_config->tree_num_points_per_octanct,
                                  sim_config->tree_chebyshev_order,
                                  sim_config->tree_max_depth,
                                  sim_config->tree_adap,
                                  sim_config->tree_tolerance,
                                  comm,
                                  (test==1)?get_gaussian_field_cylinder_atT<double,3>:
                                  get_slotted_cylinder_atT<double,3>,
                                  1,
                                  tconc_curr);
    double in_al2,in_rl2,in_ali,in_rli;
    CheckChebOutput<Tree_t>(&tconc_curr,
                            (test==1)?get_gaussian_field_cylinder_atT<double,3>:
                            get_slotted_cylinder_atT<double,3>,
                            1,
                            in_al2,in_rl2,in_ali,in_rli,
                            std::string("Input"));
    // =========================================================================
    // RUN
    // =========================================================================
    tbslas::RunSemilagSimulationInSitu(&tvel_curr,
                                       &tconc_curr,
                                       sim_config->total_num_timestep,
                                       sim_config->dt,
                                       sim_config->num_rk_step,
                                       true,
                                       true);
    // =========================================================================
    // COMPUTE ERROR
    // =========================================================================
    tcurr = sim_config->total_num_timestep*sim_config->dt;
    double al2,rl2,ali,rli;
    CheckChebOutput<Tree_t>(&tconc_curr,
                            (test==1)?get_gaussian_field_cylinder_atT<double,3>:
                            get_slotted_cylinder_atT<double,3>,
                            1,
                            al2,rl2,ali,rli,
                            std::string("Output"));
    int num_leaves = tbslas::CountNumLeafNodes(tconc_curr);
    // =========================================================================
    // REPORT RESULTS
    // =========================================================================
  typedef tbslas::Reporter<double> Rep;
  if(!myrank) {
    Rep::AddData("TOL", sim_config->tree_tolerance);
    Rep::AddData("MaxDEPTH", sim_config->tree_max_depth);
    Rep::AddData("NOCT", num_leaves);

    Rep::AddData("DT", sim_config->dt);
    Rep::AddData("TN", sim_config->total_num_timestep);

    Rep::AddData("InAL2", in_al2);
    Rep::AddData("OutAL2", al2);

    Rep::AddData("InRL2", in_rl2);
    Rep::AddData("OutRL2", rl2);

    Rep::AddData("InALINF", in_ali);
    Rep::AddData("OutALINF", ali);

    Rep::AddData("InRLINF", in_rli);
    Rep::AddData("OutRLINF", rli);
    Rep::Report();
  }
    //Output Profiling results.
    pvfmm::Profile::print(&comm);
  }

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
