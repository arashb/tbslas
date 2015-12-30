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

void (*fn_vel)(const double* , int , double*)=NULL;

int main (int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm=MPI_COMM_WORLD;
  int np;
  MPI_Comm_size(comm, &np);
  int myrank;
  MPI_Comm_rank(comm, &myrank);

  parse_command_line_options(argc, argv);

  int   test = strtoul(commandline_option(argc, argv, "-test",     "1", false,
                                          "-test <int> = (1)    : 1) Gaussian profile 2) Zalesak disk"),NULL,10);
  int   merge = strtoul(commandline_option(argc, argv, "-merge",     "1", false,
                                          "-merge <int> = (1)    : 1) no merge 2) complete merge 3) Semi-Merge"),NULL,10);

  {
    tbslas::SimConfig* sim_config       = tbslas::SimConfigSingleton::Instance();
    pvfmm::Profile::Enable(sim_config->profile);

    // =========================================================================
    // PRINT METADATA
    // =========================================================================
    if (!myrank) {
      MetaData_t::Print();
    }

    // =========================================================================
    // TEST CASE
    // =========================================================================
    int max_depth_vel=0;
    int max_depth_con=0;
    pvfmm::BoundaryType bc;
    switch(test) {
      case 1:
        fn_vel = tbslas::get_vorticity_field<double,3>;
        break;
    }

    // =========================================================================
    // SIMULATION PARAMETERS
    // =========================================================================
    sim_config->vtk_filename_prefix     = "advection";
    sim_config->vtk_filename_variable   = "conc";
    sim_config->bc = bc;

    // =========================================================================
    // INIT THE VELOCITY TREE
    // =========================================================================
    Tree_t tvel(comm);
    tbslas::ConstructTree<Tree_t>(sim_config->tree_num_point_sources,
                                  sim_config->tree_num_points_per_octanct,
                                  sim_config->tree_chebyshev_order,
                                  max_depth_vel?max_depth_vel:sim_config->tree_max_depth,
                                  sim_config->tree_adap,
                                  sim_config->tree_tolerance,
                                  comm,
                                  fn_vel,
                                  3,
                                  tvel);

    // =========================================================================
    // RUN
    // =========================================================================
    double in_al2,in_rl2,in_ali,in_rli;
    CheckChebOutput<Tree_t>(&tvel,
                            fn_vel,
                            3,
                            in_al2,in_rl2,in_ali,in_rli,
                            std::string("Input"));

    // double xinit = 0.5+1.0/pow(2.0,8.0)/sim_config->tree_chebyshev_order;
    double xinit = 0.99;
    double yinit = 0.5;
    double zinit = 0.5;

    int num_points = 1;
    std::vector<double> init_points_pos(num_points*COORD_DIM);
    std::vector<double> dprt_points_pos(num_points*COORD_DIM);

    init_points_pos[0] = xinit;
    init_points_pos[1] = yinit;
    init_points_pos[2] = zinit;

    tbslas::NodeFieldFunctor<double,Tree_t> vel_evaluator(&tvel);
    tbslas::ComputeTrajRK2(vel_evaluator,
                           init_points_pos,
                           tcurr,
                           tcurr - sim_config->total_num_timestep*sim_config->dt,
                           sim_config->total_num_timestep,
                           dprt_points_pos);

    // =========================================================================
    // COMPUTE ERROR
    // =========================================================================
    double err = 0;
    double err_sum = 0;
    double err_max = 0;
    double l2_err_sum = 0;
    for (int i = 0; i < COORD_DIM; i++) {
      err = init_points_pos[i] - dprt_points_pos[i];
      l2_err_sum += err*err;
      if (abs(err) > err_max)
        err_max = abs(err);
    }
    double l2_err = 0;
    l2_err = sqrt(l2_err_sum);

    // =========================================================================
    // REPORT RESULTS
    // =========================================================================
    int tvel_max_depth=0;
    tbslas::GetTreeMaxDepth(tvel, tvel_max_depth);

    typedef tbslas::Reporter<double> Rep;
    if(!myrank) {
      Rep::AddData("NP", np, tbslas::REP_INT);
      Rep::AddData("OMP",  sim_config->num_omp_threads, tbslas::REP_INT);

      Rep::AddData("TOL", sim_config->tree_tolerance);
      Rep::AddData("Q", sim_config->tree_chebyshev_order, tbslas::REP_INT);

      Rep::AddData("MaxD", sim_config->tree_max_depth, tbslas::REP_INT);
      Rep::AddData("VMaxD", tvel_max_depth, tbslas::REP_INT);

      Rep::AddData("CBC", sim_config->use_cubic?1:0, tbslas::REP_INT);
      Rep::AddData("CUF", sim_config->cubic_upsampling_factor, tbslas::REP_INT);

      Rep::AddData("DT", sim_config->dt);
      Rep::AddData("TN", sim_config->total_num_timestep, tbslas::REP_INT);

      Rep::AddData("InAL2", in_al2);
      Rep::AddData("InRL2", in_rl2);
      Rep::AddData("InALINF", in_ali);
      Rep::AddData("InRLINF", in_rli);

      Rep::AddData("OutL2", l2_err);
      Rep::AddData("OutLINF", err_max);

      Rep::Report();
    }

    //Output Profiling results.
    // pvfmm::Profile::print(&comm);
  }

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
