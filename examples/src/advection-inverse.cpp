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

#include <tree/tree_semilag.h>
#include <tree/tree_utils.h>

typedef pvfmm::Cheb_Node<double> Node_t;
typedef pvfmm::MPI_Tree<Node_t> Tree_t;

typedef tbslas::MetaData<std::string,
                         std::string,
                         std::string> MetaData_t;
double tcurr = 0;

void (*fn_vel)(const double* , int , double*)=NULL;
void (*fn_con)(const double* , int , double*)=NULL;


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

template<typename real_t, int sdim>
void
get_slotted_cylinder_atT(const real_t* points_pos,
                         int num_points,
                         real_t* out) {
  real_t xc = 0.5;
  real_t yc = 0.5;
  real_t zc = 0.5;
  real_t R  = 0.3;
  real_t w  = 0.1;
  real_t a  = tcurr;
  tbslas::get_slotted_cylinder<real_t, sdim>(points_pos,
                                             num_points,
                                             out,
                                             xc,
                                             yc,
                                             zc,
                                             R,
                                             w,
                                             a);
}

template <class Real_t>
void get_guassian_kernel_wraper(const Real_t* coord,
                     int n,
                     Real_t* out) {
  const Real_t xc  = 0.7;
  const Real_t yc  = 0.7;
  const Real_t zc  = 0.7;
  tbslas::gaussian_kernel(coord, n, out, xc, yc, zc);
}

template <class Real_t>
void get_multiple_guassian_kernel_wraper(const Real_t* coord,
                     int n,
                     Real_t* out) {
  // FIRST GAUSSIAN
  const Real_t xc1  = 0.7;
  const Real_t yc1  = 0.7;
  const Real_t zc1  = 0.7;
  std::vector<Real_t> out1(n);
  tbslas::gaussian_kernel(coord, n, out1.data(), xc1, yc1, zc1);
  // FIRST GAUSSIAN
  const Real_t xc2  = 0.3;
  const Real_t yc2  = 0.3;
  const Real_t zc2  = 0.3;
  std::vector<Real_t> out2(n);
  tbslas::gaussian_kernel(coord, n, out2.data(), xc2, yc2, zc2);
  // FIRST GAUSSIAN
  const Real_t xc3  = 0.3;
  const Real_t yc3  = 0.3;
  const Real_t zc3  = 0.7;
  std::vector<Real_t> out3(n);
  tbslas::gaussian_kernel(coord, n, out3.data(), xc3, yc3, zc3);
  for (int i = 0; i < n; i++) {
    out[i] = out1[i] + out2[i] + out3[i];
  }
}

template <class Real_t>
void get_hopf_field_wrapper(const Real_t* coord,
                    int n,
                    Real_t* out) {
  const Real_t xc = 0.5;
  const Real_t yc = 0.5;
  const Real_t zc = 0.5;
  tbslas::get_hopf_field(coord, n, out, xc, yc, zc);
}

template <class Real_t>
void get_taylor_green_field_wrapper(const Real_t* coord,
                            int n,
                            Real_t* out) {
  tbslas::get_taylor_green_field(coord, n, out, 0.5, 0.5, 0.5);
}

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
    int max_depth_vel = 0;
    int max_depth_con = 0;
    pvfmm::BoundaryType bc;
    switch(test) {
      // case 1:
      //   fn_vel = tbslas::get_vorticity_field<double,3>;
      //   fn_con = get_gaussian_field_cylinder_atT<double,3>;
      //   bc = pvfmm::FreeSpace;
      //   break;
      // case 2:
      //   fn_vel = tbslas::get_vorticity_field<double,3>;
      //   fn_con = get_slotted_cylinder_atT<double,3>;
      //   bc = pvfmm::FreeSpace;
      //   break;
      // case 3:
      //   fn_vel = tbslas::get_vel_field_hom_y<double,3>;
      //   fn_con = get_gaussian_field_cylinder_atT<double,3>;
      //   bc = pvfmm::Periodic;
      //   break;
      // case 4:                     // regular V, regular C
      //   max_depth_vel=6;
      //   max_depth_con=6;
      //   fn_vel = tbslas::get_vel_field_hom_x<double,3>;
      //   fn_con = tbslas::get_linear_field_y<double,3>;
      //   bc = pvfmm::Periodic;
      //   break;
      // case 5:         // regular V, irregular C
      //   max_depth_vel=6;
      //   fn_vel = tbslas::get_vorticity_field<double,3>;
      //   fn_con = get_guassian_kernel_wraper<double>;
      //   bc = pvfmm::FreeSpace;
      //   break;
      // case 6:         // irregular V, irregular C
      //   fn_vel = get_hopf_field_wrapper<double>;
      //   fn_con = get_guassian_kernel_wraper<double>;
      //   bc = pvfmm::FreeSpace;
      //   break;
      // case 7:  // scaling test case -> uniform fields
      //   fn_vel = tbslas::get_vel_field_hom_x<double,3>;
      //   fn_con = tbslas::get_linear_field_y<double,3>;
      //   bc = pvfmm::Periodic;
      //   break;
      case 1:
        fn_vel = get_taylor_green_field_wrapper<double>;
        fn_con = get_guassian_kernel_wraper<double>;
        bc = pvfmm::Periodic;
        break;
      // case 9:
      //   fn_vel = get_taylor_green_field_wrapper<double>;
      //   fn_con = get_multiple_guassian_kernel_wraper<double>;
      //   bc = pvfmm::Periodic;
      //   break;
    }

    // =========================================================================
    // SIMULATION PARAMETERS
    // =========================================================================
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
    tbslas::NodeFieldFunctor<double,Tree_t> tvel_func = tbslas::NodeFieldFunctor<double,Tree_t>(&tvel);

    // =========================================================================
    // INIT THE CONCENTRATION TREE
    // =========================================================================
    tcurr = 0;
    Tree_t tcon(comm);
    tbslas::ConstructTree<Tree_t>(sim_config->tree_num_point_sources,
                                  sim_config->tree_num_points_per_octanct,
                                  sim_config->tree_chebyshev_order,
                                  max_depth_con?max_depth_con:sim_config->tree_max_depth,
                                  sim_config->tree_adap,
                                  sim_config->tree_tolerance,
                                  comm,
                                  fn_con,
                                  1,
                                  tcon);

    if (sim_config->vtk_save_rate) {
      tvel.Write2File(tbslas::GetVTKFileName(0, "vel").c_str(),
                      sim_config->vtk_order);

      tcon.Write2File(tbslas::GetVTKFileName(0, sim_config->vtk_filename_variable).c_str(),
                      sim_config->vtk_order);
    }

    // =========================================================================
    // RUN
    // =========================================================================
    // set the input_fn to NULL -> needed for adaptive refinement
    std::vector<Node_t*>  ncurr_list = tcon.GetNodeList();
    for(int i = 0; i < ncurr_list.size(); i++) {
      ncurr_list[i]->input_fn = (void (*)(const double* , int , double*))NULL;
    }

    double in_al2,in_rl2,in_ali,in_rli;
    CheckChebOutput<Tree_t>(&tcon,
                            fn_con,
                            1,
                            in_al2,in_rl2,in_ali,in_rli,
                            std::string("Input"));


    int con_noct_sum = 0;
    int con_noct_max = 0;
    int vel_noct_sum = 0;
    int vel_noct_max = 0;
    int timestep = 1;
    int total_ts = sim_config->total_num_timestep;
    for (; timestep < total_ts/2+1; timestep++) {

        // =====================================================================
        // (SEMI) MERGE TO FIX IMBALANCE
        // =====================================================================
        switch(merge) {
          case 2:
            pvfmm::Profile::Tic("CMerge", &sim_config->comm, false, 5);
            tbslas::MergeTree(tvel, tcon);
            pvfmm::Profile::Toc();
            break;
          case 3:
            pvfmm::Profile::Tic("SMerge", &sim_config->comm, false, 5);
            tbslas::SemiMergeTree(tvel, tcon);
            pvfmm::Profile::Toc();
            break;
        }

        // =====================================================================
        // ESTIMATE THE PROBLEM SIZE -> NUMBER OF TREES' OCTANTS
        // =====================================================================
        // if (sim_config->profile) {
        //   int con_noct = tbslas::CountNumLeafNodes(tcon);
        //   con_noct_sum += con_noct;
        //   if (con_noct > con_noct_max) con_noct_max = con_noct;
        //   int vel_noct = tbslas::CountNumLeafNodes(tvel);
        //   vel_noct_sum += vel_noct;
        //   if (vel_noct > vel_noct_max) vel_noct_max = vel_noct;
        // }

        pvfmm::Profile::Tic(std::string("Solve_TN" + tbslas::ToString(static_cast<long long>(timestep))).c_str(), &comm, true);
        {
          // =====================================================================
          // SOLVE SEMILAG
          // =====================================================================
          pvfmm::Profile::Tic("SLM", &sim_config->comm, false, 5);
          tbslas::SolveSemilagInSitu(tvel_func,
                                     tcon,
                                     timestep,
                                     sim_config->dt,
                                     sim_config->num_rk_step);
          pvfmm::Profile::Toc();
        }
        pvfmm::Profile::Toc();        // solve

        // =====================================================================
        // REFINE TREE
        // =====================================================================
        pvfmm::Profile::Tic("RefineTree", &sim_config->comm, false, 5);
        tcon.RefineTree();
        pvfmm::Profile::Toc();

        pvfmm::Profile::Tic("Balance21", &sim_config->comm, false, 5);
        tcon.Balance21(sim_config->bc);
        pvfmm::Profile::Toc();

        //TODO: ONLY FOR STEADY VELOCITY TREES
        tvel.RefineTree();

      // ======================================================================
      // Write2File
      // ======================================================================
      if (sim_config->vtk_save_rate) {
        if (timestep % sim_config->vtk_save_rate == 0)
          tcon.Write2File(tbslas::GetVTKFileName(timestep, sim_config->vtk_filename_variable).c_str(), sim_config->vtk_order);
      }
    }  // end for

    // convert velocity tree
    {
      Node_t* n_curr = tvel.PostorderFirst();
      while (n_curr != NULL) {
        if(!n_curr->IsGhost() && n_curr->IsLeaf()) break;
        n_curr = tvel.PostorderNxt(n_curr);
      }
      while (n_curr != NULL) {
        if (n_curr->IsLeaf() && !n_curr->IsGhost()) {
          pvfmm::Vector<double>& cheb_data = n_curr->ChebData();
          for (int i = 0; i < cheb_data.Dim(); i++)
            cheb_data[i] = -1*cheb_data[i];
        }
        n_curr = tvel.PostorderNxt(n_curr);
      }
    }

    for (; timestep < total_ts+1; timestep++) {
      // =====================================================================
      // (SEMI) MERGE TO FIX IMBALANCE
      // =====================================================================
      switch(merge) {
        case 2:
          pvfmm::Profile::Tic("CMerge", &sim_config->comm, false, 5);
          tbslas::MergeTree(tvel, tcon);
          pvfmm::Profile::Toc();
          break;
        case 3:
          pvfmm::Profile::Tic("SMerge", &sim_config->comm, false, 5);
          tbslas::SemiMergeTree(tvel, tcon);
          pvfmm::Profile::Toc();
          break;
      }

      pvfmm::Profile::Tic(std::string("Solve_TN" + tbslas::ToString(static_cast<long long>(timestep))).c_str(), &comm, true);
      {
        // =====================================================================
        // SOLVE SEMILAG
        // =====================================================================
        pvfmm::Profile::Tic(std::string("SL_TN" + tbslas::ToString(static_cast<long long>(timestep))).c_str(), &sim_config->comm, false, 5);
        tbslas::SolveSemilagInSitu(tvel_func,
                                   tcon,
                                   timestep,
                                   sim_config->dt,
                                   sim_config->num_rk_step);
        pvfmm::Profile::Toc();

        // =====================================================================
        // REFINE TREE
        // =====================================================================
        pvfmm::Profile::Tic("RefineTree", &sim_config->comm, false, 5);
        tcon.RefineTree();
        pvfmm::Profile::Toc();

        pvfmm::Profile::Tic("Balance21", &sim_config->comm, false, 5);
        tcon.Balance21(sim_config->bc);
        pvfmm::Profile::Toc();
      }
      pvfmm::Profile::Toc();        // solve

      //TODO: ONLY FOR STEADY VELOCITY TREES
      tvel.RefineTree();

      // ======================================================================
      // Write2File
      // ======================================================================
      if (sim_config->vtk_save_rate) {
        if (timestep % sim_config->vtk_save_rate == 0)
          tcon.Write2File(tbslas::GetVTKFileName(timestep, sim_config->vtk_filename_variable).c_str(),
                          sim_config->vtk_order);
      }
    }  // end for

    // =========================================================================
    // COMPUTE ERROR
    // =========================================================================
    double al2,rl2,ali,rli;
    CheckChebOutput<Tree_t>(&tcon,
                            fn_con,
                            1,
                            al2,rl2,ali,rli,
                            std::string("Output"));

    // =========================================================================
    // REPORT RESULTS
    // =========================================================================
    int tcon_max_depth=0;
    int tvel_max_depth=0;
    tbslas::GetTreeMaxDepth(tcon, tcon_max_depth);
    tbslas::GetTreeMaxDepth(tvel, tvel_max_depth);

    typedef tbslas::Reporter<double> Rep;
    if(!myrank) {
      Rep::AddData("NP", np, tbslas::REP_INT);
      Rep::AddData("OMP",  sim_config->num_omp_threads, tbslas::REP_INT);

      Rep::AddData("TOL", sim_config->tree_tolerance);
      Rep::AddData("Q", sim_config->tree_chebyshev_order, tbslas::REP_INT);

      Rep::AddData("MaxD", sim_config->tree_max_depth, tbslas::REP_INT);
      Rep::AddData("CMaxD", tcon_max_depth, tbslas::REP_INT);
      Rep::AddData("VMaxD", tvel_max_depth, tbslas::REP_INT);

      Rep::AddData("CBC", sim_config->use_cubic?1:0, tbslas::REP_INT);
      Rep::AddData("CUF", sim_config->cubic_upsampling_factor, tbslas::REP_INT);

      Rep::AddData("DT", sim_config->dt);
      Rep::AddData("TN", sim_config->total_num_timestep, tbslas::REP_INT);
      Rep::AddData("TEST", test, tbslas::REP_INT);
      Rep::AddData("MERGE", merge, tbslas::REP_INT);

      Rep::AddData("InAL2", in_al2);
      Rep::AddData("OutAL2", al2);

      Rep::AddData("InRL2", in_rl2);
      Rep::AddData("OutRL2", rl2);

      Rep::AddData("InALINF", in_ali);
      Rep::AddData("OutALINF", ali);

      Rep::AddData("InRLINF", in_rli);
      Rep::AddData("OutRLINF", rli);

      Rep::AddData("CNOCT", con_noct_sum/(sim_config->total_num_timestep+1),
                   tbslas::REP_INT);
      Rep::AddData("CMaxNOCT", con_noct_max, tbslas::REP_INT);

      Rep::AddData("VNOCT", vel_noct_sum/(sim_config->total_num_timestep+1),
                   tbslas::REP_INT);
      Rep::AddData("VMaxNOCT", vel_noct_max, tbslas::REP_INT);

      Rep::Report();
    }

    //Output Profiling results.
    pvfmm::Profile::print(&comm);
  }

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
