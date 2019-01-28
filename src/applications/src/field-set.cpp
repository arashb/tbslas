// *************************************************************************
// Copyright (C) 2015 by Arash Bakhtiari
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
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <cheb_node.hpp>
#include <cheb_utils.hpp>
#include <mpi_tree.hpp>
#include <profile.hpp>
#include <pvfmm_common.hpp>
#include <vector.hpp>

#include <field_wrappers.h>
#include <utils.hpp>

#include <utils/common.h>
#include <utils/metadata.h>
#include <utils/reporter.h>

#include <tree/tree_semilag.h>
#include <tree/tree_utils.h>

#include <tree/tree_set_functor.h>

typedef pvfmm::Cheb_Node<double> Node_t;
typedef pvfmm::MPI_Tree<Node_t> Tree_t;

typedef tbslas::MetaData<std::string, std::string, std::string> MetaData_t;
double tcurr = 0;

void (*fn_1)(const double*, int, double*) = NULL;
void (*fn_2)(const double*, int, double*) = NULL;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int np;
  MPI_Comm_size(comm, &np);
  int myrank;
  MPI_Comm_rank(comm, &myrank);

  parse_command_line_options(argc, argv);

  int test =
      strtoul(commandline_option(
                  argc, argv, "-test", "1", false,
                  "-test <int> = (1)    : 1) Gaussian profile 2) Zalesak disk"),
              NULL, 10);

  {
    tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
    tbslas::new_nodes<Tree_t::Real_t>(sim_config->tree_chebyshev_order, 3);

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
    double data_dof = 0;
    pvfmm::BoundaryType bc;
    switch (test) {
      case 1:

        // fn_2 = tbslas::get_linear_field_y<double,3>;
        // data_dof = 1;

        // fn_2 = get_gaussian_field_atT<double,3>;
        // data_dof = 1;

        // fn_2 = get_vorticity_field_tv_wrapper<double>;
        // data_dof = 3;

        fn_2 = get_taylor_green_field_tv_ns_wrapper<double>;
        data_dof = 3;

        // fn_2 = tbslas::get_vorticity_field<double,3>;
        // data_dof = 3;

        // fn_2 = get_hopf_field_wrapper<double>;
        // data_dof = 3;

        // fn_2 = get_taylor_green_field_tv_wrapper<double>;
        // data_dof = 3;

        bc = pvfmm::Periodic;
        // bc = pvfmm::FreeSpace;

        break;
    }

    // =========================================================================
    // SIMULATION PARAMETERS
    // =========================================================================
    sim_config->vtk_filename_variable = "conc";
    sim_config->bc = bc;

    double DT = sim_config->dt;
    // double DT = 0.3925;
    // double DT = 1;

    double cheb_deg = sim_config->tree_chebyshev_order;

    std::vector<double> tree_set_times;
    std::vector<Tree_t*> tree_set_elems;

    // =========================================================================
    // INIT THE TREE 1
    // =========================================================================
    tcurr = 0;
    Tree_t tree1(comm);
    tbslas::ConstructTree<Tree_t>(
        sim_config->tree_num_point_sources,
        sim_config->tree_num_points_per_octanct,
        sim_config->tree_chebyshev_order, sim_config->tree_max_depth,
        sim_config->tree_adap, sim_config->tree_tolerance, comm, fn_2, data_dof,
        tree1);
    if (sim_config->vtk_save_rate) {
      tree1.Write2File(tbslas::GetVTKFileName(1, "tree").c_str(),
                       sim_config->vtk_order);
    }
    tree_set_times.push_back(tcurr);
    tree_set_elems.push_back(&tree1);

    // =========================================================================
    // INIT THE TREE 2
    // =========================================================================
    tcurr += DT;
    Tree_t tree2(comm);
    tbslas::ConstructTree<Tree_t>(
        sim_config->tree_num_point_sources,
        sim_config->tree_num_points_per_octanct,
        sim_config->tree_chebyshev_order, sim_config->tree_max_depth,
        sim_config->tree_adap, sim_config->tree_tolerance, comm, fn_2, data_dof,
        tree2);

    if (sim_config->vtk_save_rate) {
      tree2.Write2File(tbslas::GetVTKFileName(2, "tree").c_str(),
                       sim_config->vtk_order);
    }

    tree_set_times.push_back(tcurr);
    tree_set_elems.push_back(&tree2);

    // =========================================================================
    // INIT THE TREE 3
    // =========================================================================
    tcurr += DT;
    Tree_t tree3(comm);
    tbslas::ConstructTree<Tree_t>(
        sim_config->tree_num_point_sources,
        sim_config->tree_num_points_per_octanct,
        sim_config->tree_chebyshev_order, sim_config->tree_max_depth,
        sim_config->tree_adap, sim_config->tree_tolerance, comm, fn_2, data_dof,
        tree3);

    if (sim_config->vtk_save_rate) {
      tree3.Write2File(tbslas::GetVTKFileName(3, "tree").c_str(),
                       sim_config->vtk_order);
    }

    tree_set_times.push_back(tcurr);
    tree_set_elems.push_back(&tree3);

    // =========================================================================
    // INIT THE TREE 4
    // =========================================================================
    tcurr += DT;
    Tree_t tree4(comm);
    tbslas::ConstructTree<Tree_t>(
        sim_config->tree_num_point_sources,
        sim_config->tree_num_points_per_octanct,
        sim_config->tree_chebyshev_order, sim_config->tree_max_depth,
        sim_config->tree_adap, sim_config->tree_tolerance, comm, fn_2, data_dof,
        tree4);
    if (sim_config->vtk_save_rate) {
      tree4.Write2File(tbslas::GetVTKFileName(4, "tree").c_str(),
                       sim_config->vtk_order);
    }

    tree_set_times.push_back(tcurr);
    tree_set_elems.push_back(&tree4);

    // =========================================================================
    // MERGE
    // =========================================================================
    tcurr = 1.5 * DT;
    Tree_t merged_tree(comm);
    tbslas::ConstructTree<Tree_t>(
        sim_config->tree_num_point_sources,
        sim_config->tree_num_points_per_octanct,
        sim_config->tree_chebyshev_order, sim_config->tree_max_depth,
        sim_config->tree_adap, sim_config->tree_tolerance, comm, fn_2, data_dof,
        merged_tree);

    double in_al2, in_rl2, in_ali, in_rli;
    CheckChebOutput<Tree_t>(&merged_tree, fn_2, data_dof, in_al2, in_rl2,
                            in_ali, in_rli, std::string("Input"));

    // =========================================================================
    // COLLECT THE MERGED TREE POINTS
    // =========================================================================
    std::vector<double> merged_tree_points_pos;
    int num_leaf =
        tbslas::CollectChebTreeGridPoints(merged_tree, merged_tree_points_pos);

    // =========================================================================
    // CONSTRUCT THE FUNCTOR
    // =========================================================================
    tbslas::FieldSetFunctor<double, Tree_t> tree_set_functor(tree_set_elems,
                                                             tree_set_times);
    // tbslas::NodeFieldFunctor<double,Tree_t> tree_functor(&merged_tree);

    // int num_points = 1;
    // std::vector<double> xtmp(num_points*3);
    // std::vector<double> vtmp(num_points*data_dof);

    // xtmp[0] = 0.8;
    // xtmp[1] = 1.0;
    // xtmp[2] = 0.3;

    // tree_functor(xtmp.data(),
    //              num_points,
    //              1.5*DT,
    //              vtmp.data());
    // std::cout << "MYRANK: " << myrank << " vals: " << vtmp[0] << " " <<
    // vtmp[1] << " " << vtmp[2] << std::endl;

    // =========================================================================
    // INTERPOLATE IN TIME
    // =========================================================================
    int merged_tree_num_points = merged_tree_points_pos.size() / 3;
    std::vector<double> merged_tree_points_val(merged_tree_num_points *
                                               data_dof);

    pvfmm::Profile::Tic("EvalSet", &sim_config->comm, false, 5);
    tree_set_functor(merged_tree_points_pos.data(), merged_tree_num_points,
                     1.5 * DT, merged_tree_points_val.data());

    // tree_functor(merged_tree_points_pos.data(),
    //              merged_tree_num_points,
    //              merged_tree_points_val.data());
    pvfmm::Profile::Toc();

    // =========================================================================
    // FIX THE VALUES MEMORY LAYOUT
    // =========================================================================
    int d = cheb_deg + 1;
    int num_pnts_per_node = d * d * d;
    std::vector<double> mt_pnts_val_ml(merged_tree_num_points * data_dof);
    for (int nindx = 0; nindx < num_leaf; nindx++) {
      int input_shift = nindx * num_pnts_per_node * data_dof;
      for (int j = 0; j < num_pnts_per_node; j++) {
        for (int i = 0; i < data_dof; i++) {
          mt_pnts_val_ml[input_shift + j + i * num_pnts_per_node] =
              merged_tree_points_val[input_shift + j * data_dof + i];
        }
      }
    }

    // =========================================================================
    // SET INTERPOLATED VALUES
    // =========================================================================
    pvfmm::Profile::Tic("SetValues", &sim_config->comm, false, 5);
    tbslas::SetTreeGridValues(merged_tree, cheb_deg, data_dof, mt_pnts_val_ml);
    pvfmm::Profile::Toc();

    if (sim_config->vtk_save_rate) {
      merged_tree.Write2File(tbslas::GetVTKFileName(0, "tree_interp").c_str(),
                             sim_config->vtk_order);
    }

    double al2, rl2, ali, rli;
    CheckChebOutput<Tree_t>(&merged_tree, fn_2, data_dof, al2, rl2, ali, rli,
                            std::string("Output"));

    // =========================================================================
    // TEST
    // =========================================================================
    // tree_set_functor(xtmp.data(),
    //                  num_points,
    //                  1.5*DT,
    //                  vtmp.data());
    // std::cout << "vals: " << vtmp[0] << " " << vtmp[1] << " " << vtmp[2] <<
    // std::endl;

    // tbslas::NodeFieldFunctor<double,Tree_t> tf(&merged_tree);
    // tf(xtmp.data(),
    //    num_points,
    //    1.5*DT,
    //    vtmp.data());
    // std::cout << "vals: " << vtmp[0] << " " << vtmp[1] << " " << vtmp[2] <<
    // std::endl;

    // =========================================================================
    // REPORT RESULTS
    // =========================================================================
    int tcon_max_depth = 0;
    int tvel_max_depth = 0;
    tbslas::GetTreeMaxDepth(merged_tree, tcon_max_depth);
    // tbslas::GetTreeMaxDepth(tvel, tvel_max_depth);

    typedef tbslas::Reporter<double> Rep;
    if (!myrank) {
      Rep::AddData("NP", np, tbslas::REP_INT);
      Rep::AddData("OMP", sim_config->num_omp_threads, tbslas::REP_INT);

      Rep::AddData("TOL", sim_config->tree_tolerance);
      Rep::AddData("Q", sim_config->tree_chebyshev_order, tbslas::REP_INT);

      Rep::AddData("MaxD", sim_config->tree_max_depth, tbslas::REP_INT);
      Rep::AddData("CMaxD", tcon_max_depth, tbslas::REP_INT);
      // Rep::AddData("VMaxD", tvel_max_depth, tbslas::REP_INT);

      // Rep::AddData("CBC", sim_config->use_cubic?1:0, tbslas::REP_INT);
      // Rep::AddData("CUF", sim_config->cubic_upsampling_factor,
      // tbslas::REP_INT);

      Rep::AddData("DT", sim_config->dt);
      // Rep::AddData("TN", sim_config->total_num_timestep, tbslas::REP_INT);
      // Rep::AddData("TEST", test, tbslas::REP_INT);
      // Rep::AddData("MERGE", merge, tbslas::REP_INT);

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

    // Output Profiling results.
    // pvfmm::Profile::print(&comm);
  }

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
