// Copyright (C) 2014 by Arash Bakhtiari

#include <vector>

#include <gtest/gtest.h>

#include "semilag/common.h"
#include "semilag/utils.h"
#include "semilag/vec_field.h"
#include "semilag/traj.h"

typedef tbslas::VecField<double,3,3> VFieldD;
typedef std::vector<double> VecD;

#define PI 3.14159

class TrajTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    dN = 65;
    dx = 1.0/(dN-1);
    dt = dx;
    tf = 2*PI;
    tn = tf/dt;
    size_t tN = dN*dN*dN;

    VecD pnts_pos = tbslas::generate_reg_grid_points<double,3>(dN);
    VecD pnts_vls = tbslas::generate_vorticity_field(pnts_pos);
    vecfd.init(pnts_pos, pnts_vls);

    // initialization of constant velocity field
    for ( int i = 1; i < tn+1; i++)
      vecfd.push_back_values(pnts_vls, i*dt);
  }

  virtual void TearDown() {
  }

  VFieldD vecfd;
  size_t dN;
  double dx,dt,tf;
  int tn;
};

TEST_F(TrajTest, RK2SteadyVorticityField) {
  std::vector<double> xinit;
  xinit.push_back(0.75);
  xinit.push_back(0.5);
  xinit.push_back(0.5);

  tbslas::CubicInterpPolicy<double> cubic_interp_policy;
  std::vector<double> xsol(xinit.size());

  traj_rk2(vecfd, xinit, 0.0, tf, tn, cubic_interp_policy, xsol);

  // std::cout << xsol[0] << " " << xsol[1] << " " << xsol[2] << std::endl;
  ASSERT_NEAR(xinit[0], xsol[0], 0.001);
  ASSERT_NEAR(xinit[1], xsol[1], 0.001);
  ASSERT_NEAR(xinit[2], xsol[2], 0.001);
}


TEST_F(TrajTest, RK2SteadyVorticityFunctor) {
  std::vector<double> xinit;
  xinit.push_back(0.75);
  xinit.push_back(0.5);
  xinit.push_back(0.5);

  std::vector<double> xsol(xinit.size());

  tbslas::traj_rk2(tbslas::get_vorticity_field<double,3>,
                   xinit,
                   0.0,
                   tf,
                   tn,
                   xsol);

  // FIXME: somehow the asserts do not work precisely
  ASSERT_NEAR(xinit[0], xsol[0], 0.001);
  ASSERT_NEAR(xinit[1], xsol[1], 0.001);
  ASSERT_NEAR(xinit[2], xsol[2], 0.001);
}
