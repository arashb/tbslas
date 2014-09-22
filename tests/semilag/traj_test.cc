// Copyright (C) 2014 by Arash Bakhtiari

#include <vector>

#include <gtest/gtest.h>

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
    size_t tN = dN*dN*dN;
    VecD pnts_pos = tbslas::gen_reg_grid_points<double,3>(dN);
    VecD pnts_vls = tbslas::vorticity_field(pnts_pos );
    vecfd.init(pnts_pos, pnts_vls);
  }

  virtual void TearDown() {
  }

  VFieldD vecfd;
  size_t dN;
};

TEST_F(TrajTest, RK2) {
  std::vector<double> xinit;
  xinit.push_back(0.75);
  xinit.push_back(0.5);
  xinit.push_back(0.5);

  tbslas::CubicInterpPolicy<double> cubic_interp_policy;
  std::vector<double> xsol(xinit.size());
  double dx = 1.0/(dN-1);
  double dt = dx;
  double tf = 2*PI;
  int tn = tf/dt;

  traj_rk2(vecfd, xinit, 0.0, tf, tn, cubic_interp_policy, xsol);
  std::cout << xsol[0] << " " << xsol[1] << " " << xsol[2] << std::endl;
  ASSERT_NEAR(xinit[0], xsol[0], 0.001);
  ASSERT_NEAR(xinit[1], xsol[1], 0.001);
  ASSERT_NEAR(xinit[2], xsol[2], 0.001);
}
