// Copyright (C) 2014 by Arash Bakhtiari

#include <vector>

#include <gtest/gtest.h>

#include "semilag/vec_field.h"
#include "semilag/cubic_interp_policy.h"
#include "semilag/utils.h"

typedef slas::VecField<double,3,3> VFieldD;
typedef std::vector<double> VecD;

class VecFieldTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    size_t dN = 11;
    size_t tN = dN*dN*dN;
    VecD points_positions = slas::generate_reg_grid_points<double,3>(dN);
    VecD points_values    = slas::generate_vorticity_field(points_positions);
    vecfd.init(points_positions, points_values);
    vecfd.push_back_values(points_values);
  }

  virtual void TearDown() {
  }

  VFieldD vecfd;
};

TEST_F(VecFieldTest, WriteToFile) {
  vecfd.write2file("vecfieldtest_writetofile_");
}

TEST_F(VecFieldTest, Interpolate) {
  VecD query_points;
  query_points.push_back(0.75);
  query_points.push_back(0.75);
  query_points.push_back(0.50);
  VecD expected_res = slas::generate_vorticity_field(query_points);

  slas::CubicInterpPolicy<double> cubic_interp_policy;
  VecD actual_res(query_points.size());
  vecfd.interp(query_points, cubic_interp_policy, actual_res);

  ASSERT_DOUBLE_EQ(expected_res[0], actual_res[0]);
  ASSERT_DOUBLE_EQ(expected_res[1], actual_res[1]);
  ASSERT_DOUBLE_EQ(expected_res[2], actual_res[2]);
}

TEST_F(VecFieldTest, Save) {
  vecfd.save("vecfieldtest_save_");
}
