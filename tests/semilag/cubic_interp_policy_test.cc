// Copyright (C) 2014 by Arash Bakhtiari

#include <vector>
#include <cstdlib>
#include <cmath>

#include <gtest/gtest.h>

#include "semilag/cubic_interp_policy.h"
#include "semilag/utils.h"

class CubicInterpPolicyTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    N = 4;
    num_points = (size_t)std::pow(N,sdim);
    num_time_steps = 10;
    dx = 1.0/(N-1);
    dt = dx;

    grid_points = tbslas::gen_reg_grid_points<double,sdim>(N);
    std::vector<double> tmp_val(vdim*num_points);
    // creating the time line and values
    for ( int i = 0; i <= num_time_steps; i++) {
      for (int j = 0; j < num_points; j++) {
        tmp_val[j] = i;
        tmp_val[j+num_points] = i+1;
      }

      time.push_back(i*dt);
      grid_values.push_back(tmp_val);
    }

  }

  virtual void TearDown() {
  }

  static const int sdim = 3;
  static const int vdim = 2;

  size_t N, num_points;
  double dx,dt;
  int num_time_steps;
  std::vector<double> grid_points, time;
  std::vector< std::vector<double> > grid_values;
  tbslas::CubicInterpPolicy<double> interp_policy;
};

TEST_F(CubicInterpPolicyTest, InterpTimeIndependent) {
  int tstep = 1;
  std::vector<double> query_points;
  query_points.push_back(0.35);
  query_points.push_back(0.35);
  query_points.push_back(0.35);

  std::vector<double> query_values(1*vdim);
  interp_policy.interp(grid_points, sdim, grid_values[tstep], vdim,
                       query_points, query_values);
  // std::cout << query_values[0] << query_values[1] << std::endl;
  ASSERT_EQ(1,query_values[0]);
  ASSERT_EQ(2,query_values[1]);
}

TEST_F(CubicInterpPolicyTest, InterpTimeDependent) {
  std::vector<double> query_points;
  query_points.push_back(0.5);
  query_points.push_back(0.5);
  query_points.push_back(0.5);
  double query_time = 7.3*dt;

  std::vector<double> query_values(1*vdim);
  interp_policy.interp(grid_points,
                       sdim,
                       grid_values,
                       vdim,
                       time,
                       query_points,
                       query_time,
                       query_values);
  // std::cout << query_values[0] << query_values[1] << std::endl;
  ASSERT_DOUBLE_EQ(7.3, query_values[0]);
  ASSERT_DOUBLE_EQ(8.3, query_values[1]);
}

TEST_F(CubicInterpPolicyTest, InterpTimeDependentExtremeLeft) {
  std::vector<double> query_points;
  query_points.push_back(0.5);
  query_points.push_back(0.5);
  query_points.push_back(0.5);
  double query_time = 0;

  std::vector<double> query_values(1*vdim);
  interp_policy.interp(grid_points,
                       sdim,
                       grid_values,
                       vdim,
                       time,
                       query_points,
                       query_time,
                       query_values);
  // std::cout << query_values[0] << query_values[1] << std::endl;
  ASSERT_DOUBLE_EQ(0, query_values[0]);
  ASSERT_DOUBLE_EQ(1, query_values[1]);
}

TEST_F(CubicInterpPolicyTest, InterpTimeDependentExtremeRight) {
  std::vector<double> query_points;
  query_points.push_back(0.5);
  query_points.push_back(0.5);
  query_points.push_back(0.5);
  double query_time = num_time_steps*dt;

  std::vector<double> query_values(1*vdim);
  interp_policy.interp(grid_points,
                       sdim,
                       grid_values,
                       vdim,
                       time,
                       query_points,
                       query_time,
                       query_values);
  // std::cout << query_values[0] << query_values[1] << std::endl;
  ASSERT_DOUBLE_EQ(10, query_values[0]);
  ASSERT_DOUBLE_EQ(11, query_values[1]);
}

TEST_F(CubicInterpPolicyTest, InterpTimeDependentExtremeOutOfInterval) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  std::vector<double> query_points;
  query_points.push_back(0.5);
  query_points.push_back(0.5);
  query_points.push_back(0.5);
  double query_time = num_time_steps*dt+0.5*dt;

  std::vector<double> query_values(1*vdim);
  ASSERT_DEATH(
      {interp_policy.interp(grid_points,
                            sdim,
                            grid_values,
                            vdim,
                            time,
                            query_points,
                            query_time,
                            query_values);}
      ,"");
}
