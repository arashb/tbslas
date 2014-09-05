// Copyright (C) 2014 by Arash Bakhtiari
#include <VecField.hpp>
#include <CubicInterpolant.hpp>
#include <utils.hpp>
#include <gtest/gtest.h>
#include <vector>


TEST(VecFieldTest, GenerateRegGridPoints) {
  size_t N = 10;
  std::vector<double> coord = tbslas::gen_reg_grid_points<double,3>(N);
}

TEST(VecFieldTest, Constructor) {
  typedef tbslas::VecField<double,3,3> VecFD;
  size_t dN = 11;
  size_t tN = dN*dN*dN;
  VecFD vecf1;
  VecFD vecf2(tbslas::gen_reg_grid_points<double,3>(dN),
              tbslas::gen_reg_grid_points<double,3>(dN),
              tN);
}

TEST(VecFieldTest, WriteToFile) {
  typedef double real_t;
  typedef tbslas::VecField<real_t,3,3> VecFD;
  size_t dN = 11;
  size_t tN = dN*dN*dN;
  VecFD vecf2(tbslas::gen_reg_grid_points<real_t,3>(dN),
              tbslas::gen_reg_grid_points<real_t,3>(dN),
              tN);
  vecf2.write2file("testfile");
}

TEST(VecFieldTest, Interpolate) {
  typedef tbslas::VecField<double,3,1> VecFD;
  typedef std::vector<double> vec;

  size_t dN = 4;
  size_t tN = dN*dN*dN;
  vec pnts_vals;
  for (int i = 0; i < tN; i++)
    pnts_vals.push_back(2.0);
  VecFD vecf2(tbslas::gen_reg_grid_points<double,3>(dN),
              pnts_vals,
              tN);

  vec qry_pnts;
  qry_pnts.push_back(0.75);
  qry_pnts.push_back(0.75);
  qry_pnts.push_back(0.50);

  tbslas::CubicInterpolant<double> intrplnt;
  vec res = vecf2.interp< tbslas::CubicInterpolant<double> >(qry_pnts, intrplnt);
  ASSERT_DOUBLE_EQ(2.0, res[0]);
}
