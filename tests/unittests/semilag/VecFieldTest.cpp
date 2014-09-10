// Copyright (C) 2014 by Arash Bakhtiari
#include <VecField.hpp>
#include <CubicInterpolant.hpp>
#include <utils.hpp>
#include <gtest/gtest.h>
#include <vector>

typedef tbslas::VecField<double,3,3> VFieldD;
typedef std::vector<double> VecD;

class VecFieldTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    size_t dN = 11;
    size_t tN = dN*dN*dN;
    VecD pnts_pos = tbslas::gen_reg_grid_points<double>(dN);
    VecD pnts_vls = tbslas::vorticity_field(pnts_pos);
    vecfd.init(pnts_pos, pnts_vls);
  }

  virtual void TearDown() {
  }

  VFieldD vecfd;
};

TEST_F(VecFieldTest, WriteToFile) {
  vecfd.write2file("VecFieldTest-WriteToFile");
}

TEST_F(VecFieldTest, Interpolate) {
  VecD qry_pnts;
  qry_pnts.push_back(0.75);
  qry_pnts.push_back(0.75);
  qry_pnts.push_back(0.50);
  VecD expected_res = tbslas::vorticity_field(qry_pnts);

  tbslas::CubicInterpolant<double> intrplnt;
  VecD actual_res(qry_pnts.size());
  vecfd.interp(qry_pnts, intrplnt, actual_res);

  ASSERT_DOUBLE_EQ(expected_res[0], actual_res[0]);
  ASSERT_DOUBLE_EQ(expected_res[1], actual_res[1]);
  ASSERT_DOUBLE_EQ(expected_res[2], actual_res[2]);
}
