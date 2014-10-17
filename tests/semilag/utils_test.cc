// Copyright (C) 2014 by Arash Bakhtiari

#include <gtest/gtest.h>

#include "semilag/utils.h"

TEST(UtilsTest, FindGridIndex) {
  std::vector<double> grid = tbslas::generate_reg_grid_points<double,1>(11);
  int index                = tbslas::find_grid_index_1d(grid, 0.55);
  ASSERT_EQ(5, index);
}
