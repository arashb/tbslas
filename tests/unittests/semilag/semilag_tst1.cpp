#include <semilag.hpp>
#include <gtest/gtest.h>

TEST(SemilagTest, GenerateRegGridPoints){
  size_t N = 10;
  std::vector<double> coord;
  coord = semilag::gen_reg_grid_points<double>(N);
}
