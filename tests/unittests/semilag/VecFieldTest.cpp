#include <utils.hpp>
#include <VecField.hpp>

#include <gtest/gtest.h>

TEST(VecFieldTest, GenerateRegGridPoints) {
  size_t N = 10;
  std::vector<double> coord = semilag::gen_reg_grid_points<double,3>(N);
}

TEST(VecFieldTest, Constructor) {
  typedef semilag::VecField<double,3,3> VecFD;
  size_t dN = 11;
  size_t tN = dN*dN*dN;
  VecFD vecf1;
  VecFD vecf2(semilag::gen_reg_grid_points<double,3>(dN),
              semilag::gen_reg_grid_points<double,3>(dN),
              tN);
}

TEST(VecFieldTest, WriteToFile) {
  typedef semilag::VecField<double,3,3> VecFD;
  size_t dN = 11;
  size_t tN = dN*dN*dN;
  VecFD vecf2(semilag::gen_reg_grid_points<double,3>(dN),
              semilag::gen_reg_grid_points<double,3>(dN),
              tN);
  vecf2.write2file("testfile");
}
