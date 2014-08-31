#include <semilag.hpp>

int main(int argc, char *argv[])
{
  size_t N = 10;
  std::vector<double> coord;
  coord = semilag::gen_reg_grid_points<double>(N);
  return 0;
}
