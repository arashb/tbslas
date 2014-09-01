mkdir -p build; rm -rf build/*;
cd build;
cmake .. \
-DCMAKE_C_COMPILER=clang \
-DCMAKE_CXX_COMPILER=clang++ \
-DBUILD-GTEST-UNITTEST=ON \
# -DCMAKE_PREFIX_PATH=~/.local
