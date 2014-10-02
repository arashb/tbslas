mkdir -p build; rm -rf build/*;
cd build;
cmake .. \
-DBUILD-GTEST-UNITTEST=ON
# -DCMAKE_C_COMPILER=clang \
# -DCMAKE_CXX_COMPILER=clang++ \

# -DCMAKE_PREFIX_PATH=~/.local
