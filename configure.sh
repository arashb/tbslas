mkdir -p build; rm -rf build/*;
cd build;
cmake .. \
-DBUILD-GTEST-UNITTEST=ON \
# -DCMAKE_C_COMPILER=icc \
# -DCMAKE_CXX_COMPILER=icpc \

# -DCMAKE_PREFIX_PATH=~/.local
