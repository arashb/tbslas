mkdir -p build; rm -rf build/*;
cd build;
cmake .. \
-DCMAKE_C_COMPILER=clang \
-DCMAKE_CXX_COMPILER=clang++; \
# -DCMAKE_PREFIX_PATH=~/.local
