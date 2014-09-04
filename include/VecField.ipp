// Copyright (C) 2014 by Arash Bakhtiari

// *************************************************************************
// You may not use this file except in compliance with the License.
// You obtain a copy of the License in the LICENSE file.

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// *************************************************************************

#include <stdint.h>
#include <string>       // std::string
#include <vector>
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream
#include <fstream>      // std::ofstream
#include <cmath>        // pow()

namespace semilag {

  template<typename real_t, int sdim, int vdim>
  VecField<real_t, sdim, vdim>::VecField() {}

  template<typename real_t, int sdim, int vdim>
  VecField<real_t, sdim, vdim>::VecField(std::vector<real_t> pos,
                                         std::vector<real_t> vals,
                                         size_t size):
    _pos(pos),
    _vals(vals),
    _size(size)
  {};


  inline bool isLittleEndian() {
    int n = 1;
    return (*(char *)&n == 1);
  }

  template<typename real_t, int sdim, int vdim>
  VecField<real_t, sdim, vdim>::~VecField() {}

  template<typename real_t, int sdim, int vdim>
  void VecField<real_t, sdim, vdim>::write2file(const char* fname) {
    int pnt_cnt = _pos.size()/sdim;
    // FIXME: is correct only for the regular grid
    float dN = pow(pnt_cnt, 1./(sdim));
    int cll_cnt = pow(dN-1, sdim);
    int myrank = 0;

    // open file for writing.
    std::stringstream vtufname;
    vtufname << fname << std::setfill('0') << std::setw(6) << myrank << ".vts";
    std::ofstream vtsfile;
    vtsfile.open(vtufname.str().c_str());
    if (vtsfile.fail()) return;

    // write to file.
    size_t data_size = 0;
    vtsfile << "<?xml version=\"1.0\"?>\n";
    vtsfile << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" ";
    if (isLittleEndian())
      vtsfile<< "byte_order=\"LittleEndian\">\n";
    else
      vtsfile<< "byte_order=\"BigEndian\">\n";
    //**************************************************************************
    std::stringstream extend;
    extend << "0 " << dN-1 << " 0 " << dN-1 << " 0 " << dN-1;
    vtsfile <<"  <StructuredGrid "
            <<"WholeExtent=\""<< extend.str() << "\""">\n";
    vtsfile <<"    <Piece "
            <<"Extent=\""<< extend.str() << "\""
            <<">\n";
    //**************************************************************************
    vtsfile <<"      <Points>\n";
    vtsfile <<"        <DataArray "
            <<"type=\"Float" << sizeof(real_t)*8 << "\" "
            <<"NumberOfComponents=\"" << sdim << "\" "
            <<"Name=\"" << "Position" << "\" "
            <<"format=\"" << "appended" << "\" "
            <<"offset=\"" << data_size << "\" />\n";
    vtsfile <<"      </Points>\n";
    int32_t points_block_size = _pos.size()*sizeof(real_t);
    data_size += sizeof(uint32_t)+ points_block_size;
    //**************************************************************************
    vtsfile <<"      <PointData>\n";
    vtsfile <<"        <DataArray type=\"Float" << sizeof(real_t)*8 << "\" "
            <<"NumberOfComponents=\"" << vdim << "\" "
            <<"Name=\"Value\" "
            <<"format=\"appended\" "
            <<"offset=\"" << data_size << "\" />\n";
    vtsfile <<"      </PointData>\n";
    int32_t vals_block_size = _vals.size()*sizeof(real_t);
    data_size += sizeof(uint32_t) + vals_block_size;
    //**************************************************************************
    vtsfile <<"      <CellData>\n";
    vtsfile <<"      </CellData>\n";
    //**************************************************************************
    vtsfile <<"    </Piece>\n";
    vtsfile <<"  </StructuredGrid>\n";
    //**************************************************************************
    vtsfile <<"  <AppendedData encoding=\"raw\">\n";
    vtsfile <<"    _";
    // points
    vtsfile.write((char*)&points_block_size, sizeof(int32_t));
    vtsfile.write((char*)&_pos[0], points_block_size);
    // points data
    vtsfile.write((char*)&vals_block_size, sizeof(int32_t));
    vtsfile.write((char*)&_vals[0], vals_block_size);
    vtsfile <<"\n";
    vtsfile <<"  </AppendedData>\n";
    //**************************************************************************

    vtsfile <<"</VTKFile>\n";
    vtsfile.close();
  }

}  // namespace semilag
