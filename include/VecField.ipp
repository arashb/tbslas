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
#include <assert.h>
#include <string>       // std::string
#include <vector>
#include <iomanip>
#include <sstream>      // std::stringstream
#include <fstream>      // std::ofstream
#include <iostream>     // std::cout
#include <cmath>        // pow()


namespace tbslas {

template<typename real_t>
inline void
validate_field_data(std::vector<real_t>& pnts,
                    int sdim,
                    std::vector<real_t>& vls,
                    int vdim) {
  assert((pnts.size()/sdim)  == (vls.size()/vdim));
}

inline bool isLittleEndian() {
  int n = 1;
  return (*reinterpret_cast<char*>(&n) == 1);
}

template<typename real_t, int sdim, int vdim>
VecField<real_t, sdim, vdim>::VecField() {}


template<typename real_t, int sdim, int vdim>
VecField<real_t, sdim, vdim>::VecField(std::vector<real_t> pnts,
                                       std::vector<real_t> vls):
    _pnts(pnts),
    _vls(vls) {
  validate_field_data(_pnts, sdim, _vls, vdim);
}

template<typename real_t, int sdim, int vdim>
void
VecField<real_t, sdim, vdim>::init(std::vector<real_t> pnts,
                                   std::vector<real_t> vls) {
  _pnts = pnts;
  _vls = vls;
  validate_field_data(_pnts, sdim, _vls, vdim);
}



template<typename real_t, int sdim, int vdim>
template<typename C>  // Interpolant class
void
VecField<real_t, sdim, vdim>::interp(std::vector<real_t>& qry_pnts,
                                     C& interpolant,
                                     std::vector<real_t>& qry_vls) {
  interpolant.interp(_pnts, sdim, _vls, vdim, qry_pnts, qry_vls);
}

template<typename real_t, int sdim, int vdim>
VecField<real_t, sdim, vdim>::~VecField() {}

template<typename real_t, int sdim, int vdim>
void
VecField<real_t, sdim, vdim>::write2file(const char* fname) {
  int pnt_cnt = _pnts.size()/sdim;
  // FIXME: is correct only for the regular grid
  float dN = pow(pnt_cnt, 1./(sdim));
  // int cll_cnt = pow(dN-1, sdim);
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
  int32_t points_block_size = _pnts.size()*sizeof(real_t);
  data_size += sizeof(uint32_t)+ points_block_size;
  //**************************************************************************
  vtsfile <<"      <PointData vectors=\"values\" >\n";
  vtsfile <<"        <DataArray type=\"Float" << sizeof(real_t)*8 << "\" "
          <<"Name=\"values\" "
          <<"NumberOfComponents=\"" << vdim << "\" "
          <<"format=\"appended\" "
          <<"offset=\"" << data_size << "\" />\n";
  vtsfile <<"      </PointData>\n";
  int32_t vls_block_size = _vls.size()*sizeof(real_t);
  data_size += sizeof(uint32_t) + vls_block_size;
  //**************************************************************************
  // vtsfile <<"      <CellData>\n";
  // vtsfile <<"      </CellData>\n";
  //**************************************************************************
  vtsfile <<"    </Piece>\n";
  vtsfile <<"  </StructuredGrid>\n";
  //**************************************************************************
  vtsfile <<"  <AppendedData encoding=\"raw\">\n";
  vtsfile <<"    _";

  // points positions
  vtsfile.write(reinterpret_cast<char*>(&points_block_size),
                sizeof(points_block_size));
  for (int i = 0; i < pnt_cnt; i++) {
    for(int j = 0; j < sdim; j++)
      vtsfile.write(reinterpret_cast<char*>(&_pnts[i+j*pnt_cnt]),
                    sizeof(real_t));
  }
  // vtsfile.write(reinterpret_cast<char*>(&_pnts[0]),
  //               points_block_size);

  // points data
  // see the data memmory layout in vector field class
  vtsfile.write(reinterpret_cast<char*>(&vls_block_size),
                sizeof(vls_block_size));
  for (int i = 0; i < pnt_cnt; i++) {
    for(int j = 0; j < vdim; j++)
      vtsfile.write(reinterpret_cast<char*>(&_vls[i+j*pnt_cnt]),
                    sizeof(real_t));
  }
  // vtsfile.write(reinterpret_cast<char*>(&_vls[0]),
  //               vls_block_size);
  vtsfile <<"\n";
  vtsfile <<"  </AppendedData>\n";
  //**************************************************************************
  vtsfile <<"</VTKFile>\n";
  vtsfile.close();
}

}  // namespace tbslas
