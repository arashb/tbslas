// *************************************************************************
// Copyright (C) 2014 by Arash Bakhtiari

// You may not use this file except in compliance with the License.
// You obtain a copy of the License in the LICENSE file.

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// *************************************************************************

#include <stdint.h>
#include <cassert>
#include <string>       // std::string
#include <vector>
#include <iomanip>
#include <sstream>      // std::stringstream
#include <fstream>      // std::ofstream
#include <iostream>     // std::cout
#include <cmath>        // pow()

#include "semilag/utils.h"

namespace tbslas {

template<typename real_t>
inline void
validate_field_data(size_t num_field_points,
                    std::vector<real_t>& field_values,
                    int vdim) {
  assert(num_field_points == (field_values.size()/vdim));
}

template<typename real_t, int sdim, int vdim>
VecField<real_t, sdim, vdim>::VecField() {}

template<typename real_t, int sdim, int vdim>
VecField<real_t, sdim, vdim>::VecField(std::vector<real_t> field_points,
                                       std::vector<real_t> field_values,
                                       real_t time_init):
    field_points_(field_points)
{
  num_field_points_ = field_points_.size()/sdim;
  this->push_back_values(field_values, time_init);
}

template<typename real_t, int sdim, int vdim>
VecField<real_t, sdim, vdim>::~VecField() {}

template<typename real_t, int sdim, int vdim>
void
VecField<real_t, sdim, vdim>::init(std::vector<real_t> field_points,
                                   std::vector<real_t> field_values,
                                   real_t time_init) {
  field_points_ = field_points;
  num_field_points_ = field_points_.size()/sdim;
  field_values_.clear();
  this->push_back_values(field_values, time_init);
}


template<typename real_t, int sdim, int vdim>
void
VecField<real_t, sdim, vdim>::push_back_values(std::vector<real_t>& field_values,
                                               real_t time) {
  validate_field_data(num_field_points_, field_values, vdim);
  field_values_.push_back(field_values);
  time_.push_back(time);
}

template<typename real_t, int sdim, int vdim>
template<typename InterpPolicy>
void
VecField<real_t, sdim, vdim>::interp(const std::vector<real_t>& query_points,
                                     const InterpPolicy& interpolant,
                                     std::vector<real_t>& query_values,
                                     int timestep) const {
  assert (timestep >= 0 && timestep < field_values_.size());
  interpolant.interp(field_points_, sdim,
                     field_values_[timestep], vdim,
                     query_points,
                     query_values);
}

template<typename real_t, int sdim, int vdim>
template<typename InterpPolicy>
void
VecField<real_t, sdim, vdim>::interp(const std::vector<real_t>& query_points,
                                     const real_t query_time,
                                     const InterpPolicy& interpolant,
                                     std::vector<real_t>& query_values
                                     ) const {
  interpolant.interp(field_points_, sdim,
                     field_values_, vdim,
                     time_,
                     query_points,
                     query_time,
                     query_values);
}


template<typename real_t, int sdim, int vdim>
void
VecField<real_t, sdim, vdim>::write2file(const char* file_name, int timestep) {
  assert (timestep < field_values_.size());
  // FIXME: is correct only for the regular grid
  float dN = pow(num_field_points_, 1./(sdim));
  // int cll_cnt = pow(dN-1, sdim);

  // open file for writing.
  std::stringstream vtufname;
  vtufname << file_name << std::setfill('0') << std::setw(6) << timestep << ".vts";
  std::ofstream vtsfile;
  vtsfile.open(vtufname.str().c_str());
  if (vtsfile.fail()) return;

  // write to file.
  size_t data_size = 0;
  vtsfile << "<?xml version=\"1.0\"?>\n";
  vtsfile << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" ";
  if (is_little_endian())
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
  int32_t points_block_size = field_points_.size()*sizeof(real_t);
  data_size += sizeof(uint32_t)+ points_block_size;
  //**************************************************************************
  vtsfile <<"      <PointData>\n";
  vtsfile <<"        <DataArray type=\"Float" << sizeof(real_t)*8 << "\" "
          <<"Name=\"values\" "
          <<"NumberOfComponents=\"" << vdim << "\" "
          <<"format=\"appended\" "
          <<"offset=\"" << data_size << "\" />\n";
  vtsfile <<"      </PointData>\n";
  int32_t vls_block_size = field_values_[timestep].size()*sizeof(real_t);
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
  for (int i = 0; i < num_field_points_; i++) {
    for(int j = 0; j < sdim; j++)
      vtsfile.write(reinterpret_cast<char*>(&field_points_[i+j*num_field_points_]),
                    sizeof(real_t));
  }

  // points data
  // see the data memmory layout in vector field class
  vtsfile.write(reinterpret_cast<char*>(&vls_block_size),
                sizeof(vls_block_size));
  std::vector<real_t>& fvalues = field_values_[timestep];
  for (int i = 0; i < num_field_points_; i++) {
    for(int j = 0; j < vdim; j++) {
      vtsfile.write(reinterpret_cast<char*>(&fvalues[i+j*num_field_points_]),
                    sizeof(real_t));
    }
  }
  vtsfile <<"\n";
  vtsfile <<"  </AppendedData>\n";
  //**************************************************************************
  vtsfile <<"</VTKFile>\n";
  vtsfile.close();
}

template<typename real_t, int sdim, int vdim>
void
VecField<real_t, sdim, vdim>::save(const char* file_name) {
  int num_timesteps = field_values_.size();
  for (int i = 0 ; i < num_timesteps; i++)
    this->write2file(file_name, i);

  // save the pvd file
  // TODO
}

}  // namespace tbslas
