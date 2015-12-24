// *************************************************************************
// Copyright (C) 2015 by Arash Bakhtiari
// You may not use this file except in compliance with the License.
// You obtain a copy of the License in the LICENSE file.

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// *************************************************************************

#ifndef SRC_UTILS_REPORTER_H_
#define SRC_UTILS_REPORTER_H_

#include <string>
#include <vector>

namespace tbslas {

enum Rep_Data_Type {
  REP_INT,
  REP_DOUBLE,
};

template<typename _Real_Type>
class Reporter {
 public:
  Reporter();
  virtual ~Reporter();

  static void
  AddData(std::string header,
          _Real_Type value,
          tbslas::Rep_Data_Type type = tbslas::REP_DOUBLE,
          int alignment = 10) {
    header_vec.push_back(header);
    values_vec.push_back(value);
    aligns_vec.push_back(alignment);
    types_vec.push_back(type);
  }

  static
  void Report() {
    printf("======================================================================\n");
    printf(" REPORTER\n");
    printf("======================================================================\n");
    // PRINT HEADER
    std::cout << "#TBSLAS-HEADER: ";
    for (int i = 0; i < header_vec.size(); ++i) {
      std::cout << std::right
                << std::setw(aligns_vec[i])
                << header_vec[i];
    }
    std::cout << std::endl;
    // REPORT VALUES
    std::cout << "#TBSLAS-RESULT: ";
    for (int i = 0; i < header_vec.size(); ++i) {
      switch(types_vec[i]) {
        case tbslas::REP_DOUBLE:
          std::cout << std::right
                    << std::setw(aligns_vec[i])
                    << std::setprecision(2)
                    << std::scientific
                    << values_vec[i];
          break;
        case tbslas::REP_INT:
          std::cout << std::right
                    << std::setw(aligns_vec[i])
                    // << std::setprecision(2)
                    // << std::scientific
                    << static_cast<int>(values_vec[i]);
          break;
      }
    }
    std::cout << std::endl;
    printf("======================================================================\n");
  }

 private:
  // static members declaration
  static std::vector<int> aligns_vec;
  static std::vector<std::string> header_vec;
  static std::vector<_Real_Type> values_vec;
  static std::vector<tbslas::Rep_Data_Type> types_vec;
};

// static members definition
template<typename _Real_Type>
std::vector<std::string>
Reporter<_Real_Type>::header_vec;

template<typename _Real_Type>
std::vector<_Real_Type>
Reporter<_Real_Type>::values_vec;

template<typename _Real_Type>
std::vector<int>
Reporter<_Real_Type>::aligns_vec;

template<typename _Real_Type>
std::vector<tbslas::Rep_Data_Type>
Reporter<_Real_Type>::types_vec;

}  // namespace tbslas

#endif  // SRC_UTILS_REPORTER_H_
