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

#ifndef SRC_UTILS_METADATA_H_
#define SRC_UTILS_METADATA_H_

#include <string>
#include <map>

namespace tbslas {

template<typename _KeyType,
         typename _ValueType,
         typename _DescType>
class MetaData {
 public:
  MetaData();
  virtual ~MetaData();

  static void
  AddMetaData(_KeyType key,
              _ValueType value,
              _DescType desc
              ) {
    values[key]       = value;
    descriptions[key] = desc;
  }

  static
  void Print() {
    printf("======================================================================\n");
    printf(" METADATA\n");
    printf("======================================================================\n");
    typename std::map<_KeyType, _ValueType>::const_iterator _iter;
    _iter = values.begin();
    for (; _iter != values.end(); ++_iter) {
      _KeyType _key   = _iter->first;
      _ValueType _val = _iter->second;
      _DescType _desc = descriptions[_key];
      std::cout << std::setw(10)
                << std::left
                << _key
                << std::setw(10)
                << std::left
                << _val
                << std::setw(30)
                << std::left
                << _desc
                << std::endl;
    }
    printf("======================================================================\n");
  }

 private:
  // static members declaration
  static std::map<_KeyType, _ValueType> values;
  static std::map<_KeyType, _DescType> descriptions;
};

// static members definition
template<typename _KeyType,
         typename _ValueType,
         typename _DescType>
std::map<_KeyType, _ValueType>
MetaData<_KeyType, _ValueType, _DescType>::values;

template<typename _KeyType,
         typename _ValueType,
         typename _DescType>
std::map<_KeyType, _DescType>
MetaData<_KeyType, _ValueType, _DescType>::descriptions;

}  // namespace tbslas

#endif  // SRC_UTILS_METADATA_H_
