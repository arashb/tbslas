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

#ifndef COMMON_H
#define COMMON_H

#include <vector>

inline bool is_little_endian() {
  int n = 1;
  return (*reinterpret_cast<char*>(&n) == 1);
}

#endif /* COMMON_H */
