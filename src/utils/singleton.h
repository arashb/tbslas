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

#ifndef SRC_UTILS_SINGLETON_H_
#define SRC_UTILS_SINGLETON_H_

#include <cstddef>
#include <cassert>

template<class T>
class Singleton
{
public:
    static T* Instance() {
        if (!m_pInstance)
            m_pInstance = new T;
        assert(m_pInstance != NULL);
        return m_pInstance;
    }
protected:
    Singleton();
    ~Singleton();
private:
    Singleton(Singleton const&);
    Singleton& operator=(Singleton const&);
    static T* m_pInstance;
};

template<class T> T* Singleton<T>::m_pInstance = NULL;

#endif  // SRC_UTILS_SINGLETON_H_
