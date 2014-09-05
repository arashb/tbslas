// Created 20-May-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_RUNTIME_ERROR
#define LIKELY_RUNTIME_ERROR

#include <stdexcept>
#include <string>

namespace likely {
class RuntimeError : public std::runtime_error {
 public:
    explicit RuntimeError(std::string const &reason);
    virtual ~RuntimeError() throw ();
 private:
  }; // RuntimeError
  inline RuntimeError::RuntimeError(std::string const &reason)
      : std::runtime_error(reason) { }

  inline RuntimeError::~RuntimeError() throw () { }

}  // namespace likely

#endif // LIKELY_RUNTIME_ERROR
