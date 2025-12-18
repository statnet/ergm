#pragma once

#include "../ergm_Rutil.h"
#include <string>

namespace ergm {

// Generic proxy classes for accessing a SEXP's list elements and attributes.
// These proxies store a pointer to the SEXP directly.

class RAttrProxy {
public:
  explicit RAttrProxy(SEXP sexp): sexp_(sexp) {}
  SEXP operator[](const char* name) const { return getAttrib(sexp_, Rf_install(name)); }
  SEXP operator[](const std::string& name) const { return getAttrib(sexp_, Rf_install(name.c_str())); }
private:
  SEXP sexp_;
};

class RListProxy {
public:
  explicit RListProxy(SEXP sexp): sexp_(sexp), attr(sexp) {}
  SEXP operator[](const char* name) const { return getListElement(sexp_, name); }
  SEXP operator[](const std::string& name) const { return getListElement(sexp_, name.c_str()); }
  operator SEXP() const { return sexp_; }
  RAttrProxy attr;
private:
  SEXP sexp_;
};

} // namespace ergm
