/*  File inst/include/cpp/ergm_R_proxy.h in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#pragma once

#include "../ergm_Rutil.h"
#include <string>

namespace ergm {
inline namespace v1 {

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
  explicit RListProxy(SEXP sexp): attr(sexp), sexp_(sexp) {}
  SEXP operator[](const char* name) const { return getListElement(sexp_, name); }
  SEXP operator[](const std::string& name) const { return getListElement(sexp_, name.c_str()); }
  operator SEXP() const { return sexp_; }
  RAttrProxy attr;
private:
  SEXP sexp_;
};

} // namespace v1
} // namespace ergm
