#pragma once

#include "../ergm_Rutil.h"
#include <string>

// Generic proxy classes for accessing an R-backed struct's list elements
// and attributes. The template parameter RStruct must have a member `SEXP R`.

template<typename RStruct>
class RAttrProxy {
public:
  explicit RAttrProxy(RStruct* p): p_(p) {}
  SEXP operator[](const char* name) const { return getAttrib(p_->R, Rf_install(name)); }
  SEXP operator[](const std::string& name) const { return getAttrib(p_->R, Rf_install(name.c_str())); }
private:
  RStruct* p_;
};

template<typename RStruct>
class RListProxy {
public:
  explicit RListProxy(RStruct* p): p_(p), attr(p) {}
  SEXP operator[](const char* name) const { return getListElement(p_->R, name); }
  SEXP operator[](const std::string& name) const { return getListElement(p_->R, name.c_str()); }
  operator SEXP() const { return p_->R; }
  RAttrProxy<RStruct> attr;
private:
  RStruct* p_;
};
