#include <set>
#include <string>
#include <iostream>
#include <boost/python.hpp>
#include "thulip/cpp/thermo_loader.h"
#include "thulip/cpp/thermo.h"

using namespace std;


Thermo create_thulip(string name, boost::python::list& ns){
  ThermoModelLoader loader(name);
  set<string> comps;
  for(int i = 0; i < len(ns); ++i)
    comps.insert(boost::python::extract<string>(ns[i]));

  return Thermo(loader.init_with_comps(comps));
}

boost::python::list thermo_get(Thermo& t, const char* key){
  boost::python::list l;
  int nel = t.number_of_elements(key);
  double* ptr = t.get_ptr(key);
  double unitconv = t.unitconvfac(key);
  for(int i = 0; i < nel; ++i)
     l.append(ptr[i] / unitconv);

  return l;
}

void set_by_pylist(Thermo& t, const char* key, boost::python::list& l){
  int nel = t.number_of_elements(key);
  if(boost::python::len(l) > nel)
    throw ThermoException("set_by_pylist: list longer than number of elements");
  double* ptr = t.get_ptr(key);
  double unitconv = t.unitconvfac(key);
  for(int i = 0; i < boost::python::len(l); ++i)
    ptr[i] = boost::python::extract<double>(l[i] * unitconv);
}

BOOST_PYTHON_MODULE(thulip_ext)
{
   using namespace boost::python;
   def("create_thulip", create_thulip);
   def("thermo_get", thermo_get);
   def("set_by_pylist", set_by_pylist);

   class_<Thermo>("Thermo", init<const Thermo&>())
     .def("unitconvfac", &Thermo::unitconvfac)
     .def("number_of_elements", &Thermo::number_of_elements)
     .def("number_of_dims", &Thermo::number_of_dims)
     .def("timestamp", &Thermo::timestamp)
     .def("name", &Thermo::name)
     .def("help", &Thermo::help)
     .def("keys", &Thermo::keys)
     .def("memsize", &Thermo::memsize)
     .def("calc", &Thermo::calc)
     .def("calc_nothrow", &Thermo::calc_nothrow)
     .def("has_key", &Thermo::has_key)
     .def("thermo_element", &Thermo::thermo_element)
     .def("dim_info", &Thermo::dim_info)
     .def("unit_info", &Thermo::unit_info)
     .def("get_scalar", &Thermo::get_scalar)
     .def("set_scalar", &Thermo::set_scalar)
     .def("get", &Thermo::get)
     .def("set", &Thermo::set)
     ;
}
