#include "thulip/cpp/thermo_loader.h"
#include "thulip/cpp/thermo.h"
#include <set>
#include <string>
#include <iostream>

using namespace std;

Thermo ConstructThulip(string name, set<string> components){
  ThermoModelLoader loader(name);
  return Thermo(loader.init_with_comps(components));
}

int main(){

  set<string> comps;
  comps.insert("methane");
  comps.insert("ethane");
  
  Thermo heh = ConstructThulip("pr_c4", comps);

  cout << "heh.."<<heh.calc_nothrow()<<"\n";
  
  return 0;
}
