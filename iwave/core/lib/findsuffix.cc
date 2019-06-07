#include "findsuffix.hh"

namespace TSOpt {
  std::string findsuf(std::string fname) {
    size_t pos = fname.find_last_of(".");
    if (pos==std::string::npos || pos >= fname.size()-1) {
      RVL::RVLException e;
      e<<"Error: findsuf\n";
      e<<"  filename "<<fname<<" has no suffix\n";
      throw e;
    }
    size_t net = fname.size()-pos-1;
    return fname.substr(pos+1,net);
  }
}
