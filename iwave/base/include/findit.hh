#include "std_cpp_includes.hh"

/** find a next substring of string arg, terminated by character sep,
    convert to a thing of type T if possible. pos is beginning of
    search zone on call, next occurrence of sep orelse length of arg
    on return. If at end, returns empty string. Integer flag ierr is =
    0 for successful conversion to T, else 1.
*/
template<typename T>
std::string findit(std::string arg, char sep, size_t & pos, T & thingie, int & ierr) {
  std::string ret = "";
  // bail if at end already
  if (pos >= arg.size()) { return ret; }
  // skip over sep if currently there
  if (arg[pos]==sep) pos++;
  // save current initial position
  size_t posc=pos;
  // look for next sep, or end
  while (pos < arg.size() && arg[pos] != sep) pos++;
  // extract substring
  ret = arg.substr(posc,pos-posc);
  // convert to type T if possible
  stringstream alpha;
  if (ret != "") alpha << ret;
  alpha >> thingie;
  // success flag
  ierr=0;
  if ( (alpha.rdstate() & std::ifstream::failbit ) != 0 ) ierr=1;
  // return value is substring
  return ret;
}


