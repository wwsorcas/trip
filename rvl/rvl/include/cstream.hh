#include "std_cpp_includes.hh"

namespace RVL {
  
  class CStream {
  private:
    FILE * fp;
  public:
    CStream(FILE * _fp): fp(_fp) {}
    CStream(CStream const & cstr): fp(cstr.fp) {}
    ~CStream() {}

    void flush() { fflush(fp); }
    CStream &  operator<< ( string str ) {
      fprintf(fp,"%s",str.c_str());
      return * this;
    }
    CStream &  operator<< ( const char* str ) { 
      fprintf(fp,"%s",str);
      return * this;
    }
    CStream &  operator<< ( int i ) {
      fprintf(fp,"%d",i);
      return * this;
    }
    CStream &  operator<< ( unsigned int i ) {
      fprintf(fp, "%u", i );
      return * this;
    }
    CStream &  operator<< ( long i ) {
      fprintf(fp, "%ld", i );
      return * this;
    }
    CStream &  operator<< ( unsigned long i ) {
      fprintf(fp, "%lu", i );
      return * this;
    }
    CStream &  operator<< ( short i ) {
      fprintf(fp, "%d", i );
      return * this;
    }
    CStream &  operator<< ( unsigned short i ) {
      fprintf(fp, "%d", i );
      return * this;
    }
    CStream &  operator<< ( double d ) {
      fprintf(fp, "%g", d );
      return * this;
    }
    CStream &  operator<< ( float d ) {
      fprintf(fp, "%g", d );
      return * this;
    }
    template<class T>
    CStream &  operator<< ( complex<T> d ) {
      fprintf(fp, "(%g,%g)", d.real(),d.imag() );
      return * this;
    }
    CStream &  operator<< ( char c ) {
      fprintf(fp,"%c",c);
      return * this;
    }
  };
}
