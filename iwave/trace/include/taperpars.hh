#include "findit.hh"
#include "except.hh"
#include "write.hh"
#include "product.hh"

#define NCORNER 4

typedef struct s_TAPERDATA {
  std::string key;
  float t[NCORNER];
} TAPERDATA;

using RVL::Writeable;
using RVL::ROProduct;

class TaperPars: public RVL::Writeable,
	     public RVL::ROProduct<TAPERDATA> {

private:

  std::vector<TAPERDATA *> tdata;

public:

  TaperPars(std::string arg); 
  TaperPars(TaperPars const & t);
  ~TaperPars();
  size_t getSize() const { return tdata.size(); }
  TAPERDATA const & operator[](size_t i) const;
  ostream & write(ostream & str) const;
};

