#ifndef __INV_SELFDOC__
#define __INV_SELFDOC__

const char * sdoc[] = {
  " ------------------- IWAVE ITERATIVE LINEAR INVERSION ------------------------",
  "                               based on",
  " ",
#include "simtitle.h"
#include "usage_selfdoc.h"
  " "
  " ------------------------ CHOOSING INVERSION PARAMETERS -----------------------",
  " For all IWAVE inversion drivers, choose parameters to invert by keyword: if ",
  " foo is the keyword for one of the active parameters of the chosen model class,",
  " invert for foo by adding foo_est=<path> to the argument list of the inversion",
  " command. Further details at the end of this self-doc.",
  " ",
  " ------------------------ CHOOSING INVERSION PARAMETERS -----------------------",
  " ",
#include "params.h"
#include "fd_params.h"
#include "pml_params.h"
#include "mpi_params.h"
#include "cg_pars.h"
#include "iwopt_doc.h"
NULL };

#endif
