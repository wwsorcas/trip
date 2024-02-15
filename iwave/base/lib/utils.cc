/* 
   utils.c
   Igor Terentyev.
*/
/*============================================================================*/

#include "utils.h"

#ifndef sun
#include <mm_malloc.h>
#endif

/*----------------------------------------------------------------------------*/

static int m_rk = 0;
static int m_sz = 1;
static MPI_Comm m_cm = 0;
static FILE *m_stream = 0;

/* added 31.01.09 WWS
 */
static int m_grk = 0;
static int m_gsz = 1;
static MPI_Comm m_gcm = 0;
static int m_ts = 1;

/* added 21.11.10 WWS
 */
static int m_g = 0;
static int m_ng = 1;

/* added 17.01.11 WWS
 */
static MPI_Comm m_rcm = 0;

/*----------------------------------------------------------------------------*/

void* usermalloc_(size_t size)
{
#ifdef sun
  return malloc(size);
#else
  return _mm_malloc(size, USERMALLOC_ALIGN_BYTES);
#endif
}
/*----------------------------------------------------------------------------*/

void userfree_(void *ptr)
{
#ifdef sun
  free(ptr); 
#else
  _mm_free(ptr);
#endif
}
/*----------------------------------------------------------------------------*/

int* IASN(IPNT l, const IPNT r)
{
  if ( l == r ) return l;
  return (int*)memcpy((void*)l, (const void*)r, IARR_MAX_NDIM * sizeof(int));
}
/*----------------------------------------------------------------------------*/

ireal* RASN(RPNT l, const RPNT r)
{
  if ( l == r ) return l;
  return (ireal*)memcpy((void*)l, (const void*)r, RARR_MAX_NDIM * sizeof(ireal));
}
/*----------------------------------------------------------------------------*/

int gen_3n1(int ndim, int *n)
{
  int data[3] = {2, 8, 26};
    
  if ( (unsigned int)(--ndim) >= (unsigned int)3 ) return E_BADINPUT;
    
  *n = data[ndim];
    
  return 0;
}
/*----------------------------------------------------------------------------*/

int gen_i2pnt(int ndim, int i, IPNT p)
{
  int idim, k, imax;
    
  if ( gen_3n1(ndim, &imax) ) return E_BADINPUT;
  if ( (unsigned int)i > (unsigned int)imax ) return E_BADINPUT;
    
  for ( idim = 0; idim < ndim; ++idim )
    {
      k = i % 3;
      i /= 3;
      if      ( k == 0 ) p[idim] = -1;
      else if ( k == 2 ) p[idim] = 0;
      else if ( k == 1 ) p[idim] = 1;
      else return E_INTERNAL;
    }

  return 0;
}
/*----------------------------------------------------------------------------*/

int gen_pnt2i(int ndim, const IPNT p, int *i)
{
  int idim, k;
    
  if ( (unsigned int)(--ndim) >= (unsigned int)3 ) return E_BADINPUT;

  *i = 0;
  for ( idim = ndim; idim >= 0; --idim )
    {
      if      ( p[idim] == -1 ) k = 0;
      else if ( p[idim] ==  0 ) k = 2; /* code 0 as 2, so [0,...,0] is last */
      else if ( p[idim] ==  1 ) k = 1;
      else return E_BADINPUT;
      *i = 3 * (*i) + k;
    }
    
  return 0;
}
/*----------------------------------------------------------------------------*/

void storeRank(int rk)
{
  m_rk = rk;
}
/*----------------------------------------------------------------------------*/

int retrieveRank()
{
  return m_rk;
}
/*----------------------------------------------------------------------------*/

void storeSize(int sz)
{
  m_sz = sz;
}
/*----------------------------------------------------------------------------*/

int retrieveSize()
{
  return m_sz;
}
/*----------------------------------------------------------------------------*/

void storeComm(MPI_Comm cm)
{
  m_cm = cm;
}
/*----------------------------------------------------------------------------*/

MPI_Comm retrieveComm()
{
  return m_cm;
}
/*----------------------------------------------------------------------------*/

void storeOutstream(FILE *stream)
{
  m_stream = stream;
}
/*----------------------------------------------------------------------------*/

FILE* retrieveOutstream()
{
  return m_stream;
}
/* added 31.01.09 WWS */
/*----------------------------------------------------------------------------*/

void storeGlobalRank(int rk)
{
  m_grk = rk;
}
/*----------------------------------------------------------------------------*/

int retrieveGlobalRank()
{
  return m_grk;
}
/*----------------------------------------------------------------------------*/

void storeGlobalSize(int sz)
{
  m_gsz = sz;
}
/*----------------------------------------------------------------------------*/

int retrieveGlobalSize()
{
  return m_gsz;
}
/*----------------------------------------------------------------------------*/

void storeGlobalComm(MPI_Comm cm)
{
  m_gcm = cm;
}
/*----------------------------------------------------------------------------*/

MPI_Comm retrieveGlobalComm()
{
  return m_gcm;
}
/*----------------------------------------------------------------------------*/

void storeThreadSupp(int ts)
{
  m_ts = ts;
}
/*----------------------------------------------------------------------------*/

int retrieveThreadSupp()
{
  return m_ts;
}

/* added 21.11.10 WWS */
/*----------------------------------------------------------------------------*/
void storeGroupID(int g) { m_g=g; }

/*----------------------------------------------------------------------------*/
int retrieveGroupID() { return m_g; }

/*----------------------------------------------------------------------------*/
void storeNumGroups(int ng) { m_ng = ng; }

/*----------------------------------------------------------------------------*/
int retrieveNumGroups() { return m_ng; }

/* added 17.01.11 WWS */
/*----------------------------------------------------------------------------*/
void storeRemComm(MPI_Comm cm) { m_rcm = cm; }

/*----------------------------------------------------------------------------*/
MPI_Comm retrieveRemComm() { return m_rcm; }

/*----------------------------------------------------------------------------*/
int getMachineEndianness()
{
  short int word = 0x0001;
  char *byte = (char *) &word;
  return (byte[0] ? SEAMX_LITTLE_ENDIAN : SEAMX_BIG_ENDIAN);
}
/*----------------------------------------------------------------------------*/

void swapBytes(unsigned char *arr, int arrsize, int atomsize)
{
  int i, j, k, atomsize1;
  unsigned char tmp;
	
  atomsize1 = atomsize - 1;
	
  for ( i = 0; i < arrsize; ++i )
    {
      j = 0;
      k = atomsize1;
      while (j < k)
	{
	  tmp = arr[j];
	  arr[j++] = arr[k];
	  arr[k--] = tmp;
	}
      arr += atomsize;
    }
}

/*----------------------------------------------------------------------------*/

int getBoundaryFlags(IPNT lbc, IPNT rbc, int dim) {

  IPNT cdims;
  IPNT crank;
  
  IASN(cdims, IPNT_1); /* default grid size */
  IASN(crank, IPNT_0); /* default cartisian ranks */
  
#ifdef IWAVE_USE_MPI
  IPNT cpers;          /* workspace for periodic wrap info  - currently not used */
  
  MPI_Comm cm=retrieveComm();
  
  if ( MPI_Cart_get(cm, RARR_MAX_NDIM, cdims, cpers, crank)
       != MPI_SUCCESS )  return 1;
#endif
  
  /* set boundary flags */
  IASN(lbc, IPNT_0); /* default left bc flag */
  IASN(rbc, IPNT_0); /* default right bc flag */
  
  for (int i=0;i<dim;i++) {
    if (crank[i]==0) lbc[i]=1;
    if (crank[i]==cdims[i]-1) rbc[i]=1;
  }

  return 0;
}

/*----------------------------------------------------------------------------*/

float cosfun(float t) {
  return 0.5 + cos(M_PI * t)/2.0f;
}

// mutefun = 0, t \le 0; mutefun = 1, t \ge 1
float mutefun(float t) {
  t = iwave_min(1.0f,iwave_max(t,0.0f));
  return 3*t*t-2*t*t*t;
}

// returns true if a < x < b, sets fac=hermite cubic,
// else returns false and fac=1
int taperfun(int lr, float a, float b, float x, float & fac) {
  if (x<=a) {
    if (lr==0) {
      fac*=0.0f;
      return 0;
    }
    else {
      fac*=1.0f;
      return 1;
    }
  }
  else if (x>=b) {
    if (lr==0) {
      fac*=1.0f;
      return 1;
    }
    else {
      fac*=0.0f;
      return 0;
    }
  }
  else {
    if (lr==0) fac*=mutefun((x-a)/(b-a));
    else fac*=mutefun((b-x)/(b-a));
    return 0;
  }
}

