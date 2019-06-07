#include <su.h>
#include <segy.h>
#include <header.h>
#include <parser.h>

#include "utils.h"
#include "std_cpp_includes.hh"

typedef struct s_Rtuple{
  RPNT coor;
} Rtuple;

using RVL::RVLException;

// Extracting receiver/source positions from given sufile.
//----------------------------------------------------------------------//
vector<Rtuple> ex_pos( string sufile, bool rcv ){
//----------------------------------------------------------------------//
  try{
    
    vector<Rtuple> pos;
    
    /*
      #ifdef IWAVE_USE_MPI
      if( retrieveGlobalRank()==0 ){
      #endif
    */
    //Reading in file
    FILE *fp = NULL;
    if (!(fp=fopen(sufile.c_str(),"r"))) {
      RVLException e;
      e << "failed to open su file = "<<sufile<<"\n";
      throw e;
    }
    segy tr;
    
   //extracting positions 
    fvgettr(fp,&tr);
    Rtuple curr;

    if(rcv){
      curr.coor[0] = tr.gelev;
      curr.coor[1] = tr.gx;
      curr.coor[2] = tr.gy;
    }
    else{
      curr.coor[0] = tr.selev;
      curr.coor[1] = tr.sx;
      curr.coor[2] = tr.sy;
    }
    if(curr.coor[0]<0) curr.coor[0]*=-1;

    pos.push_back(curr);
    int s = pos.size();
    
    while(fvgettr(fp,&tr)){
      if(rcv){
	curr.coor[0] = tr.gelev;
	curr.coor[1] = tr.gx;
	curr.coor[2] = tr.gy;
      }
      else{
	curr.coor[0] = tr.selev;
	curr.coor[1] = tr.sx;
	curr.coor[2] = tr.sy;
      }
      if(curr.coor[0]<0) curr.coor[0]*=-1;
      
      if( curr.coor[0]!=pos[s-1].coor[0] ||
	  curr.coor[1]!=pos[s-1].coor[1] ||
	  curr.coor[2]!=pos[s-1].coor[2] ){
	pos.push_back(curr);
      }
      s = pos.size();
    }
    fclose(fp);
    
    /*
      #ifdef IWAVE_USE_MPI
      }
      
      int sz = pos.size();
      MPI_Bcast(&sz,1,MPI_INT,0,retrieveGlobalComm()); 
      
      //broadcasting array containing pos to other ranks
      double * arr = new double[sz*3];
      if( retrieveGlobalRank()==0 ){
      for( int i=0; i<sz; i++ ){
      for( int j=0; j<3; j++ ){
      arr[j+i*3] = pos[i].coor[j];
      }
      }
      }
      MPI_Bcast(arr,sz*3,MPI_DOUBLE,0,retrieveGlobalComm()); 
	  
      //reading into pos
      if( retrieveGlobalRank()!=0 ){
      pos.resize(sz);
      for( int i=0; i<sz; i++ ){
      for( int j=0; j<3; j++ ){
      pos[i].coor[j] = arr[j+i*3];
      }
      }
      }
      delete[] arr;
      #endif
    */      
    return pos;
  }
  catch(RVLException &e){
    e << "ERROR from ex_pos()!\n";
    throw e;
  }
}


/* 
Analytical solution for point trace of 3D point source pressure field, 
based on repro/iwavecal paper.

ASG PDE:
dp/dt = bulk[-div v + f(t) delta(x-x_s)]
dv/dt = -buoy grad p

output is p(x,t) where |x-x_s| = r
*/
//----------------------------------------------------------------------//
void exact(int nt_trc, float t0_trc,
	   int nt_src, float t0_src,
	   float dt,   float r,
	   float bulk, float buoy,
	   float * f,
	   float * p) {
//----------------------------------------------------------------------//
  // compute df
  float * df = (float *)(usermalloc_(nt_src*sizeof(float)));
  df[0]=(f[1]-f[0])/dt;
  df[nt_src-1]=(f[nt_src-1]-f[nt_src-2])/dt;
  for (int i=1;i<nt_src-1;i++)
    df[i]=(f[i+1]-f[i-1])/(2.0f*dt);

  float c = sqrt(bulk*buoy); //sound speed
  float xoc = r/c; //traveltime
  float T_src = (nt_src-1)*dt + t0_src; //final source time
  float pi = 4.0f*atan(1.0);

  //time loop
  for (int i=0;i<nt_trc;i++) {

    p[i]=0.0f;
    float tau = i*dt - t0_trc - xoc;

    // only nonzero for t-r/c in (t0_src,T_src)
    if( tau>=t0_src && tau<=T_src ) {
      int j = int((tau+t0_src)/dt);
      p[i] = df[j]/(4*pi*c*c*r);
    }

  }
  // cleanup
  userfree_(df);
}


int xargc;
char **xargv;

//----------------------------------------------------------------------//
int main(int argc, char ** argv) {
//----------------------------------------------------------------------//
  try {

    PARARRAY *pars=NULL;
    char *str=NULL;
    Value val;
    FILE *fp_src=NULL;
    FILE *fp_trc=NULL;
    float bulk, buoy;
    string src, trc;
    segy tr_src, tr_trc;
    int nt_src, nt_trc;
    float dt_src, dt_trc, dt;
    float delrt_src, delrt_trc;
    float *pbuf=NULL;
    float r, tmp;
    fpos_t pos_trc;
    
    xargc=argc; xargv=argv;
    //requestdoc(1);

    pars = ps_new();
    if (ps_createargs(pars,argc-1,&(argv[1]))) {
      RVLException e;
      e<<"ERROR: exact_3D.x from ps_creatargs \n";
      e<<"  called with args:\n";
      e<<"  argc = "<<argc-1<<"\n";
      for (int i=0;i<argc-1;i++) 
	e<<"  argv["<<i<<"] = "<<argv[i+1]<<"\n";
      throw e;
    }

   
    // fish constants out of command line
    bulk = RVL::valparse<float>(*pars,"bulk");
    buoy = RVL::valparse<float>(*pars,"buoy");

    // filenames for source and pressure trace
    src = RVL::valparse<std::string>(*pars,"source");
    trc = RVL::valparse<std::string>(*pars,"pressure");

    // open source file
    if (!(fp_src=fopen(src.c_str(),"r"))) {
      RVLException e;
      e<<"ERROR: exact_3D from fopen\n";
      e<<"  failed to open source trace file = "<<src<<"\n";
      throw e;
    }

    // extract nt, dt, data samples from src
    if (!(fvgettr(fp_src,&tr_src))) {
      RVLException e;
      e<<"ERROR: exact_3D from fgettr\n";
      e<<"  failed to read first trace from source file = "<<src<<"\n";
      throw e;
    }
    fclose(fp_src);
    
    str = (char *)usermalloc_(128*sizeof(char));

    strcpy(str,"ns");
    gethdval(&tr_src,str,&val);
    nt_src = vtoi(hdtype(str),val);

    strcpy(str,"dt");
    gethdval(&tr_src,str,&val);
    dt_src = vtof(hdtype(str),val);

    strcpy(str,"delrt");
    gethdval(&tr_src,str,&val);
    delrt_src = vtof(hdtype(str),val); //in ms

    // open output trace file for reading    
    if (!(fp_trc=fopen(trc.c_str(),"r"))) {
      RVLException e;
      e<<"ERROR: exact_3D from fopen\n";
      e<<"  failed to open pressure trace file = "<<trc<<"\n";
      throw e;
    }

    // extract nt, dt, data samples from trc
    if (!(fvgettr(fp_trc,&tr_trc))) {
      RVLException e;
      e<<"ERROR: exact_3D from fgettr\n";
      e<<"  failed to read first trace from pressure file = "<<trc<<"\n";
      throw e;
    }
    fclose(fp_trc);

    strcpy(str,"ns");
    gethdval(&tr_trc,str,&val);
    nt_trc = vtoi(hdtype(str),val);

    strcpy(str,"dt");
    gethdval(&tr_trc,str,&val);
    dt_trc = vtof(hdtype(str),val);

    strcpy(str,"delrt");
    gethdval(&tr_trc,str,&val);
    delrt_trc = vtof(hdtype(str),val); //in ms

    // sanity check
    if(dt_src!=dt_trc){
      RVLException e;
      e<<"ERROR: exact_3D from fgettr\n"
	"  dt for source and pressure files differ!\n"
	"  dt_src="<<dt_src<<", dt_trc="<<dt_trc<<"\n";
      throw e;
    }

    // convert dt to ms 
    dt = (0.001)*dt_src;

    // extracting receiver positions and source position
    vector<Rtuple> g_pos = ex_pos(trc,true);
    vector<Rtuple> s_pos = ex_pos(trc,false);

    /*
    cerr << " rcv positions: \n";
    for( int i=0; i<g_pos.size(); i++ ){
      cerr << "   g_pos["<<i<<"] = (";
      for( int d=0; d<3; d++ ){
	cerr << g_pos[i].coor[d] <<",";
      }
      cerr << "\n";
    }

    cerr << " src positions: \n";
    for( int i=0; i<s_pos.size(); i++ ){
      cerr << "   s_pos["<<i<<"] = (";
      for( int d=0; d<3; d++ ){
	cerr << s_pos[i].coor[d] <<",";
      }
      cerr << "\n";
    }
    */

    // data buffer for output pressure trace
    pbuf = (float *)usermalloc_(g_pos.size()*nt_trc*sizeof(float));

    // loop over receivers
    for( int i=0; i<g_pos.size(); i++){

      //computing r = |x_r-x_s|
      r=0.0f;
      for(int d=0; d<3; d++){
	tmp = g_pos[i].coor[d] - s_pos[0].coor[d];
	tmp *= tmp;
	r += tmp;
      }
      r = sqrt(r);

      //computing pressure trace
      exact( nt_trc, delrt_trc, 
	     nt_src, delrt_src, 
	     dt, r, bulk, buoy,
	     tr_src.data, pbuf+i*nt_trc);
    }

    // copying buffer to trace and writing out
    if (!(fp_trc=fopen(trc.c_str(),"r+"))) {
      RVLException e;
      e<<"ERROR: exact_3D from fopen\n";
      e<<"  failed to open pressure trace file = "<<trc<<" for writing\n";
      throw e;
    }
    for( int i_tr=0; i_tr<g_pos.size(); i_tr++ ){
      fgetpos(fp_trc,&pos_trc);
      fvgettr(fp_trc,&tr_trc);
      for( int i_t=0; i_t<nt_trc; i_t++ ){
	tr_trc.data[i_t] = pbuf[ i_t+i_tr*nt_trc ];
      }
      fsetpos(fp_trc,&pos_trc);
      fputtr(fp_trc,&tr_trc);
    }
    fflush(fp_trc);
    fclose(fp_trc);
    
    // clean up
    userfree_(str);
    userfree_(pbuf);

    exit(0);
  }
  catch (RVL::RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
