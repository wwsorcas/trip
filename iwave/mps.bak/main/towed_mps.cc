// towed_mps.cc
// Author: Mario J. Bencomo
// last modified: 04/02/18
// Based on towed_array.cc

#include <MPS_Space.hh>

using TSOpt::Rtuple;

const char * sdoc[] = { 
  "Usage: towed_mps.x data= src= towed_src=",
  "",
  "Purpose: create the effect of a towed multipole source by replicating given",
  "source traces with source coordinates specified by data file.",
  "",
  "Required parameters:",
  "  data           = filename containing source coordinates",
  "  src            = source traces to be replicated",
  "  towed          = output file",
  "",
  NULL};


int main(int argc, char ** argv) {
  
  try{
      
    // ARGUMENTS 
    string data;     // trace data file name
    string src;      // source file name 
    string towed;    // towed source file name
    
    // INTERNAL VARIABLES 
    PARARRAY * par;  // param array 
    FILE * fp_src;   // source input file pointer
    FILE * fp_towed; // towed output file pointer
    segy tr;         // trace workspace 
    
    xargc=argc; xargv=argv;
    requestdoc(1);
    
    // extract input parameters
    par=ps_new();
    if ( ps_createargs(par, argc - 1, argv + 1) ) {
      printf("Error parsing input data. ABORT.\n");
      exit(1);
    }
    
    data  = RVL::valparse<string>(*par,"data");
    src   = RVL::valparse<string>(*par,"source");
    towed = RVL::valparse<string>(*par,"towed");

    // open data files
    if (!(fp_src=fopen(src.c_str(),"r"))) {
      RVLException e;
      e << "Failed to open source input file="<<src<<".\n";
      throw e;
    }
    if (!(fp_towed=fopen(towed.c_str(),"w"))) {
      RVLException e;
      e << "Failed to open towed output file="<<towed<<".\n";
      throw e;
    }

    // extracting source locations from data file
    vector<Rtuple> spos = TSOpt::ex_pos_serial(data,false);
    int NSRC    = spos.size(); //number of sources
    //number of traces in source file = MPS dimension
    int MPS_dim = TSOpt::get_ntr_serial(src); 

    for( int i_src=0; i_src<NSRC; i_src++ ){
      fseek(fp_src,0L,SEEK_SET);      

      for( int i_b=0; i_b<MPS_dim; i_b++ ){
	fgettr(fp_src,&tr);
	
	tr.selev = -spos[i_src].coor[0];
	tr.sx = spos[i_src].coor[1];
	tr.sy = spos[i_src].coor[2];
	fputtr(fp_towed,&tr);
      }
    }

    ps_delete(&par);
    fclose(fp_src);
    fclose(fp_towed);
  }
  catch(RVLException &e){
    e << "Exiting with error from towed_mps.x!\n";
    e.write(cerr);
    exit(1);
  }
}
