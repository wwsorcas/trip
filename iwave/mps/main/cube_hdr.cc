// cube_hdr.cc
// Author: Mario J. Bencomo
// last modified: 02/20/17

#include <MPS_Space.hh>

using TSOpt::Rtuple;

const char * sdoc[] = { 
  "Usage: cube_hdr.x infile= outfile= x0_o= x1_o= x0_e= x1_e=",
  "",
  "Purpose: create a SEGY header file with receivers placed in a rectangular prism.",
  "",
  "Required parameters:",
  "  x#_o    = starting point in x#-coordinate, #=0,1,2",
  "  x#_e    = end point in x#-coordinate, #=0,1,2",
  "  x#_d    = distance between receivers in x#-coordinate",
  "  infile  = input filename containing trace info (ns,dt,delrt source positions) for output",
  "  outfile = output filename",
  "",
  "Optional parameters:",
  "  full = flag for type of data coverage",
  "         1=full cube (default value is 1)",
  "         0=sparse cube",
  "",
  NULL};


int main(int argc, char ** argv) {
  
  try{
    /* ARGUMENTS */
    string outfile;
    string infile;
    RPNT x_o, x_e, x_d, x_m;
    IPNT x_n;
    int full;

    /* INTERNAL VARIABLES */
    PARARRAY * par;  
    FILE * fp_in;    
    FILE * fp_out;   
    segy tr;         
    
    xargc=argc; xargv=argv;
    requestdoc(1);
    
    // extract input parameters
    par=ps_new();
    if ( ps_createargs(par, argc - 1, argv + 1) ) {
      printf("Error parsing input data. ABORT.\n");
      exit(1);
    }
    
    infile  = RVL::valparse<string>(*par,"infile");
    outfile = RVL::valparse<string>(*par,"outfile");
    full    = RVL::valparse<int>(*par,"full",1);
    
    x_o[0]  = RVL::valparse<float>(*par,"x0_o",0.0);
    x_o[1]  = RVL::valparse<float>(*par,"x1_o",0.0);
    x_o[2]  = RVL::valparse<float>(*par,"x2_o",0.0);

    x_e[0]  = RVL::valparse<float>(*par,"x0_e",0.0);
    x_e[1]  = RVL::valparse<float>(*par,"x1_e",0.0);
    x_e[2]  = RVL::valparse<float>(*par,"x2_e",0.0);

    x_d[0]  = RVL::valparse<float>(*par,"x0_d",1.0);
    x_d[1]  = RVL::valparse<float>(*par,"x1_d",1.0);
    x_d[2]  = RVL::valparse<float>(*par,"x2_d",1.0);
    
    // number receivers in each direction
    x_n[0]  = int((x_e[0]-x_o[0])/x_d[0])+1;
    x_n[1]  = int((x_e[1]-x_o[1])/x_d[1])+1;
    x_n[2]  = int((x_e[2]-x_o[2])/x_d[2])+1;

    // midpoints in each direction
    x_m[0]  = (x_e[0]+x_o[0])*0.5;
    x_m[1]  = (x_e[1]+x_o[1])*0.5;
    x_m[2]  = (x_e[2]+x_o[2])*0.5;

    cerr << " x_o =("<<x_o[0]<<","<<x_o[1]<<","<<x_o[2]<<")\n"
	 << " x_e =("<<x_e[0]<<","<<x_e[1]<<","<<x_e[2]<<")\n"
	 << " x_d =("<<x_d[0]<<","<<x_d[1]<<","<<x_d[2]<<")\n"
	 << " x_n =("<<x_n[0]<<","<<x_n[1]<<","<<x_n[2]<<")\n"
	 << " x_m =("<<x_m[0]<<","<<x_m[1]<<","<<x_m[2]<<")\n";

    // open data files
    if (!(fp_in=fopen(infile.c_str(),"r"))) {
      RVLException e;
      e << "Failed to open input file="<<infile<<".\n";
      throw e;
    }
    if (!(fp_out=fopen(outfile.c_str(),"w"))) {
      RVLException e;
      e << "Failed to open output file="<<outfile<<".\n";
      throw e;
    }

    // extracting and zeroing sample trace
    fgettr(fp_in,&tr);
    for( int i_t=0; i_t<tr.ns; i_t++ ){
      tr.data[i_t] = 0.0;
    }

    // extracting source locations from data file
    vector<Rtuple> spos = TSOpt::ex_pos(infile,false);
    int NSRC = spos.size();
    tr.tracl = 1;

    // looping over sources
    for( int i_src=0; i_src<NSRC; i_src++ ){

      tr.selev = -spos[i_src].coor[0];
      tr.sx    = spos[i_src].coor[1];
      tr.sy    = spos[i_src].coor[2];


      ////////////////////////////
      // faces x0=x_0[0],x_e[0] //
      ////////////////////////////
      RPNT gx0; 
      gx0[0] = x_o[0];
      gx0[1] = x_e[0];

      for( int i=0; i<2; i++){
	if( i==1 && x_n[0]<=1 ) break;

	tr.gelev = gx0[i];
	
	//full face
	if(full){
	  for( int i1=0; i1<x_n[1]; i1++ ){
	    tr.gx = x_o[1] + x_d[1]*i1;
	    for( int i2=0; i2<x_n[2]; i2++ ){
	      tr.gy = x_o[2] + x_d[2]*i2;
	      fputtr(fp_out,&tr);
	      tr.tracl++;
	    }
	  }
	}
	//sparse face
	else{
	  if(x_n[2]>1){
	    tr.gx = x_m[1];
	    for( int i2=0; i2<x_n[2]; i2++ ){
	      tr.gy = x_o[2] + x_d[2]*i2;
	      fputtr(fp_out,&tr);
	      tr.tracl++;
	    }
	  }
	  if(x_n[1]>1){
	    tr.gy = x_m[2];
	    for( int i1=0; i1<x_n[1]; i1++ ){
	      tr.gx = x_o[1] + x_d[1]*i1;
	      fputtr(fp_out,&tr);
	      tr.tracl++;
	    }
	  }
	}
      }
      
      ///////////////////////////
      // face x1=x_o[1],x_e[1] //
      ///////////////////////////
      RPNT gx1; 
      gx1[0] = x_o[1];
      gx1[1] = x_e[1];

      for( int i=0; i<2; i++){
	if( i==1 && x_n[1]<=1 ) break;

	tr.gx = gx1[i];
	
	//full face
	if(full){
	  for( int i0=0; i0<x_n[0]; i0++ ){
	    tr.gelev = x_o[0] + x_d[0]*i0;
	    for( int i2=0; i2<x_n[2]; i2++ ){
	      tr.gy = x_o[2] + x_d[2]*i2;
	      fputtr(fp_out,&tr);
	      tr.tracl++;
	    }
	  }
	}
	//sparse face
	else{
	  if(x_n[2]>1){
	    tr.gelev = x_m[0];
	    for( int i2=0; i2<x_n[2]; i2++ ){
	      tr.gy = x_o[2] + x_d[2]*i2;
	      fputtr(fp_out,&tr);
	      tr.tracl++;
	    }
	  }
	  if(x_n[0]>1){
	    tr.gy = x_m[2];
	    for( int i0=0; i0<x_n[0]; i0++ ){
	      tr.gelev = x_o[0] + x_d[0]*i0;
	      fputtr(fp_out,&tr);
	      tr.tracl++;
	    }
	  }
	}
      }

      ///////////////////////////
      // face x2=x_o[2],x_e[2] //
      ///////////////////////////
      RPNT gx2; 
      gx2[0] = x_o[2];
      gx2[1] = x_e[2];

      for( int i=0; i<2; i++){
	if( i==1 && x_n[2]<=1 ) break;

	tr.gy = gx2[i];
	
	//full face
	if(full){
	  for( int i0=0; i0<x_n[0]; i0++ ){
	    tr.gelev = x_o[0] + x_d[0]*i0;
	    for( int i1=0; i1<x_n[1]; i1++ ){
	      tr.gx = x_o[1] + x_d[1]*i1;
	      fputtr(fp_out,&tr);
	      tr.tracl++;
	    }
	  }
	}
	//sparse face
	else{
	  if(x_n[1]>1){
	    tr.gelev = x_m[0];
	    for( int i1=0; i1<x_n[1]; i1++ ){
	      tr.gx = x_o[1] + x_d[1]*i1;
	      fputtr(fp_out,&tr);
	      tr.tracl++;
	    }
	  }
	  if(x_n[0]>1){
	    tr.gx = x_m[1];
	    for( int i0=0; i0<x_n[0]; i0++ ){
	      tr.gelev = x_o[0] + x_d[0]*i0;
	      fputtr(fp_out,&tr);
	      tr.tracl++;
	    }
	  }
	}
      }

    }
    
    fclose(fp_in);
    fclose(fp_out);
    ps_delete(&par);
    
  }
  catch(RVLException &e){
    e << "Exiting with error from cube_hdr.x!\n";
    e.write(cerr);
    exit(1);
  }
}
