// MPS_to_RHS.cc
// Author: Mario J. Bencomo
// last modified: 11/10/16

#include "MPS_to_RHS.hh"

//#define VERBOSE_MJB

namespace TSOpt {


  //----------------------------------------------------------------------//
  void write_RHS_SEGY( vector<MPS_stencil> stencils,
		       vector<Rtuple> s_pos,
		       string ref_file,
		       vector<string> RHS_files ){
  //----------------------------------------------------------------------//
    try{

#ifdef VERBOSE_MJB
      cerr << "Inside write_RHS_SEGY\n";
#endif

#ifdef IWAVE_USE_MPI
      if( retrieveGlobalRank()!=0 ){
	RVLException e;
	e << "write_RHS_SEGY should only be called by rank 0!\n";
      }
#endif

      //extrating prototype trace from ref_file
      FILE *fp_ref = iwave_const_fopen(ref_file.c_str(),"r",NULL,stderr); 
      segy tr;
      fvgettr(fp_ref,&tr);
      iwave_fclose(fp_ref);

      //zeroing trace
      for( int it=0; it<tr.ns; it++)
	tr.data[it]=0;
      
      ////////////////////////////
      // looping over RHS files //
      ////////////////////////////
      for( int i_rhs=0; i_rhs<RHS_files.size(); i_rhs++){
	
	//writing output
	FILE *fp_rhs = iwave_const_fopen(RHS_files[i_rhs].c_str(),"w",NULL,stderr); 
	int i_tr_rhs = 0;

	///////////////////////
	// loop over sources //
	///////////////////////
	for( int i_src=0; i_src<stencils.size(); i_src++){
	  
	  tr.selev = -s_pos[i_src].coor[0];
	  tr.sx    = s_pos[i_src].coor[1];
	  tr.sy    = s_pos[i_src].coor[2];

#ifdef VERBOSE_MJB
	  cerr <<"  inside source loop, i_src="<<i_src<<"\n"
	       <<"    selev = "<<tr.selev<<"\n"
	       <<"    sx    = "<<tr.sx<<"\n"
	       <<"    sy    = "<<tr.sy<<"\n";
#endif

	  
	  /////////////////////////////////////
	  // loop over source stencil points //
	  /////////////////////////////////////
	  for( int i_pt=0; i_pt<stencils[i_src].points.size(); i_pt++){
	    
	    tr.gelev =-stencils[i_src].points[i_pt].loc[0];
	    tr.gx    = stencils[i_src].points[i_pt].loc[1];
	    tr.gy    = stencils[i_src].points[i_pt].loc[2];

	    tr.tracl = i_tr_rhs+1;
	    tr.tracr = i_tr_rhs+1;

#ifdef VERBOSE_MJB
	  cerr <<"    inside stencil point loop, i_pt="<<i_pt<<"\n"
	       <<"      gelev = "<<tr.gelev<<"\n"
	       <<"      gx    = "<<tr.gx<<"\n"
	       <<"      gy    = "<<tr.gy<<"\n"
	       <<"      tracl = "<<tr.tracl<<"\n";
#endif	    
	    fvputtr(fp_rhs,&tr);
	    i_tr_rhs++;
	  }//i_pt loop      
	}//i_src loop

	fflush(fp_rhs);
	iwave_fclose(fp_rhs);      
      }//i_rhs loop

    }
    catch(RVLException &e){
      e << "ERROR from wwrite_RHS_SEGY\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void MPS_to_RHS_kern( int a_ord,
			Ituple d_ord,
			grid g,
			vector<Rtuple> s_pos,
			vector<MPS_stencil> stencils,
			vector<MPS_base> MPS_basis,
			string MPS_file,
			vector<string> RHS_files ){
  //----------------------------------------------------------------------//
    try{

#ifdef IWAVE_USE_MPI
      if( retrieveGlobalRank()!=0 ){
	RVLException e;
	e << "MPS_to_RHS_kern should only be called by rank 0!\n";
	throw e;
      }
#endif

#ifdef VERBOSE_MJB      
      cerr << "*************************************************\n"
	   << "  Inside MPS_to_RHS_kern\n"
	   << "*************************************************\n\n";
      cerr << "MPS file = " << MPS_file << "\n"
	   << "RHS_files = "<< RHS_files[0]<<"\n";
#endif
            
      FILE *fp_mps = iwave_const_fopen(MPS_file.c_str(),"r",NULL,stderr); 
      segy tr_mps;
      fpos_t pos_rhs;

      //////////////////////////
      // loop over RHS source //
      //////////////////////////
      for( int i_rhs=0; i_rhs<RHS_files.size(); i_rhs++){
	
#ifdef VERBOSE_MJB
	cerr << "\nInside RHS loop: i_rhs="<< i_rhs <<"\n";
	cerr << "RHS file = " <<RHS_files[i_rhs] <<"\n";
#endif

	FILE *fp_rhs = iwave_const_fopen(RHS_files[i_rhs].c_str(),"r+",NULL,stderr);
	segy tr_rhs;

	//trace index for rhs
	int i_tr_rhs=0;
	
	//reset file to beginning
	if (fseeko(fp_rhs, 0L, SEEK_SET)){
	  RVLException e;
	  e << "could not find trace 0 in file "<< RHS_files[i_rhs]<<"\n";
	  throw e;
	}

	///////////////////////
	// loop over sources //
	///////////////////////
	for( int i_src=0; i_src<s_pos.size(); i_src++){

#ifdef VERBOSE_MJB
	  cerr <<"\n  Inside MPS loop: i_src="<< i_src<<"\n";
#endif

	  //trace index for mps
	  int i_tr_mps=0;
	  for( int i=0; i<i_src; i++)
	    i_tr_mps += MPS_basis.size();

	  ///////////////////////////////
	  //looping over stencil points//
	  ///////////////////////////////
	  for( int i_pt=0; i_pt<stencils[i_src].points.size(); i_pt++){

#ifdef VERBOSE_MJB				  
	    cerr <<"\n    Inside point loop: i_pt="<<i_pt<<"\n";
#endif		  
	    //Getting current trace to be written
	    fgetpos(fp_rhs,&pos_rhs);
	    fgettr(fp_rhs,&tr_rhs);

	    ////////////////////////////////
	    //looping over multipole basis//
	    ////////////////////////////////
	    for( int i_base=0; i_base<MPS_basis.size(); i_base++){
	
#ifdef VERBOSE_MJB				  
	      cerr <<"\n      Inside multipole basis loop: basis i_base="<<i_base<<"\n";
#endif		  
	      if( MPS_basis[i_base].RHS_flags[i_rhs]==0 ) continue;
	      
	      //////////////////////////////////
	      //looping over MPS terms in base//
	      //////////////////////////////////
	      for( int i_trm=0; i_trm<MPS_basis[i_base].NT; i_trm++){
	
#ifdef VERBOSE_MJB				    
		cerr <<"\n        Inside multipole terms in bases loop: i_trm="<<i_trm<<"\n";
#endif
		//tensor prod weight for point source approx
		float w=1;
		float x_in;
		
		//weight in z-dim
		if( stencils[i_src].n[0]>1 ){
		  x_in = -stencils[i_src].points[i_pt].glo[0] + s_pos[i_src].coor[0];
		  w *= ssappx( x_in, 
			       g.axes[0].d,
			       MPS_basis[i_base].derv[i_trm].coor[0],
			       a_ord);
		  w *= stencils[i_src].points[i_pt].ext[0];
		}

		//weight in x-dim
		if( stencils[i_src].n[1]>1 ){		
		  x_in = stencils[i_src].points[i_pt].glo[1] - s_pos[i_src].coor[1];
		  w *= ssappx( x_in,
			       g.axes[1].d,
			       MPS_basis[i_base].derv[i_trm].coor[1],
			       a_ord);
		  w *= stencils[i_src].points[i_pt].ext[1];
		}

		//weight in y-dim
		if( stencils[i_src].n[2]>1 ){
		  x_in = stencils[i_src].points[i_pt].glo[2] - s_pos[i_src].coor[2];
		  w *= ssappx( x_in,
			       g.axes[2].d,
			       MPS_basis[i_base].derv[i_trm].coor[2],
			       a_ord);
		  w *= stencils[i_src].points[i_pt].ext[2];
		}

#ifdef VERBOSE_MJB		
		cerr << "        *** w="<<w<<"\n";	
		cerr << "        i_tr_mps+i_base="<<i_tr_mps+i_base<<"\n";		    
#endif		    

		//getting trace for current multipole term
		fgettra(fp_mps,&tr_mps,i_tr_mps+i_base);
		for( int it=0; it<tr_mps.ns; it++)
		  tr_rhs.data[it] += tr_mps.data[it]*w;
		    
	      }//loop i_trm
	    }//loop i_base
		
	    //writting out trace
#ifdef VERBOSE_MJB				
	    cerr <<"        writing out at trace i_tr_rhs="<<i_tr_rhs<<"\n";
#endif		
		
	    fsetpos(fp_rhs,&pos_rhs);
	    fvputtr(fp_rhs,&tr_rhs);
	    i_tr_rhs++;
	  }//loop i_pt 
	}//loop MPSs
     
	fflush(fp_rhs);
	iwave_fclose(fp_rhs);
      }//loop RHSs
     
      iwave_fclose(fp_mps);

#ifdef VERBOSE_MJB
      cerr <<"EXITING MPS_to_RHS_kern()\n";
#endif
    }
    catch(RVLException &e){
      e << "ERROR from MPS_to_RHS_kern()\n";
      throw e;
    }
  }


  //----------------------------------------------------------------------//
  void RHS_to_MPS_kern( int a_ord,
			Ituple d_ord,
			grid g,
			vector<Rtuple> s_pos,
			vector<MPS_stencil> stencils,
			vector<MPS_base> MPS_basis,
			string MPS_file,
			vector<string> RHS_files ){
  //----------------------------------------------------------------------//
    try{

#ifdef IWAVE_USE_MPI
      if( retrieveGlobalRank()!=0 ){
	RVLException e;
	e << "RHS_to_MPS_kern should only be called by rank 0!\n";
	throw e;
      }
#endif

#ifdef VERBOSE_MJB      
      cerr << "*************************************************\n"
	   << "  Inside RHS_to_MPS_kern\n"
	   << "*************************************************\n\n";
#endif
      
      FILE *fp_mps = iwave_const_fopen(MPS_file.c_str(),"r+",NULL,stderr); 
      segy tr_mps;
      fpos_t pos_mps;

      //reset file to beginning
      if (fseeko(fp_mps, 0L, SEEK_SET)){
	RVLException e;
	e << "could not find trace 0 in file "<< MPS_file<<"\n";
	throw e;
      }

      ///////////////////////
      // loop over sources //
      ///////////////////////
      for( int i_src=0; i_src<s_pos.size(); i_src++){
	
#ifdef VERBOSE_MJB
	cerr <<"\nInside source loop: i_src="<< i_src<<"\n";
#endif

	//////////////////////////
	//looping over MPS bases//
	//////////////////////////
	for( int i_base=0; i_base<MPS_basis.size(); i_base++){
	  
#ifdef VERBOSE_MJB				  
	  cerr <<"\nInside multipole basis loop: basis i_base="<<i_base<<"\n";

	  //computing mps trace index
	  int i_tr_mps = i_base + i_src*MPS_basis.size();
	  cerr <<"\ni_tr_mps="<<i_tr_mps<<"\n";
#endif

	  //Getting current trace to be written
	  fgetpos(fp_mps,&pos_mps);
	  fgettr(fp_mps,&tr_mps);
	  
	  ///////////////////////
	  //looping over points//
	  ///////////////////////
	  for( int i_pt=0; i_pt<stencils[i_src].points.size(); i_pt++){
	
#ifdef VERBOSE_MJB				    
	    cerr <<"\nInside point loop: i_pt="<<i_pt<<"\n";
#endif
	    //computing weight W = sum_{i_trm} w_{i_trm}
	    float W = 0.0;

	    //////////////////////////////////
	    //looping over MPS terms in base//
	    //////////////////////////////////
	    for( int i_trm=0; i_trm<MPS_basis[i_base].NT; i_trm++){
	
#ifdef VERBOSE_MJB				    
	      cerr <<"\nInside multipole terms in bases loop: i_trm="<<i_trm<<"\n";
#endif
	      //tensor prod weight for point source approx
	      float w=1;
	      float x_in;

	      //weight in z-dim
	      if( stencils[i_src].n[0]>1 ){
		x_in = -stencils[i_src].points[i_pt].glo[0] + s_pos[i_src].coor[0];
		w *= ssappx( x_in,
			     g.axes[0].d,
			     MPS_basis[i_base].derv[i_trm].coor[0],
			     a_ord);
		w *= stencils[i_src].points[i_pt].ext[0];
	      }

	      //weight in x-dim
	      if( stencils[i_src].n[1]>1 ){
		x_in = stencils[i_src].points[i_pt].glo[1] - s_pos[i_src].coor[1];
		w *= ssappx( x_in,
			     g.axes[1].d,
			     MPS_basis[i_base].derv[i_trm].coor[1],
			     a_ord);
		w *= stencils[i_src].points[i_pt].ext[1];
	      }

	      //weight in y-dim
	      if( stencils[i_src].n[2]>1 ){
		x_in = stencils[i_src].points[i_pt].glo[2] - s_pos[i_src].coor[2];
		w *= ssappx( x_in,
			     g.axes[2].d,
			     MPS_basis[i_base].derv[i_trm].coor[2],
			     a_ord); 
		w *= stencils[i_src].points[i_pt].ext[2];
	      }

#ifdef VERBOSE_MJB		
	      cerr << "*** w="<<w<<"\n";	
	      cerr << "i_tr_mps+i_base="<<i_tr_mps+i_base<<"\n";		    
#endif		    

	      //update weight
	      W += w;
		  
	    }//loop i_trm

	    //computing trace index for rhs traces
	    int i_tr_rhs = i_pt + i_src*( stencils[i_src].n[0]
					 *stencils[i_src].n[1]
					 *stencils[i_src].n[2] );

#ifdef VERBOSE_MJB
	    cerr <<"\ncomputed weight W="<<W<<"\n";
	    cerr <<"\ni_tr_rhs = "<<i_tr_rhs<<"\n";
#endif

	    //////////////////////////
	    //looping over RHS terms//
	    //////////////////////////
	    for( int i_rhs=0; i_rhs<RHS_files.size(); i_rhs++){
	      
#ifdef VERBOSE_MJB
	      cerr << "\nInside RHS loop: i_rhs="<< i_rhs <<"\n";
	      cerr << "\nMPS_basis["<<i_base<<"].RHS_flags["<<i_rhs<<"]="
		   << MPS_basis[i_base].RHS_flags[i_rhs] <<"\n";
#endif
	      if( MPS_basis[i_base].RHS_flags[i_rhs]==0 ){
		continue;
	      }
	      
	      FILE *fp_rhs = iwave_const_fopen(RHS_files[i_rhs].c_str(),"r",NULL,stderr);
	      segy tr_rhs;
	      
	      //getting trace for current RHS term
	      fgettra(fp_rhs,&tr_rhs,i_tr_rhs);
	      
#ifdef VERBOSE_MJB
	      cerr << "\nUpdating mps trace\n";
#endif
	      //updating over mps traces
	      for( int it=0; it<tr_rhs.ns; it++)
		tr_mps.data[it] += tr_rhs.data[it]*W;
	      
	      iwave_fclose(fp_rhs);
	      
	    }//loop i_rhs  
	  }//loop i_pt
	  
	  fsetpos(fp_mps,&pos_mps);
	  fputtr(fp_mps,&tr_mps);
	}//loop MPS bases
     
      }//loop sources
      
      fflush(fp_mps);
      iwave_fclose(fp_mps);
      
#ifdef VERBOSE_MJB
      cerr <<"EXITING RHS_to_MPS_kern()\n";
#endif
      //exit(1);
    }
    catch(RVLException &e){
      e << "ERROR from RHS_to_MPS_kern()\n";
      throw e;
    }
  }


  //////////////////////////////////////////////////////////////////////////
  // MPS_to_RHS class //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  //----------------------------------------------------------------------//
  MPS_to_RHS::MPS_to_RHS( MPS_to_RHS const &MtoR )
    : dom(MtoR.dom),
      rng_ptr(MtoR.rng_ptr),
      pars(NULL), 
      mks(MtoR.mks),
      make_RHS(MtoR.make_RHS),
      a_ord(MtoR.a_ord),
      g(MtoR.g),
      stencils(MtoR.stencils){
  //----------------------------------------------------------------------//
    try{

      if( !(rng_ptr) ){
	RVLException e;
	e << "empty range space pointer\n";
	throw e;
      }
      
      //copying pars
      pars = ps_new();
      if(ps_copy(&pars,*MtoR.pars)){
	RVLException e;
	e <<"failed to copy parameter table\n";
	throw e;
      }
    }
    catch(RVLException &e){
      e << "ERROR from MPS_to_RHS constructor\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  MPS_to_RHS::MPS_to_RHS( MPS_KEYS _mks, 
			  PARARRAY const &_pars, 
			  MPS_Space const & mps,
			  bool _make_RHS )
    : dom(mps),
      rng_ptr(),
      pars(NULL),
      mks(_mks),
      make_RHS(_make_RHS){
  //----------------------------------------------------------------------//
    try{

      //copying pars
      pars = ps_new();
      if(ps_copy(&pars,_pars)){
	RVLException e;
	e <<"failed to copy parameter table\n";
	throw e;
      }
      
      //extracting a_ord
      a_ord = valparse<int>(*pars,mks.appx_ord);

      //setting up spatial grid info
      string grid_file = valparse<string>(*pars,mks.grid_file);
      init_default_grid(&g);
      read_grid(&g,grid_file.c_str(),stderr);

      if(g.axes[0].o<0) g.axes[0].o*=-1;
      if(g.axes[0].d<0) g.axes[0].d*=-1;

      //set stencils
      set_stencils();

      //set range
      set_rng_ptr();

      //sanity checks
      sanity_checks();

    }
    catch(RVLException &e){
      e << "ERROR from MPS_to_RHS constructor!\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  MPS_to_RHS::~MPS_to_RHS(){
  //----------------------------------------------------------------------//
    ps_delete(&pars);
  }

  //----------------------------------------------------------------------//
  void MPS_to_RHS::set_stencils(){
  //----------------------------------------------------------------------//
    try{
      int Nsrc = dom.get_s_pos().size();
      stencils.resize(Nsrc);

      //////////////////////////
      // Looping over sources //
      //////////////////////////
      for( int i_src=0; i_src<Nsrc; i_src++ ){
	RPNT Gstn_o;
	RPNT g_o, g_e, g_d;

	//getting global stencil info
	for(int d=0; d<3; d++){

	  g_d[d] = g.axes[d].d;
	  g_o[d] = g.axes[d].o;
	  g_e[d] = g_o[d] + (g.axes[d].n-1)*g_d[d];
	  
	  stencils[i_src].n[d] = (d<g.dim)? 
	    ( a_ord + dom.get_d_ord().coor[d] ):(1);
	  
	  float sto = dom.get_s_pos()[i_src].coor[d] - stencils[i_src].n[d]*g_d[d]/2;	         
	  int idx = (int)((sto - g_o[d])/g_d[d]);
	  Gstn_o[d] = idx*g_d[d] + g_o[d];
	  float res = sto - Gstn_o[d];
	  
#ifdef VERBOSE_MJB
	  cerr << "idx_float = "<< (sto - g_o[d])/g_d[d] <<"\n";
	  cerr << "n-sten="<<stencils[i_src].n[d]<<"\n"
	       << "sto="<<sto<<"\n"
	       << "idx="<<idx<<"\n"
	       << "res="<<res<<"\n"
	       << "Gstn_o["<<d<<"]="<<Gstn_o[d]<<"\n";
#endif
	  //res==0 implies sto is a grid point, otherwise increment idx
	  if( sto==0 ){
	    if( fabs(res)>=1e-9 ) Gstn_o[d]+=g_d[d];
	  }
	  else{
	    if( fabs(res/sto)>=1e-9 && sto>0 ) Gstn_o[d]+=g_d[d];
	  }

#ifdef VERBOSE_MJB
	  cerr << "inside setting up global stencil info, d="<<d<<"\n";
	  cerr << "  g_o["<<d<<"]="<<g_o[d]<<"\n"
	       << "  g_d["<<d<<"]="<<g_d[d]<<"\n"
	       << "  g_e["<<d<<"]="<<g_e[d]<<"\n"
	       << "  Gstn_o["<<d<<"]="<<Gstn_o[d]<<"\n";
#endif
	  //Gstn_o[d] = s_idx[i_src].coor[d] - (int)((stencils[i_src].n[d]-1)/2);
	  //Gstn_o[d] *= g_d[d];
	  //Gstn_o[d] += g_o[d];
	}
	
	//setting stencil points
	vector<MPS_base> MPS_basis = dom.get_MPS_basis();
	MPS_point pt; 
	pt.ext[0]=1; 
	pt.ext[1]=1; 
	pt.ext[2]=1;
	
	//MOD
	//pt.Lidx.resize(MPS_basis.size());

	/////////////////////////
	// Looping over z-axis //
	/////////////////////////
	for( int iz=0; iz<stencils[i_src].n[0]; iz++ ){
	  
	  pt.glo[0] = Gstn_o[0] + iz*g_d[0];
	  pt.loc[0] = pt.glo[0];

	  //checking if pt.loc[0] is not an interior point
	  if( pt.loc[0]-g_o[0] < 1e-5 ) {
	    //boundary case
	    if( fabs(pt.loc[0]-g_o[0])< 1e-5 ){
	      pt.ext[0] = 1;
	    }
	    //exterior case
	    else{
	      pt.loc[0] = 2*g_o[0] - pt.loc[0];
	      pt.ext[0] = -1;	      
	    }
	  }
	  else 
	  if( pt.loc[0]-g_e[0] > -1e-5 ){
	    //boundary case
	    if( fabs(pt.loc[0]-g_o[0])< 1e-5 ){
	      pt.ext[0] = 1;
	    }
	    //exterior case
	    else{
	      pt.loc[0] = pt.loc[0] - 2*g_e[0];
	      pt.ext[0] = -1;
	    }
	  }
	

	  /////////////////////////
	  // Looping over x-axis //
	  /////////////////////////
	  for( int ix=0; ix<stencils[i_src].n[1]; ix++ ){
	    
	    pt.glo[1] = Gstn_o[1] + ix*g_d[1];
	    pt.loc[1] = pt.glo[1];
	    
	    //checking if pt.loc[1] is not an interior point
	    if( pt.loc[1]-g_o[1] < 1e-5 ) {
	      //boundary case
	      if( fabs(pt.loc[1]-g_o[1])< 1e-5 ){
		pt.ext[1] = 1;
	      }
	      //exterior case
	      else{
		pt.loc[1] = 2*g_o[1] - pt.loc[1];
		pt.ext[1] = -1;	      
	      }
	    }
	    else 
	    if( pt.loc[1]-g_e[1] > -1e-5 ){
	      //boundary case
	      if( fabs(pt.loc[1]-g_o[1])< 1e-5 ){
		pt.ext[1] = 1;
	      }
	      //exterior case
	      else{
		pt.loc[1] = pt.loc[1] - 2*g_e[1];
		pt.ext[1] = -1;
	      }
	    }


	    /////////////////////////
	    // Looping over y-axis //
	    /////////////////////////
	    for( int iy=0; iy<stencils[i_src].n[2]; iy++ ){

	      pt.glo[2] = Gstn_o[2] + iy*g_d[2];
	      pt.loc[2] = pt.glo[2];

	      //checking if pt.loc[2] is not an interior point
	      if( pt.loc[2]-g_o[2] < 1e-5 ) {
		//boundary case
		if( fabs(pt.loc[2]-g_o[2])< 1e-5 ){
		  pt.ext[2] = 1;
		}
		//exterior case
		else{
		  pt.loc[2] = 2*g_o[2] - pt.loc[2];
		  pt.ext[2] = -1;	      
		}
	      }
	      else 
	      if( pt.loc[2]-g_e[2] > -1e-5 ){
		//boundary case
		if( fabs(pt.loc[2]-g_o[2])< 1e-5 ){
		  pt.ext[2] = 1;
		}
		//exterior case
		else{
		  pt.loc[2] = pt.loc[2] - 2*g_e[2];
		  pt.ext[2] = -1;
		}
	      }
	      
	      
	      stencils[i_src].points.push_back(pt);
	      pt.ext[2]=1;
	    }
	    pt.ext[1]=1;
	  }
	  pt.ext[0]=1;
	}

#ifdef VERBOSE_MJB
	cerr << "--- From MPS_to_RHS::set_stencils(), printing out stencil["<<i_src<<"]!\n";
	for( int i=0; i<stencils[i_src].points.size(); i++ ){
	  cerr << "  pt.loc = (" 
	       << stencils[i_src].points[i].loc[0] << ","
	       << stencils[i_src].points[i].loc[1] << ","
	       << stencils[i_src].points[i].loc[2] << ")\n"
	       << "  pt.glo = (" 
	       << stencils[i_src].points[i].glo[0] << ","
	       << stencils[i_src].points[i].glo[1] << ","
	       << stencils[i_src].points[i].glo[2] << ")\n"
	       <<"  ext = (" 
	       << stencils[i_src].points[i].ext[0] << ","
	       << stencils[i_src].points[i].ext[1] << ","
	       << stencils[i_src].points[i].ext[2] << ")\n";
	}
#endif
	
      }//src loop
    }
    catch(RVLException &e){
      e << "ERROR from MPS_to_RHS::set_stencils()\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void MPS_to_RHS::set_rng_ptr(){
  //----------------------------------------------------------------------//
    try{

      vector<string> RHS_files;
      
      ///////////////////////////////////////////////
      // case where RHS files need to be generated //
      ///////////////////////////////////////////////
      if( make_RHS ){

	try{

	  //Extracting RHS filenames, if empty will generate filenames
	  for( int i=0; i<mks.RHS_files.size(); i++){
	    string fn = valparse<string>(*pars,mks.RHS_files[i]);
	    if(fn.compare("empty")==0){
	      RVLException e;
	      throw e;
	    }
	    RHS_files.push_back(valparse<string>(*pars,mks.RHS_files[i]));
	  }
	  
	}
	//Generating RHS filenames
	catch(RVLException &e){
	  cerr << "WARNING: from MPS_to_RHS::set_rng_ptr\n"
	       << "No RHS filenames given, will generate names: ";
	  
	  RHS_files.clear();
	  
	  for( int i=0; i<mks.RHS_files.size(); i++){
	    //string pre = "RHS_"+ to_string(i) +"_";
	    std::stringstream ss;
	    ss << "RHS_" << i << "_";
	    string pre = ss.str();
	    string tmp = dom.get_coreFilename();
	    size_t pos = tmp.find_last_of("/");
	    string path = tmp.substr(0,pos+1);
	    string file = tmp.substr(pos+1);
	    
	    string fn = path + pre + file;
	    RHS_files.push_back(fn);
	    cerr << fn << ",";
	  }
	  cerr <<"\n";
	}
	
	//generating RHS SEGY files
#ifdef IWAVE_USE_MPI
	  if( retrieveGlobalRank()==0 ){
#endif

	    write_RHS_SEGY( get_stencils(),
			    dom.get_s_pos(),
			    dom.get_coreFilename(),
			    RHS_files );

	    
#ifdef IWAVE_USE_MPI
	  }
	  MPI_Barrier(retrieveGlobalComm());
#endif
      }
      ////////////////////////////////
      // case where RHS files exist //
      ////////////////////////////////
      else{
	for( int i=0; i<mks.RHS_files.size(); i++)
	  RHS_files.push_back(valparse<string>(*pars,mks.RHS_files[i]));
      }
      
      ///////////////////////////
      // generating RHS spaces //
      ///////////////////////////
      
      //updating parameter table
      add_to_pars( *pars, mks.RHS_files, RHS_files );
      //cerr << "Inside gen_RHS_Space\n";
      //ps_printall(*pars,stderr);
      
      //generating space
      //will replace with:
      //shared_ptr<IWaveSpace> sp_ptr 
      //  = make_shared<IWaveSpace>(*pars,ic,true,2,cerr)
      IWaveSpace iw_sp(*pars,ic,true,2,cerr);
      shared_ptr< Space<ireal> > tmp = Space<ireal>::clonePtr(iw_sp);
      rng_ptr = dynamic_pointer_cast< IWaveSpace >(tmp);
	  
      if( !(rng_ptr) ){
	RVLException e;
	e << "failed to make rng_ptr\n";
	throw e;
      }
    }
    catch(RVLException &e){
      e << "ERROR from MPS_to_RHS::set_rng_ptr\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void MPS_to_RHS::sanity_checks() const{
  //----------------------------------------------------------------------//
    try{

      //cheking that s_pos is inside spatial grid
      vector<Rtuple> s_pos = dom.get_s_pos();

      for( int i=0; i<s_pos.size(); i++){
	for( int j=0; j<3;j++){//g.dim; j++){
	  float g_end = g.axes[j].n*g.axes[j].d + g.axes[j].o;
	  if( s_pos[i].coor[j]<g.axes[j].o ||
	      s_pos[i].coor[j]>g_end ){
	    RVLException e;
	    e << "source position is out of bounds!\n"
	      << "s_pos["<<i<<"].coor["<<j<<"]="<<s_pos[i].coor[j]<<"\n";
	    e << "Grid interval is ["<<g.axes[j].o<<","<<g_end<<"]\n";
	    e << "g.axes["<<j<<"].o="<<g.axes[j].o<<"\n"
	      << "g.axes["<<j<<"].d="<<g.axes[j].d<<"\n"
	      << "g.axes["<<j<<"].n="<<g.axes[j].n<<"\n";
	    throw e;
	  }
	}
      }
      
      //checking spatial dimension of grid vs MPS_Space
      if( g.dim!=dom.get_s_dim() ){
	RVLException e;
	e << "spatial dimension of MPS_Space and rsf grid does not match!\n";
	throw e;
      }     
      
      //checking that ntr in RHS is compatible with MPS.
      //loop over RHS sources
      for( int i_rhs=0; 
	   i_rhs < rng_ptr->getKeys().size(); 
	   i_rhs++ ){
	
	string rhs_file = valparse<string>(*pars,rng_ptr->getKeys()[i_rhs]);
	int ntr = get_ntr(rhs_file);

	//expected number of traces
	int NTR = dom.get_s_pos().size(); //number of sources
	for(int i=0; i<g.dim; i++) 
	  //MOD
	  NTR *= a_ord + dom.get_d_ord().coor[i];
	  //NTR *= a_ord + dom.get_MPS_ord();

	if(ntr!=NTR){
	  RVLException e;
	  e << "number of tracs in rhs_file="<<rhs_file<<" is not compatible!\n"
	    << "counted ntr  = "<< ntr <<"\n"
	    << "expected ntr = "<< NTR <<"\n";
	  throw e;
	}
      }
    }
    catch(RVLException &e){
      e << "ERROR from MPS_to_RHS::sanity_checks()\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void MPS_to_RHS::param_set( PARARRAY &locpars, 
			      Vector<ireal> const &x, 
			      Vector<ireal> const &y ) const{
  //----------------------------------------------------------------------//
    try{
      //setting x
      AssignParams apx(locpars,
		       get_MPS_keys().MPS_file);
      x.eval(apx);
      
      //setting y
      Components<ireal> cy(y);
      vector<string> RHS_keys = rng_ptr->getKeys();

      for( int i=0; i<RHS_keys.size(); i++) {
        AssignParams apy(locpars,RHS_keys[i]);
        cy[i].eval(apy);
      }
    }
    catch(RVLException &e){
      e << "ERROR from MPS_to_RHS::param_set!\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  string MPS_to_RHS::get_dom_filename() const{
  //----------------------------------------------------------------------//
    try{
      string fn = dom.get_coreFilename();
      return fn;
    }
    catch(RVLException &e){
      e << "ERROR from MPS_to_RHS::get_dom_filename\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  vector<string> MPS_to_RHS::get_rng_filenames() const{
  //----------------------------------------------------------------------//
    try{
      vector<string> keys = rng_ptr->getKeys();
      vector<string> fns(keys.size());
      for( int i=0; i<keys.size(); i++)
	fns[i] = valparse<string>(*pars,keys[i]);
      return fns;
    }
    catch(RVLException &e){
      e << "ERROR from MPS_to_RHS::get_dom_filename\n";
      throw e;
    }
  }


  //----------------------------------------------------------------------//
  void MPS_to_RHS::print_info(ostream &str) const{
  //----------------------------------------------------------------------//
    try{
      str <<"======================================================\n"
	  <<"Printing MPS_to_RHS info.\n"
	  <<"printing domain: ";
      dom.print_info(str);
      
      str <<"\n"
	  <<"make_RHS = "<< make_RHS <<"\n"
	  <<"printing range: ";
      rng_ptr->write(str);

      str <<"\n"
	  <<"printing grid:\n"
	  <<"    g.dim="<<g.dim<<"\n";
      for(int i=0; i<3;i++){//g.dim; i++){
	str <<"    g.axes["<<i<<"].o="<<g.axes[i].o<<"\n"
	    <<"    g.axes["<<i<<"].d="<<g.axes[i].d<<"\n"
	    <<"    g.axes["<<i<<"].n="<<g.axes[i].n<<"\n";
      }
      str <<"======================================================\n";
    }
    catch(RVLException &e){
      e <<"ERROR from MPS_to_RHS::print_info\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void MPS_to_RHS::apply( const Vector<ireal> &x,
			  Vector<ireal> &y ) const {
  //----------------------------------------------------------------------//
    try{
      SpaceTest(get_MPS_Domain(),x,"MPS_to_RHS::apply (dom)");
      SpaceTest(get_IWave_Range(),y,"MPS_to_RHS::apply (rng)");

      PARARRAY * locpars = ps_new();
      ps_copy(&locpars,*pars);
      param_set(*locpars,x,y);

      // ps_printall(*pars,stderr);
      //fflush(stderr);
     
      //zero output
      y.zero();

      //extracting filenames
      string MPS_file = valparse<string>( *locpars, mks.MPS_file );

      vector<string> RHS_files( rng_ptr->getKeys().size() );
      for( int i=0; i<RHS_files.size(); i++){
	RHS_files[i] = valparse<string>( *locpars, (rng_ptr->getKeys())[i] );
      }
      
#ifdef IWAVE_USE_MPI
      if( retrieveGlobalRank()==0 ){
#endif

	//function writes out RHS, given pars and iokeys
	MPS_to_RHS_kern( a_ord,
			 dom.get_d_ord(),
			 g,
			 dom.get_s_pos(),
			 stencils,
			 dom.get_MPS_basis(),
			 MPS_file,
			 RHS_files );

#ifdef IWAVE_USE_MPI
      }
      MPI_Barrier(retrieveGlobalComm());
#endif
    }
    catch( RVLException &e ){
      e << "ERROR from MPS_to_RHS::apply!\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void MPS_to_RHS::applyAdj( const Vector<ireal> &y, 
			     Vector<ireal> &x ) const {
  //----------------------------------------------------------------------//
    try{
      SpaceTest(get_MPS_Domain(),x,"MPS_to_RHS::applyAdj (dom)");
      SpaceTest(get_IWave_Range(),y,"MPS_to_RHS::applyAdj (rng)");
      
      PARARRAY * locpars = ps_new();
      ps_copy(&locpars,*pars);
      param_set(*locpars,x,y);
     
      //zero output
      x.zero();

      //extracting filenames
      string MPS_file = valparse<string>( *locpars, mks.MPS_file );

      vector<string> RHS_files( rng_ptr->getKeys().size() );
      for( int i=0; i<RHS_files.size(); i++){
	RHS_files[i] = valparse<string>( *locpars, (rng_ptr->getKeys())[i] );
      }
      
#ifdef IWAVE_USE_MPI
      if( retrieveGlobalRank()==0 ){
#endif
	//function writes out RHS, given pars and iokeys
	RHS_to_MPS_kern( a_ord,
			 dom.get_d_ord(),
			 g,
			 dom.get_s_pos(),
			 stencils,
			 dom.get_MPS_basis(),
			 MPS_file,
			 RHS_files );
#ifdef IWAVE_USE_MPI
      }
      MPI_Barrier(retrieveGlobalComm());
#endif
    }
    catch( RVLException &e ){
      e << "ERROR from MPS_to_RHS::applyAdj!\n";
      throw e;
    }
  }
  
}//end TSOpt
