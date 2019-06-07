// MPS_Space.cc
// Author: Mario J. Bencomo
// last modified: 12/07/16

#include "MPS_Space.hh"
//#define VERBOSE_MJB

namespace TSOpt {
  
  //----------------------------------------------------------------------//
  void add_to_pars( PARARRAY const &pars, 
		    vector<string> keys, 
		    vector<string> vals){
  //----------------------------------------------------------------------//
    try{
     
      if( keys.size()!=vals.size() ){
	RVLException e;
	e << "number of keys does not match number of values ...\n";
	throw e;
      }
      
      for( int i=0; i<keys.size(); i++ ){
	if( ps_slcstring(pars, keys[i].c_str(), vals[i].c_str()) ){
	  RVLException e;
	  e << "trouble setting (key,val) pair "
	    << "= ("<<keys[i]<<","<<vals[i]<<")\n";
	}
      }
      
    }
    catch(RVLException &e){
      e << "ERROR from add_to_pars!\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void add_to_pars( PARARRAY const &pars, 
		    string key, 
		    string val ){
  //----------------------------------------------------------------------//
    try{
      
      if( ps_slcstring(pars, key.c_str(), val.c_str()) ){
	RVLException e;
	e << "trouble setting (key,val) pair "
	    << "= ("<<key<<","<<val<<")\n";
      }
      
    }
    catch(RVLException &e){
      e << "ERROR from add_to_pars!\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void add_to_pars( PARARRAY const &pars, 
		    string key, 
		    int val ){
  //----------------------------------------------------------------------//
    try{
      
      if( ps_slint(pars, key.c_str(), val) ){
	RVLException e;
	e << "trouble setting (key,val) pair "
	    << "= ("<<key<<","<<val<<")\n";
      }
      
    }
    catch(RVLException &e){
      e << "ERROR from add_to_pars!\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  vector<Rtuple> ex_pos( string sufile, bool rcv ){
  //----------------------------------------------------------------------//
    try{

      vector<Rtuple> pos;

#ifdef IWAVE_USE_MPI
      if( retrieveGlobalRank()==0 ){
#endif
	//Reading in file
	FILE *fp = NULL;
	if(!(fp=iwave_const_fopen(sufile.c_str(),"r",NULL,stderr))){
	  RVLException e;
	  e << "failed to open core file:"<<sufile<<"\n";
	  throw e;
	}
	segy tr;
	
	//extracting source positions                                                                                                                                                            
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
	
	while(fgettr(fp,&tr)){
	  
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
	iwave_fclose(fp);
	
#ifdef IWAVE_USE_MPI
      }

      int sz = pos.size();
      MPI_Bcast(&sz,1,MPI_INT,0,retrieveGlobalComm()); 
   
      //broadcasting array containing s_pos to other ranks
      double * arr = new double[sz*3];
      if( retrieveGlobalRank()==0 ){
	for( int i=0; i<sz; i++ ){
	  for( int j=0; j<3; j++ ){
	    arr[j+i*3] = pos[i].coor[j];
	  }
	}
      }
      MPI_Bcast(arr,sz*3,MPI_DOUBLE,0,retrieveGlobalComm()); 

      //reading into s_pos
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
      
      return pos;
    }
    catch(RVLException &e){
      e << "ERROR from ex_pos()!\n";
      throw e;
    }
  }


  //----------------------------------------------------------------------//
  vector<Rtuple> ex_pos_serial( string sufile, bool rcv ){
  //----------------------------------------------------------------------//
    try{

      vector<Rtuple> pos;

      //Reading in file
      FILE *fp = NULL;
      if(!(fp=iwave_const_fopen(sufile.c_str(),"r",NULL,stderr))){
	RVLException e;
	e << "failed to open core file:"<<sufile<<"\n";
	throw e;
      }
      segy tr;
      
      //extracting source positions                                                                                                                                                            
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
      
      while(fgettr(fp,&tr)){
	
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
      iwave_fclose(fp);
   
      return pos;
    }
    catch(RVLException &e){
      e << "ERROR from ex_pos()!\n";
      throw e;
    }
  }

  

  //----------------------------------------------------------------------//
  int get_ntr(string sufile){
  //----------------------------------------------------------------------//
    try{
      
      int ntr=0;

#ifdef IWAVE_USE_MPI      
      if( retrieveGlobalRank()==0 ){
#endif	
	FILE *fp = NULL;
	if( !(fp=iwave_const_fopen(sufile.c_str(),"r",NULL,stderr)) ){
	  RVLException e;
	  e << "failed to open su file "<<sufile<<"\n";
	  throw e;
	}
	
	segy tr;
	fvgettr(fp,&tr);
	ntr=tr.ntr;
	
	if(ntr==0){
	  ntr++;
	  while(fvgettr(fp,&tr)) ntr++;
	}
	iwave_fclose(fp);

#ifdef IWAVE_USE_MPI
      }
      MPI_Bcast(&ntr,1,MPI_INT,0,retrieveGlobalComm());
#endif
      return ntr;
    }
    catch(RVLException &e){
      e << "ERROR from ntr_segy\n";
      throw e;
    }
  }


  //----------------------------------------------------------------------//
  int get_ntr_serial(string sufile){
  //----------------------------------------------------------------------//
    try{
      
      int ntr=0;

      FILE *fp = NULL;
      if( !(fp=iwave_const_fopen(sufile.c_str(),"r",NULL,stderr)) ){
	RVLException e;
	e << "failed to open su file "<<sufile<<"\n";
	throw e;
      }
      
      segy tr;
      fvgettr(fp,&tr);
      ntr=tr.ntr;
      
      if(ntr==0){
	ntr++;
	while(fvgettr(fp,&tr)) ntr++;
      }
      iwave_fclose(fp);
      
      return ntr;
    }
    catch(RVLException &e){
      e << "ERROR from ntr_segy\n";
      throw e;
    }
  }


  //----------------------------------------------------------------------//
  void write_MPS_SEGY( vector<Rtuple> s_pos,
		       int N_basis,
		       string ref_file,
		       string MPS_file ){
  //----------------------------------------------------------------------//
    try{

#ifdef IWAVE_USE_MPI
      if( retrieveGlobalRank()!=0 ){
	RVLException e;
	e << "write_MPS_SEGY should only be called by rank 0!\n";
      }
#endif

      //opening MPS file
      FILE *fp_mps;
      fp_mps = iwave_const_fopen(MPS_file.c_str(),"w",NULL,stderr); 
      if(fp_mps==NULL){
	RVLException e;
	e << "failed to open MPS file!\n";
	throw e;
      }

      //opening ref file
      FILE *fp_ref = iwave_const_fopen(ref_file.c_str(),"r",NULL,stderr); 
      if(fp_ref==NULL){
	RVLException e;
	e << "failed to open reference file!\n";
	throw e;
      }

      segy tr;
      fgettra(fp_ref,&tr,0);
      
      //zeroing trace
      for( int it=0; it<tr.ns; it++)
	tr.data[it]=0;

      int i=0;
      //looping over sources
      for( int i_src=0; i_src<s_pos.size(); i_src++){

	//looping over MPS bases
	for( int i_base=0; i_base<N_basis; i_base++){
	  
	  tr.tracr=++i;
	  tr.tracl=i;

	  tr.gelev=0;
	  tr.gx=0;
	  tr.gy=0;

	  tr.selev=-s_pos[i_src].coor[0];
	  tr.sx   =s_pos[i_src].coor[1];
	  tr.sy   =s_pos[i_src].coor[2];
	  
	  fvputtr(fp_mps,&tr);
	}
      }
	
      fflush(fp_mps);

      iwave_fclose(fp_mps);
      iwave_fclose(fp_ref);
    }
    catch(RVLException &e){
      e << "ERROR from write_MPS_SEGY\n";
      throw e;
    }
  }

  //////////////////////////////////////////////////////////////////////////
  // MPS_Space ////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////  

  //----------------------------------------------------------------------//
  void MPS_Space::gen_SEGY_core( PARARRAY const &pars ){
  //----------------------------------------------------------------------//
    try{

      string MPS_file;
      string ref_file;
      string hdr_file;
      
      try{
	  MPS_file = valparse<string>(pars,mks.MPS_file);
	  ref_file = valparse<string>(pars,mks.MPS_ref);
	  hdr_file = valparse<string>(pars,mks.hdr_files[0]);
      }
      catch(RVLException &e){
	cerr << "Missing MPS, reference, or header filename.\n";
	throw e;
      }
	
      //extract source positions from header file
      s_pos = ex_pos(hdr_file,false);

#ifdef IWAVE_USE_MPI
      if( retrieveGlobalRank()==0 ){
#endif	
	//generate MPS SEGY file
	write_MPS_SEGY( s_pos,
			get_MPS_basis().size(),
			ref_file,
			MPS_file );
#ifdef IWAVE_USE_MPI
      }
      MPI_Barrier(retrieveGlobalComm());
#endif

      //allocating core
#ifdef IWAVE_USE_MPI
      MPISEGYSpace segy_sp(MPS_file,"MPS_file",retrieveGlobalComm(),cerr);
      core = RVL::Space<float>::export_clone(segy_sp);
#else
      SEGYSpace segy_sp(MPS_file,"MPS_file");
      core = RVL::Space<float>::export_clone(segy_sp);
#endif      
      
      if(core==NULL){
	RVLException e;
	e << "failed to allocate core space\n";
	throw e;
      }

      //sanity checks
      //sanity_checks();

    }
    catch(RVLException &e){
      e << "ERROR from MPS_Space::gen_SEGY_core()\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void MPS_Space::sanity_checks() const{
  //----------------------------------------------------------------------//
    try{
      string corefn = this->get_coreFilename();
      int ntr = get_ntr(corefn);
      int NTR = this->get_dim()*this->get_s_pos().size(); 

      if( ntr!=NTR ){
	RVLException e;
	e << "Number of traces in core file ("<<corefn<<") "
	  << "does not match what is expected!\n"
	  << "ntr in file = "<< ntr << "\n"
	  << "ntr expected =" << NTR <<"\n";
      }
    }
    catch(RVLException &e){
      e << "ERROR from MPS_Space::sanity_checks!\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  MPS_Space::MPS_Space( MPS_Space const &sp)
    : core(NULL),
      mks(sp.mks),
      s_pos(sp.s_pos),
      s_dim(sp.s_dim) {
  //----------------------------------------------------------------------//
    try{

      //copying of core space
#ifdef IWAVE_USE_MPI
      MPISEGYSpace const* segy_sp = NULL;
      segy_sp = dynamic_cast<MPISEGYSpace const*>(sp.core);
      core = RVL::Space<float>::export_clone(*segy_sp);
#else
      SEGYSpace const* segy_sp = NULL;
      segy_sp = dynamic_cast<SEGYSpace const*>(sp.core);
      core = RVL::Space<float>::export_clone(*segy_sp);
#endif
      
      if( core==NULL ){
	RVLException e;
	e << "failed to copy core\n";
	throw e;
      }
    }
    catch(RVLException &e){
      e << "ERROR from MPS_Space copy constructor\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  MPS_Space::MPS_Space( MPS_KEYS _mks,
			PARARRAY const& pars, 
			bool make )
    : core(NULL), mks(_mks) {
  //----------------------------------------------------------------------//
    try{

      //setting spatial dimension
      string grid_file = valparse<string>(pars,mks.grid_file);
      grid g; init_default_grid(&g);
      read_grid(&g,grid_file.c_str(),stderr);
      s_dim = g.dim;

      //case where MPS file exists
      if(!make){

	//allocating core
	string MPS_file = valparse<string>(pars,mks.MPS_file);
#ifdef IWAVE_USE_MPI
	MPISEGYSpace segy_sp(MPS_file,"MPS_file",retrieveGlobalComm(),cerr);
	core = RVL::Space<float>::export_clone(segy_sp);
#else
	SEGYSpace segy_sp(MPS_file,"MPS_file");
	core = RVL::Space<float>::export_clone(segy_sp);
#endif      

	if(core==NULL){
	  RVLException e;
	  e << "failed to allocate core space\n";
	  throw e;
	}

	//setting source positions
	s_pos = ex_pos(MPS_file,false);
  
	if(s_pos.size()==0){
	  RVLException e;
	  e << "no source positions extracted!\n";
	  throw e;
	}
	
	//sanity checks
	//sanity_checks();
      }
      //case where SEGY core is made from scratch
      else{

#ifdef IWAVE_USE_MPI
	if( retrieveGlobalRank()==0 ){
#endif
	  cerr << "WARNING from MPS_Space constructor, \n"
	       << "MPS file is assumed to not exist, \n"
	       << "will make from scratch in subclass constructor!\n";
#ifdef IWAVE_USE_MPI
	}
#endif
      }
    }
    catch(RVLException &e){
      e << "ERROR from MPS_Space constructor\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  MPS_Space::~MPS_Space(){
  //----------------------------------------------------------------------//
    if(core) delete core;
  }

  /*
  //----------------------------------------------------------------------//
  DataContainer * MPS_Space::buildDataContainer() const {
  //----------------------------------------------------------------------//  
    StdProductDataContainer * d = new StdProductDataContainer();
    for (size_t i=0; i<getSize(); i++) {
      SpaceDCF<ireal> f(*(core));
      d->push(f);
    }
    return d;
  }

  //----------------------------------------------------------------------//
  size_t MPS_Space::getSize() const{
  //----------------------------------------------------------------------//
    size_t s = (core==NULL)? 0:1;
    return s;
  }

  //----------------------------------------------------------------------//
  Space<ireal> const & MPS_Space::operator[](size_t i) const{
  //----------------------------------------------------------------------//
    try{
      if( i!=0 ){
	RVLException e;
	e << "index i="<< i <<" must be equal to zero\n";
	throw e;
      }
      return *(core); 
    }
    catch(RVLException &e){
      e << "ERROR from MPS_Space::operator[]\n";
      throw e;
    }
  }
  */

  //----------------------------------------------------------------------//
  DataContainerFactory const & MPS_Space::getDCF() const{
  //----------------------------------------------------------------------//  
    StdSpace<ireal,ireal> const* stdsp =
      dynamic_cast<StdSpace<ireal,ireal> const*>(core);
    return stdsp->getDCF();
  }

  //----------------------------------------------------------------------//
  LinearAlgebraPackage<ireal> const & MPS_Space::getLAP() const{
  //----------------------------------------------------------------------//  
    StdSpace<ireal,ireal> const* stdsp =
      dynamic_cast<StdSpace<ireal,ireal> const*>(core);
    return stdsp->getLAP();
  }

  //----------------------------------------------------------------------//
  string MPS_Space::get_coreFilename() const{ 
  //----------------------------------------------------------------------//
    try{
#ifdef IWAVE_USE_MPI
      MPISEGYSpace const* segy_sp = NULL;
      segy_sp = dynamic_cast<MPISEGYSpace const*>(core);
#else
      SEGYSpace const* segy_sp = NULL;
      segy_sp = dynamic_cast<SEGYSpace const*>(core);
#endif
    
      return segy_sp->getPrototypeFilename();
    }
    catch(RVLException &e){
      e << "ERROR from MPS_Space::get_coreFilename()\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  Ituple MPS_Space::get_d_ord() const{
  //----------------------------------------------------------------------//
    try{
      
      vector<MPS_base> MPS_b = get_MPS_basis();
      Ituple d_ord;
      d_ord.coor[0] = 0;
      d_ord.coor[1] = 0;
      d_ord.coor[2] = 0;
      
      //looping over bases
      for( int i_b=0; i_b<MPS_b.size(); i_b++ ){
	//looping over terms
	for( int i_t=0; i_t<MPS_b[i_b].NT; i_t++ ){
	  //looping of spatial dim
	  for( int d=0; d<get_s_dim(); d++ ){
	    if( MPS_b[i_b].derv[i_t].coor[d]>d_ord.coor[d] )
	      d_ord.coor[d] = MPS_b[i_b].derv[i_t].coor[d];
	  }
	}
      }
      
      return d_ord;
    }
    catch(RVLException &e){
      e << "ERROR from MPS_Space::get_d_ord()\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  int MPS_Space::get_MPS_ord() const{
  //----------------------------------------------------------------------//
    try{
      
      vector<MPS_base> MPS_b = get_MPS_basis();
      int d_ord=0;
      int tmp_ord;

      //looping over bases
      for( int i_b=0; i_b<MPS_b.size(); i_b++ ){
	//looping over terms
	for( int i_t=0; i_t<MPS_b[i_b].NT; i_t++ ){
	  tmp_ord=0;
	  //looping of spatial dim
	  for( int d=0; d<get_s_dim(); d++ ){
	    tmp_ord+=MPS_b[i_b].derv[i_t].coor[d];
	  }
	  if( tmp_ord>d_ord ) d_ord = tmp_ord;
	}
      }
      
      return d_ord;
    }
    catch(RVLException &e){
      e << "ERROR from MPS_Space::get_d_ord()\n";
      throw e;
    }
  }
  
  //----------------------------------------------------------------------//
  void MPS_Space::print_info(ostream &str) const{
  //----------------------------------------------------------------------//
    try{
      str << "\n"
	  << "-----------------------------------------\n"
	  << "Printing out MPS_Space info:\n"
	  << "-----------------------------------------\n"
	  << "     type = "<< get_type() <<"\n"
	  << "    s_dim = "<< s_dim <<"\n"
	  << "  MPS dim = "<< get_dim() <<"\n"
	  << "core segy = "<< get_coreFilename() <<"\n"
	  << "    d_ord = ";
      for(int i=0; i<s_dim; i++){
	str <<get_d_ord().coor[i]<<",";
      }
      str << "\n"
	  << " src pos = \n";
      for( int j=0; j<s_pos.size(); j++){
	str << "     src["<<j<<"]=(";
	for( int i=0; i<s_dim; i++){
	  str << s_pos[j].coor[i] << ",";
	}
	str <<")\n";
      }
      str << "\n";

      //printing each basis
      for(int i=0; i<get_dim(); i++){
	str << ">>>>>>>>>> basis_"<< i <<", NT="<< get_MPS_basis()[i].NT <<"\n";
	for(int j=0; j<get_MPS_basis()[i].NT; j++){
	  str << "           derv=";
	  for(int k=0; k<s_dim; k++){
	    str << get_MPS_basis()[i].derv[j].coor[k] <<",";
	  }
	  str <<"\n";
	} 
	str <<"RHS_flags=";
	for(int j=0; j<get_MPS_basis()[i].RHS_flags.size(); j++)
	  str << get_MPS_basis()[i].RHS_flags[j] <<",";
	str <<"\n";
      }
      str << "-----------------------------------------\n";
    }
    catch(RVLException &e){
      e << "ERROR from MPS_Space::print_info!\n";
      throw e;
    }
  }

}//end TSOpt
