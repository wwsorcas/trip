// my_waveform.hh
// Author: Mario J. Bencomo
// last modified 11/30/17

#include <su.h>


//--------------------------------------------------------------------
float my_cinf(float t, float w){
//--------------------------------------------------------------------
  float w2 = w*0.5;
  float val = 0.0f;
  if( t>-w2 && t<w2 ){
    val  = exp( ( -w2*w2 + (t+w2)*(w2-t) )/( w2*w2*(t+w2)*(w2-t) ) );
  }
  return val;
}

//--------------------------------------------------------------------
float my_gauss(float t, float fpeak){
//--------------------------------------------------------------------
  float s2 = PI*PI*fpeak*fpeak;
  return exp( -s2*t*t );
}

//--------------------------------------------------------------------
float my_dgauss(float t, float fpeak){
//--------------------------------------------------------------------
  float s2 = PI*PI*fpeak*fpeak;
  return -(2*s2*t)*my_gauss(t,fpeak);
}

//--------------------------------------------------------------------
float my_ricker(float t, float fpeak){
//--------------------------------------------------------------------
  float s2 = PI*PI*fpeak*fpeak;
  return ( 1 - 2*s2*t*t ) * my_gauss(t,fpeak);
}

//--------------------------------------------------------------------
float my_ddgauss(float t, float fpeak){
//--------------------------------------------------------------------
  float s2 = PI*PI*fpeak*fpeak;
  return -(2*s2)*my_gauss(t,fpeak) - (2*s2*t)*my_dgauss(t,fpeak);
}

//--------------------------------------------------------------------
float my_dddgauss(float t, float fpeak){
//--------------------------------------------------------------------
  float s2 = PI*PI*fpeak*fpeak;
  return -(4*s2)*my_dgauss(t,fpeak) - (2*s2*t)*my_ddgauss(t,fpeak);
}

//--------------------------------------------------------------------
float my_ramp(float t, float w, float bot, float top){
//--------------------------------------------------------------------
  float w2 = w*0.5;
  float m = (top-bot)/w;
  float b = bot + m*w2;
  float val;

  if( t<=-w2 ){
    val = bot;
  }
  else 
  if( t>=w2 ){
    val = top;
  }
  else{
    val = m*t + b;
  }
  return val;
}

//--------------------------------------------------------------------
float my_cinf_ramp(float t, float w, float bot, float top){
//--------------------------------------------------------------------
  float w2 = w*0.5;
  float val = 0.0f;
  if( t>-w2 ){
    val = top-bot;    
    if( t<0 )
      val  *= exp( ( -w2*w2 + (t+w2)*(w2-t) )/( w2*w2*(t+w2)*(w2-t) ) );
  }
  val += bot;
  return val;
}

//--------------------------------------------------------------------
float my_wpack(float t, float w, float fpeak){
//--------------------------------------------------------------------
  float val = my_cinf(t,w);
  val *= cos(2*PI*fpeak*t);
  return val;
}


//--------------------------------------------------------------------
float my_waveform(int type, float t, float fpeak){
//--------------------------------------------------------------------
  float output;
  if( type==0 ){
    output = my_gauss(t,fpeak);
  }
  else 
  if( type==1 ){
    output = my_dgauss(t,fpeak);
  }
  else 
  if( type==2 ){
    output = my_ddgauss(t,fpeak);
  }
  else 
  if( type==3 ){
    output = my_dddgauss(t,fpeak);
  }
  else{
    output = -666;
  }
  return output;
}


/*
//--------------------------------------------------------------------
float my_cinf(float dt, int it, int it0, int itF){
//--------------------------------------------------------------------
  if( it<=it0 || it>=itF )
    return 0.0;
  else{
    float h = (itF+it0)*0.5;
    float tmp = exp( -1/( dt*dt*(it-it0)*(itF-it) ) );
    tmp *= exp( 1/( dt*dt*(h-it0)*(itF-h) ) );
    return tmp;
  }
}

//--------------------------------------------------------------------
float my_gauss(float dt, int it, int cit, float fpeak){
//--------------------------------------------------------------------
  float t = dt*(it-cit);
  return my_gauss(t,fpeak);
}

//--------------------------------------------------------------------
float my_dgauss(float dt, int it, int cit, float fpeak){
//--------------------------------------------------------------------
  float t = dt*(it-cit);
  return my_dgauss(t,fpeak);
}

//--------------------------------------------------------------------
float my_ricker(float dt, int it, int cit, float fpeak){
//--------------------------------------------------------------------
  float t = dt*(it-cit);
  float s2 = PI*PI*fpeak*fpeak;
  return ( 1 - 2*s2*t*t ) * my_gauss(t,fpeak);
}

//--------------------------------------------------------------------
float my_ddgauss(float dt, int it, int cit, float fpeak){
//--------------------------------------------------------------------
  float t = dt*(it-cit);
  return my_ddgauss(t,fpeak);
}

//--------------------------------------------------------------------
float my_dddgauss(float dt, int it, int cit, float fpeak){
//--------------------------------------------------------------------
  float t  = dt*(it-cit);
  return my_dddgauss(t,fpeak);
}
*/


