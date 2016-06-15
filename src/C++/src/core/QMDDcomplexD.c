/***************************************************************

Complex number defnitions and routines for QMDD using
doubles for the real and imaginary part of a complex number.

January 28, 2008
Michael Miller
University of Victoria
Victoria, BC 
CANADA V8W 3P6
mmiller@cs.uvic.ca

****************************************************************/

/****************************************************************

The basic idea is that the required complex values are 
stored in a lookup table.	

The value 0 is always in slot 0 and the value 1 is always in slot
1 so that for those two values the index corresponds to the value.

Current implementation uses simple linear searching to find a value.

QMDDinit (in QMDDpackage.c which is the initialization routine that
must be called before the other package routines are used) invokes 
QMDDinitCtable()

*****************************************************************/

#define DEFINE_COMPLEX_H_VARIABLES
#include "QMDDcomplex.h"
#include <julia.h>
#include <julia_threads.h>

void QMDDpause(void);

complex gCvalue(int ai)
{
  complex ret;

  jl_function_t *f = jl_get_function(jl_current_module,"gCvalue");
  jl_value_t *a1 = jl_box_int32(ai);
  jl_array_t *v1 = (jl_array_t *) jl_call1(f,a1);
  long double *xval = (long double *) jl_array_data(v1);
  ret.r = xval[0];
  ret.i = xval[1];

  return ret; 
}

int Clookup(complex c)
{
  int ret;
  jl_value_t *ar_t = jl_apply_array_type(jl_float64_type, 1);
  jl_array_t *x = jl_alloc_array_1d(ar_t, 2);
  JL_GC_PUSH1(&x);

  double *xData = (double *) jl_array_data(x);
  
  xData[0] = c.r;
  xData[1] = c.i;

  jl_function_t *f = jl_get_function(jl_current_module,"Clookup1");
  jl_value_t *rr = jl_call1(f, (jl_value_t *)x);
  ret = jl_unbox_int32(rr);

  JL_GC_POP();

  return ret;
} 

/**************************************

    Routines
    
**************************************/

void Cprint(complex c, std::ostream &os)
{
        if(c.r >= 0)
	  os << " ";
	os << c.r;
	if (c.i > 0)
	   os << "+" << c.i << "i";
	if (c.i < 0)
	  os << c.i << "i";
}

void Cprint(complex c)
// print a complex value
{
	std::ostringstream oss;
	Cprint(c, oss);
	std::cout << oss.str();  
  /*if(c.r==1.0) printf("1"); else if(c.r==0.0) printf("0"); else printf("%f",c.r);
  if(c.i>0.0) printf("+");
  if(c.i!=0.0) printf("%fi ",c.i);*/
}

double QMDDcos(double a, double b)
{
  jl_function_t *f = jl_get_function(jl_current_module,"QMDDcos");
  jl_value_t *v1 = jl_box_float64(a);
  jl_value_t *v2 = jl_box_float64(b);
  jl_value_t *ret = jl_call2(f,v1,v2);
  return jl_unbox_float64(ret);
}

double QMDDsin(double a, double b)
{
  jl_function_t *f = jl_get_function(jl_current_module,"QMDDsin");
  jl_value_t *v1 = jl_box_float64(a);
  jl_value_t *v2 = jl_box_float64(b);
  jl_value_t *ret = jl_call2(f,v1,v2);
  return jl_unbox_float64(ret);
}

int Ceq(int a, int b)
{  
  jl_function_t *f = jl_get_function(jl_current_module,"Ceq");
  jl_value_t *v1 = jl_box_int32(a);
  jl_value_t *v2 = jl_box_int32(b);
  jl_value_t *ret = jl_call2(f,v1,v2);
  return jl_unbox_float64(ret);
}

int Cgt(int a, int b)
// returns 1 if |a|>|b|
// returns 0 if |b|>|a|
// returns angle(a)<angle(b)
// where angle is the angle in polar coordinate representation
{  
  jl_function_t *f = jl_get_function(jl_current_module,"Cgt");
  jl_value_t *v1 = jl_box_int32(a);
  jl_value_t *v2 = jl_box_int32(b);
  jl_value_t *ret = jl_call2(f,v1,v2);
  return jl_unbox_float64(ret);
}

int Cgt_new(int a, int b)
{  
  jl_function_t *f = jl_get_function(jl_current_module,"Cgt_new");
  jl_value_t *v1 = jl_box_int32(a);
  jl_value_t *v2 = jl_box_int32(b);
  jl_value_t *ret = jl_call2(f,v1,v2);
  return jl_unbox_float64(ret);
}

int Clt(int a, int b)
// analogous to Cgt
{
  jl_function_t *f = jl_get_function(jl_current_module,"Clt");
  jl_value_t *v1 = jl_box_int32(a);
  jl_value_t *v2 = jl_box_int32(b);
  jl_value_t *ret = jl_call2(f,v1,v2);
  return jl_unbox_float64(ret);
}

complex Cmake(long double r,long double i)
// make a complex value
{
  complex c;
  
  c.r=r;
  c.i=i;
  return(c);
}

complex CmakeOne(void)
{
  return(Cmake(1.0,0.0));
}

complex CmakeZero(void)
{
  return(Cmake(0.0,0.0));
}

complex CmakeMOne(void)
{
  return(Cmake(-1.0,0.0));
}

long double Qmake(int a, int b,int c)
// returns the complex number equal to (a+b*sqrt(2))/c
// required to be compatible with quadratic irrational-based 
// complex number package
{
  jl_function_t *f = jl_get_function(jl_current_module,"Qmake");
  jl_value_t *v1 = jl_box_int32(a);
  jl_value_t *v2 = jl_box_int32(b);
  jl_value_t *ret = jl_call2(f,v1,v2);
  return jl_unbox_float64(ret);
}

void QMDDcomplexInit(void)
// initialization
{
  jl_eval_string("QMDDinitComplex()");
}

void QMDDcvalue_table_list(void)
// print the complex value table entries
{
}

complex Conj(complex c)
// return complex conjugate
{
  c.i=-c.i;
  return(c);
}

// basic operations on complex values
// meanings are self-evident from the names
// NOTE arguments are the indices to the values 
// in the complex value table not the values themselves

int Cnegative(int a)
{
  jl_function_t *f = jl_get_function(jl_current_module,"Cnegative");
  jl_value_t *v1 = jl_box_int32(a);
  jl_value_t *ret = jl_call1(f,v1);
  return jl_unbox_float64(ret);
}

int Cadd(int ai,int bi)
{
  jl_function_t *f = jl_get_function(jl_current_module,"Cadd");
  jl_value_t *v1 = jl_box_int32(ai);
  jl_value_t *v2 = jl_box_int32(bi);
  jl_value_t *ret = jl_call2(f,v1,v2);
  return jl_unbox_float64(ret);
}

int Csub(int ai,int bi)
{
  jl_function_t *f = jl_get_function(jl_current_module,"Csub");
  jl_value_t *v1 = jl_box_int32(ai);
  jl_value_t *v2 = jl_box_int32(bi);
  jl_value_t *ret = jl_call2(f,v1,v2);
  return jl_unbox_float64(ret);
}

int Cmul(int ai,int bi)
{
  jl_function_t *f = jl_get_function(jl_current_module,"Cmul");
  jl_value_t *v1 = jl_box_int32(ai);
  jl_value_t *v2 = jl_box_int32(bi);
  jl_value_t *ret = jl_call2(f,v1,v2);
  return jl_unbox_float64(ret);
}

int CintMul(int a,int bi)
{
  jl_function_t *f = jl_get_function(jl_current_module,"CintMul");
  jl_value_t *v1 = jl_box_int32(a);
  jl_value_t *v2 = jl_box_int32(bi);
  jl_value_t *ret = jl_call2(f,v1,v2);
  return jl_unbox_float64(ret);
}

int Cdiv(int ai,int bi)
{
  jl_function_t *f = jl_get_function(jl_current_module,"Cdiv");
  jl_value_t *v1 = jl_box_int32(ai);
  jl_value_t *v2 = jl_box_int32(bi);
  jl_value_t *ret = jl_call2(f,v1,v2);
  return jl_unbox_float64(ret);
}

void QMDDmakeRootsOfUnity(void)
{
  jl_eval_string("QMDDmakeRootsOfUnity()");
}

/// by PN: returns the absolut value of a complex number
int CAbs(int a)
{
  jl_function_t *f = jl_get_function(jl_current_module,"CAbs");
  jl_value_t *v1 = jl_box_int32(a);
  jl_value_t *ret = jl_call1(f,v1);
  return jl_unbox_float64(ret);
}

///by PN: returns whether a complex number has norm 1
int CUnit(int a)
{
  jl_function_t *f = jl_get_function(jl_current_module,"CUnit");
  jl_value_t *v1 = jl_box_int32(a);
  jl_value_t *ret = jl_call1(f,v1);
  return jl_unbox_float64(ret);
}

