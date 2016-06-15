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

void QMDDpause(void);

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

#define Ceq(x,y) ((fabs((x.r)-(y.r))<Ctol)&&(fabs((x.i)-(y.i))<Ctol))

long double angle(int a)
// computes angle for polar coordinate representation of Cvalue(a)
{
  complex ca;
  ca=Cvalue(a);
  if(ca.i>=0-Ctol) return(acos(ca.r/Cmag[a]));
  else return(2*Pi-acos(ca.r/Cmag[a]));
}

int Cdiv(int ai,int bi); /* prototype */

int Cgt(int a, int b)
// returns 1 if |a|>|b|
// returns 0 if |b|>|a|
// returns angle(a)<angle(b)
// where angle is the angle in polar coordinate representation
{  
  complex ca,cb;
  if(a==b) return(0);
  ca=Cvalue(a);
  cb=Cvalue(b);
  
  /// BETA: 121017
  // returns the smaller nonzero value
  if (a == 0)
    return(1);
  if (b == 0)
    return(0);
  
  /*int c=a;
  a=b;
  b=c;*/
  
  /// BETA END */
  
  if(Cmag[a]>(Cmag[b]+Ctol)) return(1);
  if(Cmag[b]>(Cmag[a]+Ctol)) return(0);
  //CHANGED by pN 120831
  return((Cangle[a]+Ctol)<Cangle[b]);
  return (0);
}

int Cgt_new(int a, int b)
{  
  complex ca,cb;
  if(a==b) return(0);
  ca=Cvalue(a);
  cb=Cvalue(b);
  if((Cangle[a]+Ctol)<Cangle[b]) return(1);
  return(Cmag[a]>(Cmag[b]+Ctol));
}

int Clt(int a, int b)
// analogous to Cgt
{
  complex ca,cb;
  if(a==b) return(0);
  ca=Cvalue(a);
  cb=Cvalue(b);
  if(Cmag[a]<(Cmag[b]+Ctol)) return(1);
  if(Cmag[b]<(Cmag[a]+Ctol)) return(0);
  return((angle(a)+Ctol)>angle(b));
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
  return(((float)a+((float)b)*sqrt(2.0))/(float)(c));
}

void QMDDinitCtable(void)
// initialize the complex value table and complex operation tables to empty
{
  int i,j;
  
  Ctentries=0;
    
  for(i=0;i<COMPLEXTSIZE;i++)
  {
    for(j=0;j<COMPLEXTSIZE;j++)
      cta[i][j]=cts[i][j]=ctm[i][j]=ctd[i][j]=-1;
  }
  
  if(VERBOSE) printf("\nDouble complex number package initialized\n\n");
}

void QMDDcomplexInit(void)
// initialization
{
  Ctol=1.0e-10;
  Cmag[0]=0;
  Cmag[1]=1;
  QMDDinitCtable();
}

void QMDDcvalue_table_list(void)
// print the complex value table entries
{
  int i;
  
  printf("\nComplex value table: %d entries\n",Ctentries);
  std::cout << "index value Magnitude Angle 1) radian 2) degree" << std::endl;
  for(i=0;i<Ctentries;i++)
  {
    std::cout << i << " ";
    //printf("%d ",i);
    Cprint(Ctable[i]);
    //printf("  || %f %f\n",Cmag[i],angle(i));
    if(i!=0&&i%100==0) QMDDpause();
    std::cout << " || " << Cmag[i]<< " " << angle(i);
    std::cout << " " <<(angle(i)*180/3.141592654) << std::endl;
 }
}

int Clookup(complex c)
// lookup a complex value in the complex value table
// if not found add it
// this routine uses linear searching
{
  int i;

  for(i=0;i<Ctentries;i++)
    if(Ceq(Ctable[i],c)) 
    return(i);
  Ctentries++;
  if(Ctentries>COMPLEXTSIZE)
  {
    printf("Complex mapping table overflow.");
    Cprint(c);
    printf("\n");
    QMDDpause();
    QMDDcvalue_table_list();
    exit(0);
  }
  i=Ctentries-1;
  Ctable[i]=c;
  Cmag[i]=sqrt(c.r*c.r+c.i*c.i);
  Cangle[i]=angle(i);
  return(i);
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
  complex c;
  c=Cvalue(a);
  c.r=-c.r;
  c.i=-c.i;
  return(Clookup(c));
}

int Cadd(int ai,int bi)
{
  complex a,b,r;
  int t;
  
  if(ai==0) return(bi); // identity cases
  if(bi==0) return(ai);

  if(0<=(t=cta[ai][bi])) return(t); // look in computation table
  
  a=Cvalue(ai); // if new compute result
  b=Cvalue(bi); 
  r.r=a.r+b.r;
  r.i=a.i+b.i;
  
  t=cta[ai][bi]=cta[bi][ai]=Clookup(r); // save result
  return(t);
}

int Csub(int ai,int bi)
{
  complex a,b,r;
  int t;
  
  if(bi==0) return(ai); // identity case
  
  if(0<=(t=cts[ai][bi])) return(t); // look in computation table
  
  a=Cvalue(ai);  // if new compute result
  b=Cvalue(bi);
  r.r=a.r-b.r;
  r.i=a.i-b.i;
  
  t=cts[ai][bi]=Clookup(r); // save result
  return(t);
}

int Cmul(int ai,int bi)
{
  complex a,b,r;
  int t;
  
  if(ai==1) return(bi); // identity cases
  if(bi==1) return(ai);
  if(ai==0||bi==0) return(0);
  
  if(0<=(t=ctm[ai][bi])) return(t); // look in computation table
  
  a=Cvalue(ai); // if new compute result
  b=Cvalue(bi);
  r.r=a.r*b.r-a.i*b.i;
  r.i=a.r*b.i+a.i*b.r;
  
  t=ctm[ai][bi]=ctm[bi][ai]=Clookup(r); // save result
  return(t);
}

int CintMul(int a,int bi)
{
  complex r;
  r=Cvalue(bi);
  r.r*=a;
  r.i*=a;
  return(Clookup(r));
}

int Cdiv(int ai,int bi)
{
  complex a,b,r;
  int t;
  long double d;
  
  if(ai==bi) return(1); // equal case
  if(ai==0) return(0); // identity cases
  if(bi==1) return(ai);
  
  if(0<=(t=ctd[ai][bi])) return(t); // check computation table
  
  a=Cvalue(ai); // if new compute result
  b=Cvalue(bi);
  if(b.i==0.0)
  {
    r.r=a.r/b.r;
    r.i=a.i/b.r;
  } else {
    d=b.r*b.r+b.i*b.i;
    r.r=(a.r*b.r+a.i*b.i)/d;
    r.i=(a.i*b.r-a.r*b.i)/d;
  }
  t=ctd[ai][bi]=Clookup(r); // save result
  return(t);
}

void QMDDmakeRootsOfUnity(void)
{
  int i;
  CTa[0]=1;
  CTa[1]=Clookup(Cmake(cosl(2*Pi/Radix),sinl(2*Pi/Radix)));
  for(i=2;i<Radix;i++)
    CTa[i]=Cmul(CTa[i-1],CTa[1]);
}

/// by PN: returns the absolut value of a complex number
int CAbs(int a)
{
  int b;
  complex r, s;
  
  if (a<2) return a; // trivial cases 0/1
  
    s=Cvalue(a);
  //printf("CAbs: "); Cprint(s); printf(" is ");
  r.r=Cmag[a];  // Cmag table contains magnitude of complex table entries
  r.i=0;
  b = Clookup(r);
  //Cprint(r);   printf("\n");
  return(b);
}

///by PN: returns whether a complex number has norm 1
int CUnit(int a)
{
 /// BETA 121017
 
 if (a<2)
   return a;
 
 if (Cmag[a] < 1-Ctol)
   return 0;
 else
   return 1;
}

  
