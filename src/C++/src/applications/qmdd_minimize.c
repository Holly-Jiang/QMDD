/** QMDD minimization using sifting. 
 *  Reads .real-file, builds its QMDD representation
 *  and applies sifting until no more improvement is achieved.
 */

  /**************************************************************************/
  /*  output is:  name of benchmark, initial node count,
   *	          CPU time for first run of sifting,
   *              # of changes to vertex weights during sifting,
   *              # of non-trivial vertex weights when sifting is completed,
   *              max # of nodes during sifting, 
   * 		  min # of nodes during sifting, , ,
   *	          # of sifting runs until no more improvement,
   *              min # of nodes at the end, 
   ***************************************************************************/	
   
#define __LINUX__     // must be __LINUX__ or __WINDOWS__ depending on platform

#ifdef __LINUX__

#include <csignal>                                     
#include <stdexcept>  
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdint.h>
//#include <fpu_control.h>

#define VERBOSE_MAIN 0

static inline int memReadStat(int field)
{
   char    name[256];
   pid_t pid = getpid();
   sprintf(name, "/proc/%d/statm", pid);
   FILE*   in = fopen(name, "rb");
   if (in == NULL) return 0;
   int     value;
   for (; field >= 0; field--)
     fscanf(in, "%d", &value);
   fclose(in);
   return value;}

   uint64_t memUsed() { return (uint64_t)memReadStat(0) * (uint64_t)getpagesize(); }

#endif


#ifdef __WINDOWS__

#define gethostname(x,y) strcpy(x,"unknown")
#include <sstream>
#include <fstream>

#define memUsed() 0
//typedef long uint64_t;
#endif

/**************************************************************************/

#include "../core/QMDDcore.h"  // include the QMDDpackage and related files

/**************************************************************************/
/*   Read a Specification file (RevLib real or spec formats)              */
/**************************************************************************/

QMDDrevlibDescription ReadSpecification(char fname[],QMDDrevlibDescription c,int match)
{
  QMDDrevlibDescription spec;
  char ftype[8];
  int i,j;
  
  // Reverse input filename to check the ending: real/spec
  ftype[0]='\0';
  i=strlen(fname)-1;
  j=0;
  while(i>=0&&fname[i]!='.')
    ftype[j++]=fname[i--];
  ftype[j]='\0';
  
  if(strcmp(ftype,"laer")==0)
    spec=QMDDcircuitRevlib(fname,c,match);
  else {
    printf("*** Invalid file type: %s\n",fname);
    exit(8);
  }
  
  return(spec);
}

/**************************************************************************/
/*   Main Procedure                                                       */
/**************************************************************************/   

int main(int argc, char *argv[] ) 
{
  
  QMDDrevlibDescription circ[2];
  
  long otime;
  int i;
  char argvn[2][64];
  
  /**************************************************************************/
  /*   Argument setup                                                       */
  /**************************************************************************/        
  if(argc==1)
  {
    printf("Plase enter name of input circuit file: ");
    scanf("%s",argvn[1]);
  }
  else for(i=0;i<2;i++) strcpy(argvn[i],argv[i]);
  
  /**************************************************************************/
  /*   Initialize QMDD package                                              */
  /**************************************************************************/
  otime=cpuTime();
  
  QMDDinit(VERBOSE_MAIN);   // 0 = no verbose, 1 = verbose 
  
  if (VERBOSE_MAIN) {
    printf("Initialization time (QMDD package): ");
    printCPUtime(cpuTime()-otime);
    printf("\n");
  }
  
  /**************************************************************************/
  /*   Read the input specification                                         */
  /**************************************************************************/
  
  circ[0]=ReadSpecification(argvn[1],circ[0],0);
  if (circ[0].n == 0) {
    printf("Failed to open file.\n");
    return 1; 
  } else {
    if (VERBOSE_MAIN) printf("Reading circuit %s with %d lines was successful.\n",argvn[1], circ[0].n);
  }
  
  QMDDedge basic = circ[0].e;
  
   char outputbuffer[50];
   int largest, lastResult = ActiveNodeCount;
   int numberOfSifts = 0;
   
   std::ostringstream oss, ossSift;
   
   // Remove "Benchmarks/" from filename
   
   if (strncmp("Benchmarks/", argvn[1], 11) == 0)
     sprintf(outputbuffer, "%20s; %5d; ", &argvn[1][11], ActiveNodeCount);
   else
     sprintf(outputbuffer, "%30s; %5d; ", argvn[1], ActiveNodeCount);
   oss << outputbuffer; 
   
   /****** FIRST SIFTING ********/
   ossSift.str("");
   largest = QMDDsift(circ[0].n, &basic, &circ[0], ossSift);
   oss << ossSift.str();
   sprintf(outputbuffer, " %5d; %5d;", largest, ActiveNodeCount); 
   oss << outputbuffer;
   
   while(ActiveNodeCount != lastResult) {
     numberOfSifts++;
     lastResult = ActiveNodeCount; 
     ossSift.str("");
     largest = QMDDsift(circ[0].n, &basic, &circ[0], ossSift);
   }
   sprintf(outputbuffer, ";; %5d; %5d;", numberOfSifts,ActiveNodeCount); 
   oss << outputbuffer;
   
   oss <<  std::endl;
   std::cout << oss.str();
   exit(0);
   
}	
	
