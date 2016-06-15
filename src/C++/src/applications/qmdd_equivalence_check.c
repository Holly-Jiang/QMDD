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
  int i,j;
  char argvn[3][64],hostName[32],equivalent,date[12],curr_time[12];
  
  uint64_t mem_used;
  
  /**************************************************************************/
  /*   Argument setup                                                       */
  /**************************************************************************/        
  if(argc<3)
  {
    if(argc>1)
      printf("*** Not enough arguments.\n*** You may provide circuit files as parameters as follows: equivalence_check FIRST_CIRCUIT SECOND_CIRCUIT\n");
    printf("Please enter name of first circuit file: ");
    scanf("%s",argvn[1]);
    printf("Please enter name of second circuit file: ");
    scanf("%s",argvn[2]);
  }
  else for(i=0;i<3;i++) strcpy(argvn[i],argv[i]);
  
  
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
  
  otime=cpuTime();
  
  /**************************************************************************/
  /*   Input the two specifications                                         */
  /**************************************************************************/
  
  
  circ[0]=ReadSpecification(argvn[1],circ[0],0);
  if (circ[0].n == 0) {
    printf("*** Failed to open file '%s'.\n", argvn[1]);
    return 1; 
  } else {
    if (VERBOSE_MAIN) printf("Reading circuit %s with %d lines was successful.\n",argvn[1], circ[0].n);
  }
  
  /**************************************************************************/
  /*   Same procedure for the 2nd circuit                 */
  /**************************************************************************/
  
  circ[1]=ReadSpecification(argvn[2],circ[0],0);
  if (circ[1].n == 0) {
    printf("*** Failed to open file '%s'.\n", argvn[2]);
    return 1; 
  } else {
    if (VERBOSE_MAIN) printf("Reading circuit %s with %d lines was successful.\n",argvn[2], circ[1].n);
  }
  
  
  
  /**************************************************************************/
  /*   Do equivalence check                                                 */
  /**************************************************************************/
  
  equivalent=QMDDedgeEqual(circ[0].e,circ[1].e);
  
  printf("\n\n************************************************\n\n");
  if(equivalent) printf("The two specifications match.\n");
  else printf("The two specifications do not match.\n");
  printf("\n\n************************************************\n\n");
  
  otime=cpuTime()-otime;
  
  gethostname(hostName,32);
  
  dateToday(date);
  wallTime(curr_time);
  
  
  printf("OUTPUT1: -DATE %s -TIME %s -SOURCE %s -HOST %s\n",date,curr_time,__FILE__,hostName);
  printf("OUTPUT2: -MAXN %d -UTBCKTS %d -CTSIZE %d -TTSIZE %d -GCLIMIT %d -CNTSIZE %d\n",MAXN,NBUCKET,CTSLOTS,TTSLOTS,GCLIMIT1,COMPLEXTSIZE);
  printf("OUTPUT3: ");
  printf("-PEAK %ld -CPU ",QMDDpeaknodecount);
  printCPUtime(otime);
  
  /**********************************************************/
  for(i=strlen(argvn[1])-1;i>=0&&argvn[1][i]!='/';i--);
  if(argvn[1][i]=='/') i++;
  
  for(j=0;argvn[1][i]!='.'&&argvn[1][i]==argvn[2][i];i++,j++)
    circ[0].name[j]=argvn[1][i];
  
  while(circ[0].name[j-1]>='0'&&circ[0].name[j-1]<='9') {j--;i--;}
  
  if(circ[0].name[j-1]=='-'||circ[0].name[j-1]=='_') j--;
  circ[0].name[j]=0;
  
  for(j=0;circ[0].name[j]!=0&&circ[0].name[j]!='_';j++);
  circ[0].name[j]=0;
  
  strcpy(circ[1].name,circ[0].name);
  
  for(j=i;argvn[1][j]!='-'&&argvn[1][j]!='_'&&argvn[1][j]!='.';j++)
    circ[0].no[j-i]=argvn[1][j];
  circ[0].no[j-i]=0;
  circ[0].modified=argvn[1][j]=='_';
  
  for(j=i;argvn[2][j]!='-'&&argvn[2][j]!='_'&&argvn[2][j]!='.';j++)
    circ[1].no[j-i]=argvn[2][j];
  circ[1].no[j-i]=0;
  circ[1].modified=argvn[2][j]=='_';
  /**********************************************************/
  
  printf(" -EXAMPLE ");
  if(circ[1].modified) printf("%s(%s)\n",circ[0].name,circ[0].no);
  else printf("%s(%s,%s)\n",circ[0].name,circ[0].no,circ[1].no);
  if(strcmp(circ[0].name,circ[1].name)) printf("warning circuit name mismatch %s %s\n",circ[0].name,circ[1].name);
  
  #ifdef __LINUX__
  mem_used = memUsed();
  printf("OUTPUT4: -STAT-Memory %.2f MB\n", mem_used / 1048576.0);
  #endif
}
