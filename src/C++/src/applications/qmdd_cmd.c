/** Command-Line interface for the QMDD-package. 
 *  Reads .real-files and builds its QMDD-representation.
 *  Gates can be added on the fly, output to DOT-file or matrix / nodelist on stdout.
 *  Typical commands like transpose / norm...
 *  See README.
 * 
 *  extended by Philipp Niemann, August-November 2012
 *  <pniemann@informatik.uni-bremen.de>
 */



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
    printf("*** ERROR: invalid file type: %s\n",fname);
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
  int autoshow = 0;				//  flag that indicates whether a dot-output is created after each applied gate 
  QMDDedge temp_dd, oldbasic; 
  char buffer[50];
  char bak_buffer[50];
  
  
  printf("*** QMDD command line interface. ***\n");   
  
  /**************************************************************************/
  /*   Argument setup                                                       */
  /**************************************************************************/        
  if(argc==1)
  {
    printf("Please enter name of input circuit file: ");
    scanf("%s",argvn[1]);
    fgets(buffer, 50, stdin);
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
    printf("*** Failed to open file.\n");
    return 1; 
  } else {
    printf("*** Reading circuit '%s' with %d lines was successful.\n",argvn[1], circ[0].n);
  }
  
  QMDDedge basic = circ[0].e;
  
  
  /**************************************************************************/
  /*   Preparing EXPORT-2-DOT                                               */
  /**************************************************************************/
  
  // DOT-Output Filename konstruieren,
  
  
  char outputFilename[64], siftFilename[64+5];
  char siftCall[70];
  char sysCall[70];
  strcpy (outputFilename, "dots/");
  char *pch;
  pch=strrchr(argvn[1],'.'); 		//Input-Filename bis zum letzten '.' (ohne Dateiendung)
  strncat (outputFilename,argvn[1], pch-argvn[1]);     
  strncat (outputFilename, ".dot", 4);
  
  printf("*** DOT outputs will be written to file: '%s'\n", outputFilename);
  
  strcpy(siftFilename, outputFilename);
  strncat(siftFilename, ".sift", 5);
  
  strcpy(siftCall, "diff -s ");  // -s == --report-identical-files
  strcat(siftCall, outputFilename);
  strcat(siftCall, " ");	
  strcat(siftCall, siftFilename);
  
  // sysCall if XDOT is installed
  /* strcpy(sysCall, "xdot ");
   *	strcat(sysCall, outputFilename);
   *	strcat(sysCall, " &");
   */
  

  //if(VERBOSE_MAIN) printf("Calling QMDDdotExport...\n");
  
  //QMDDdotExport(basic, 10000, siftFilename, circ[0],0);
  
  //if(VERBOSE_MAIN) printf("Done. Printing QMDD...\n");
  
  /**************************************************************************/
  /*   Print QMDD nodelist and matrix representation	                  */
  /**************************************************************************/
  
  //	QMDDprint(basic, 1000);
  
  //	QMDDmatrixPrint2(basic); 
  
  //	QMDDdotExport(basic, 10000, outputFilename, circ[0],autoshow);
  
  /**************************************************************************/
  /*   Successively read and add gates                                      */
  /**************************************************************************/
  
  
  bool LeftRightFlag = false;			// determines whether new gates are applied at the beginning/end of the circuit
  char siftFlag = 0, redoFlag;
  int gateCount = 0;				// counts the number of gates applied to the original circuit
  
  oldbasic = basic;				// oldbasic shall always point to a backup of the original QMDD
  
  printf("\n*** Available circuit lines: ");
    
    for(i=circ[0].n-1;i>=0;i--) 			// lookup up label in table
      printf("%s ",circ[0].line[QMDDorder[i]].variable);
  
  printf("\n*** Type 'help' for a complete list of commands or start adding gates to the circuit.\n");
  
  while(true) {
    /// print a list of line descriptors and read command
    if(siftFlag){
      printf("\n*** New variable order is: ");
      siftFlag = 0;
       for(i=circ[0].n-1;i>=0;i--) 			// lookup up label in table
      printf("%s ",circ[0].line[QMDDorder[i]].variable);
    }
    
      
    printf("\n>");
    
    fgets(buffer, 50, stdin);
    
    printf("\n");
    
    
    if (strncmp("help", buffer,4) == 0){
    printf("*** =============== List of Commands =============\n\n");
 printf("* General Commands:\n");
     printf("* - help 		prints this list of available commands\n");
     printf("* - exit / quit		quits the program\n");
     printf("\n"); 
     printf("* Modify Circuit (add gates 'on the fly'):\n");
     printf("* - lines		prints a list of available circuit lines\n");
     printf("* - undo		repeats the latest gate\n");
     printf("* - chdir		change whether gates are applied at the beginning/end of the circuit\n");
     printf("* - To add a specific gate, type in\n\n");
     printf("		Ln c1 c2 ... cn t\n\n");
printf("		where\n");
printf("			L	Type of gate, e.g. 'T' for Toffoli, 'H' for Hadamard, 'V' for controlled V\n");
printf("			n	Number of connected lines (controls + target)\n");
printf("			c1,...	Labels of control lines (add '-' for a negative control)\n");
printf("			t	Label of target line\n");
printf("		Note: n is omitted for V/V+ gates and implicitly taken as '2', e.g. V x0 x1, V+ x1 x0\n");
     printf("\n");
     printf("* Generate Output:\n");
     printf("* - show		produces DOT-file and makes system-call to produce SVG-file\n");
     printf("* - stat		print statistics of the QMDD package\n");
     printf("* - print		print QMDD in matrix format\n");
     printf("* - printnodes		print QMDD in nodelist format\n");
     printf("* - ashow		en/disable automatic SVG-expor after adding a new gate\n");
     	 printf("* - active		print number of currently active QMDD nodes\n");
	 printf("* - weight		print weight of QMDD root edge\n");
     printf("\n"); 
     //printf("* QMDD operations:\n");
     //printf("* - transpose\n");
     //printf("* - trace		call 'traceX', traces out level X resp. computes general trace if X>5\n");
     //printf("* - norm		take element-wise norm of QMDD\n");
     //printf("\n"); 
     printf("* Checking for Reversibility:\n");
     printf("* - Ucheck		check whether QMDD represents a unitary matrix\n");
     printf("* - Rcheck		check whether QMDD represents a 0/1-matrix (reversible)\n");
     printf("\n");
     printf("* Variable Reordering:\n");
     printf("* - sift		apply sifting to minimize QMDD size [ heuristic reordering according to cost function ]\n");
     printf("* - presift		(before sifting:) store backup QMDD \n");
     printf("* - compsift		(after sifting:) compare resulting QMDD with backup\n");
     printf("* - exactmin		fexhaustive search for best variable ordering (very slow!)\n");
     printf("* - reorder		perform reordering to initial variable order\n");
     printf("* - moveDIR VARIABLE	move VARIABLE one position up/down (DIR=up,down) or to the top/bottom (DIR=top,bottom)\n");
     
      continue;
    }
       
    if (strncmp("exit", buffer,4) == 0 || strncmp("quit", buffer,4) == 0){
      //printf("let's exit, %d", strncmp(exit, buffer, 4));
      break;
    }
    if (strncmp("make", buffer,4) == 0){
      //printf("let's exit, %d", strncmp(exit, buffer, 4));
      system("make&");
      break;
    }
    
    /// some display/output commands
    else if(strncmp("lines", buffer,5) == 0) {
    printf("*** Available circuit lines: ");
    
    for(i=circ[0].n-1;i>=0;i--) 			// lookup up label in table
      printf("%s ",circ[0].line[QMDDorder[i]].variable);
    
      continue;
    }
    
    else if(strncmp("stat", buffer,4) == 0) {
      QMDDstatistics();
      continue;
    }		
    else if(strncmp("show", buffer,4) == 0) {		// dot/svg export
      printf("*** Exporting QMDD to DOT-file '%s' and convert to SVG...",outputFilename);
      QMDDdotExport(basic, 10000, outputFilename, circ[0],1);
      printf("done.\n");
      continue;
    }
    else if(strncmp("printnodes", buffer,9) == 0) {		// print nodelist
      QMDDprint(basic, 500);
      continue;
    }
    else if(strncmp("print", buffer,5) == 0) {		// print matrix
      QMDDmatrixPrint2(basic); 
      continue;
    }		
    else if(strncmp("ashow", buffer,5) == 0) {		// enable/disable auto graph display
      autoshow= 1-autoshow;
      printf("*** Automatic DOT output %s.\n", (autoshow)?"enabled":"disabled");
      continue;
    }
    else if(strncmp("active", buffer,6) == 0) {		
		QMDDprintActive(circ[0].n);
  
		  continue;
		}
    else if(strncmp("weight", buffer,6) == 0) {	// print current top edge weight
		  printf("Top edge weight is ");
		  Cprint(gCvalue(basic.w));
		  printf("\n");
		  continue;  
		}		
    
    /// some CHECKS
    else if(strncmp("Rcheck", buffer,6) == 0) {		// reversible matrix?
      if(basic.w==1&&basic.p->c01) printf("*** QMDD represents a 0-1 matrix\n");
      else printf("** QMDD does NOT represent a 0-1 matrix\n");
      
      continue;
    }
    else if(strncmp("Ucheck", buffer,6) == 0) {		// unitary matrix?
      QMDDedge ut=QMDDmultiply(basic,QMDDconjugateTranspose(basic));
      if(ut.p->ident) printf("*** QMDD represents unitary matrix\n");
      else { 
	printf("*** QMDD does NOT represent a unitary matrix\n");
	QMDDprint(ut, 80);
	oldbasic = basic;
	basic = ut;
      }
      continue;
    }
    else if(strncmp("Bcheck", buffer,6) == 0) {		//check and count block matrices
      printf("Checking for block matrices...\n");
      printf("There are %d block matrices in the current QMDD (global counter is %d).\n", checkBlockMatrices(basic,0), blockMatrixCounter);
      printf("Now there are %d block matrices in the current QMDD.\n", checkBlockMatrices(basic,1));
      continue;
    }
    /// general QMDD operations		
    else if(strncmp("transpose", buffer,5) == 0) {		// transpose QMDD
      printf("Transposing QMDD (basic only, AND oldbasic!)...\n");
      temp_dd = QMDDtranspose(basic);
      QMDDincref(temp_dd);
      QMDDdecref(basic);
      basic = temp_dd;
      oldbasic = QMDDtranspose(oldbasic);
      QMDDincref(temp_dd);
      QMDDdecref(oldbasic);
      oldbasic = temp_dd;
      QMDDprint(basic, 1000);
      continue;
    }
    else if(strncmp("trace", buffer,5) == 0) {		// call 'traceX', traces out level X resp. computes general trace if X>5
      printf("Tracing out QMDD...\n");
      char rem[] = {0,0,0,0,0,0,0,0,0,0};
      rem[buffer[5]-'0'] = 1;
      temp_dd = QMDDtrace(basic, basic.p->v, rem, (buffer[5]>'5'));
      //oldbasic = QMDDtranspose(oldbasic);
      printf("Root weight is: "); Cprint(gCvalue(temp_dd.w)); printf("\n");
      QMDDprint(temp_dd, 1000);
      QMDDmatrixPrint2(temp_dd);
      oldbasic=basic;
      basic=temp_dd;
      continue;
    }
    /// SIFTING commands			      
    else if(strncmp("exactmin", buffer,8) == 0) {  // perform Steinhaus-Johnson-Trotter algorithm for complete survey of all possible variable orders
      SJTalgorithm(basic, circ[0].n);
      continue;
    }
    else if(strncmp("reorder", buffer,7) == 0) { // restores initial variable order
      printf("*** Reordering QMDD back to normal... start with %d nodes...\n", ActiveNodeCount);
      int order[MAXN];
      for (int i=0; i<circ[0].n; i++)
	order[i]=i;
      
      QMDDreorder(order, circ[0].n, &basic);
      printf("done...\n");
      //QMDDprint(basic, 1000);
      
      oldbasic=basic;
      siftFlag = 1;
      continue;
    }
    else if(strncmp("move", buffer,4) == 0) {		// move{up/down/bottom/top} VAR    moves variable VAR to the respective position
      siftFlag = QMDDmoveVariable(&basic, &buffer[4], &circ[0]);
      oldbasic = basic;
      continue;
    }
    else if(strncmp("presift", buffer,7) == 0) {		// store current QMDD in backup file .sift
      printf("Prepare sifting. Store current QMDD in siftfile.\n");
      QMDDdotExport(basic, 10000, siftFilename, circ[0],0);
      continue;
    }
    else if(strncmp("compsift", buffer,8) == 0) {
      printf("Compare sifting. Store current QMDD in outputfile and run diff.\n");
      QMDDdotExport(basic, 10000, outputFilename, circ[0],0);
      if(system(siftCall))
	printf("\nend of diff.\n");
      continue;
    }
    else if(strncmp("sift", buffer,4) == 0) {
      printf("Start sifting... (currently %d nodes)...\n", ActiveNodeCount);
      int oldmin = siftingCostFunction(basic);
      int largest = QMDDsift(circ[0].n, &basic, &circ[0]);
      printf("Sifting completed. CostFunction: Actual %d, Start %d, Largest %d.\n", siftingCostFunction(basic), oldmin, largest);
      printf("done...\n");
      siftFlag = 1;
      oldbasic=basic;
      continue;
    }
    
    /// a few command-line interface specific calls
    else if(strncmp("chdir", buffer,5) == 0) {		// where shall new gates be added to the circuit?
      LeftRightFlag = !LeftRightFlag;
      if (LeftRightFlag)
	printf("*** Add new gates at the END of the circuit (default).\n");
      else
	printf("*** Add new gates at the BEGINNING of the circuit.\n");
      continue;
    }
    else if (strncmp("undo", buffer,4) == 0){			// undo/repeat last gate
      if (gateCount == 0) {
	printf("*** ERROR: undo is impossible at this stage (no gate applied yet).\n");
	continue;
      }
      gateCount -= 2;
      //
      redoFlag = 1;
      //printf("XredoFlag has value %d.\n", redoFlag);
    }		
    else {
      printf("*** in: %s", buffer);
      if (strlen(buffer)<5) {
	printf("input too short (at least 4 characters needed)!\n");
	continue;
      }
    }
    /// ADD new gates
    //printf("redoFlag has value %d.\n", redoFlag);
    // write gate to tempFile
    if (redoFlag == 1) {
      printf("Undo/repeat last gate...\n");
      strcpy(buffer, bak_buffer);
      redoFlag = 0;
      
    }
    
    // read and build gate from tempFile
    //tempFile=openTextFile(tempfileName,'r');
    temp_dd = QMDDreadGateFromString(buffer, &circ[0]);
    //fclose(tempFile);
    if (temp_dd.p == NULL) {
      printf("*** Syntax error.\n");
      continue;
    }
    
    // multiply gate from the right(front) or left (end)
    if (LeftRightFlag)
      basic = QMDDmultiply(temp_dd, oldbasic);
    else
      basic = QMDDmultiply(oldbasic,temp_dd);
    QMDDincref(basic);
    QMDDdecref(oldbasic);
    //}
    
    //QMDDmatrixPrint2(basic); 
    
    
    
    gateCount++;
    
    printf("*** Added gate (#%d): %s", gateCount, buffer);
    printf("*** New QMDD size: %d nodes. \n", ActiveNodeCount);
    
    oldbasic=basic;
    strcpy(bak_buffer, buffer);
    
    //QMDDdotExport(basic, 10000, outputFilename, circ[0],autoshow);
    
  }
  
  printf("\n");
  //QMDDmatrixPrint2(circ[0].e);     
}
