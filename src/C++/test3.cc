#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "src/core/QMDDcore.h"
#include "src/core/QMDDpackage.h"
#include "src/core/QMDDcircuit.h"

int main() {
   clock_t tStart = clock();

   QMDDrevlibDescription spec;

   QMDDinit(1);

   int col[] = {0,1,0}; // |0000010>

   QMDDedge id = QMDDident(0,2); // n - 1 where n is qubit length
   QMDDedge k0 = QMDDextractColumn(col,id);

   //QMDDmatrixPrint2(k0);

   spec = QMDDcircuitRevlib("qFou3.real",spec,0);

   if (spec.n == 0) {
     printf("*** Failed to open file.\n");
     return 1; 
   } else {
     printf("*** Reading circuit with %d lines was successful.\n",
       spec.n);
   }

   QMDDedge qft = spec.e;
   //QMDDmatrixPrint2(qft);

   //QMDDedge HAD = QMDDgate(Hm,1,0,0);
   //QMDDmatrixPrint2(HAD);

   QMDDedge hh = QMDDmultiply(qft,k0);
   QMDDcolumnPrint(hh,3);
  
   std::cout << "Time: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << std::endl;

   return 0;
}

