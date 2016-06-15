/***********************************************************************

This file contains routines for constructing QMDD for gates and
for reading reversible and quantum circuits in extended RevLib real 
format.


    Michael Miller
    University of Victoria
    mmiller@uvic.ca
    
    Date: October 2008

***********************************************************************/

#include "QMDDcircuit.h"
/*******************************************************************
    Routines            
*****************************************************************/

int getlabel(char lab[],QMDDrevlibDescription circ,int *cont)
// get a label and return its index in a global Label table
// cont is returned as 1 (positive control/target), 0 (negative control), or -1 if character after label is '\n'

// this is a really quick and dirty implementation that uses
// linear searching
{
  int i,k;
  char ch = lab[strlen(lab)-1]; // last character of lab
  
  
  *cont = 1;
  if (lab[0] == '-'){   // if negative control
  *cont = 0;
    lab = &lab[1];
  //for (int i=0; i<strlen(lab); i++)
  //    lab[i] = lab[i+1];
  }
  
  if (ch == '\n')
    *cont = -1;				//set continuation control
    
    if(lab[0]==' ') return(-1);			// no label
    
    for(i=0;i<circ.n;i++) 			// lookup up label in table
    {
      if(0==strcmp(lab,circ.line[i].variable)) 
      {
	return(i); // return if found
      }
    }
    printf("label not found: %s\n",lab);
    return -1;
}


QMDDedge QMDDreadGateFromString(char *str, QMDDrevlibDescription *circ)
{
	int cont,i,j,k,m,n,t,line[MAXN];
	QMDDedge f,f2;
	char ch1,ch2,ch3;
	char token[MAXSTRLEN];

	int div,pc[MAXN];
	int sign;

	cont = 1;
	n=(*circ).n;
	f.p=NULL;
	f.w=0;

	k=0;
	ch1=processChar(str[k++]);



	if(ch1=='E'||(ch1=='.')) return(f);
	else
	{
		(*circ).ngates++;

		if(ch1=='V'||ch1=='P'||ch1=='R') ch2=str[k++]; // get gate subtype designation

		// m is number of gate lines
		if(ch1=='N') m=1;
		else if(ch1=='C'||ch1=='V') m=2;
		else if(ch1=='P') m=3;
		else { // read number of lines
			ch3=str[k++];
			m=ch3-'0';
			ch3=str[k++];
			while(ch3>='0'&&ch3<='9')
			{
				m=m*10+ch3-'0';
				ch3=str[k++];
			}
		}

		// for R or Q gate get divisor
		if(ch1=='R'||ch1=='Q')
		{
			sign = 1;
			if(ch3!=':')
			{
				printf("Error when reading line: %s", str);
				printf("ch1: %c, ch2: %c, ch3: %c",ch1, ch2, ch3);
				throwException("error in R/Q gate spec (missing or misplaced :)\n",0);
			}
			ch3=str[k++];
			if(ch3=='-') {sign = -1; ch3=str[k++]; }
			div=0;
			while(ch3>='0'&&ch3<='9')
			{
				div=div*10+ch3-'0';
				ch3=str[k++];
			}
			div *= sign;
			if(VERBOSE) printf("divisor of R/Q gate: %d/%d\n", m,div);
		}


		if (m > circ->n && m!=2) { // too much lines, more than supported by the circuit
			printf("Error when reading line: %s", str);
			printf("Too much lines (%d)! Circuit only supports %d lines.", m, circ->n);
			f.p=NULL; return(f);
		}
		// define line controls
		for(i=0;i<n;i++) line[i]=-1;
		for(i=0;i<m-1;i++)
		{;
			k=k+getstr(&str[k],token);
			j=getlabel(token,*circ,&cont);
			if (j==-1) {f.p=NULL; return(f);}
			if (cont == -1) {f.p=NULL; printf("Too few variables."); return(f);}
			line[j]=cont; // control line  NOTE embedded assignment to j
			pc[i]=j;
		}
		k=k+getstr(&str[k],token);
		t=getlabel(token,*circ,&cont);
		if (t==-1) {f.p=NULL; return(f);}
		line[t]=2;  // target line  NOTE embedded assignment to t

		// set f to point to QMDD for gate
		if(ch1=='T'||ch1=='C'||ch1=='N') // T, C or N gate
		{
			if(m==1||ch1=='N') circ->ngate=1;
			else if(m==2||ch1=='C') circ->cgate=1;
			else circ->tgate=1;
			f=TTlookup(n,m,t,line);
			if(f.p==NULL)
			{
				f=QMDDmvlgate(Nm,n,line);
				TTinsert(n,m,t,line,f);
			}
			(*circ).qcost+=gate_qcost(m,n,TOFFOLI_GATE);
		}
		else if(ch1=='F') // Fredkin gate
		{
			circ->fgate=1;
			f=QMDDmvlgate(Nm,n,line);
			for(i=0;i<n;i++) line[i]=-1;
			line[t]=1;
			line[pc[m-2]]=2;
			f2=QMDDmvlgate(Nm,n,line);
			f=QMDDmultiply(f2,QMDDmultiply(f,f2));
			(*circ).qcost+=gate_qcost(m,n,FREDKIN_GATE);
		}
		else if(ch1=='P') // Peres gate
		{
			circ->pgate=1;
			f=QMDDmvlgate(Nm,n,line);
			line[t]=-1;
			line[pc[1]]=2;
			f2=QMDDmvlgate(Nm,n,line);
			if(ch2==' ') f=QMDDmultiply(f2,f);
			else if(ch2=='I') f=QMDDmultiply(f,f2);
			else printf("invalid subtype for Peres gate\n");
			(*circ).qcost+=4;  // fixed cost for PERES gate
		}
		else if(ch1=='H') f=QMDDmvlgate(Hm,n,line); // Hadamard gate
		else if(ch1=='Z') f=QMDDmvlgate(Zm,n,line); // Pauli-Z gate
		else if(ch1=='S') f=QMDDmvlgate(Sm,n,line); // Phase gate
		else if(ch1=='0') f=QMDDmvlgate(ZEROm,n,line); // zero pseudo gate
		else if(ch1=='V')	// V or V+ gate
		{
			circ->vgate=1;
			if(ch2==' ') f=QMDDmvlgate(Vm,n,line);
			else if(ch2=='P'||ch2=='+') f=QMDDmvlgate(VPm,n,line);
			else {
				printf("invalid V subtype  '%c'\n",ch2);
				throwException("",0);
			}
			(*circ).qcost+=1;
		}
		else if(ch1=='Q')
		{
			Qm[1][1]=Cmake(QMDDcos(1,div),QMDDsin(1,div));
			f=QMDDmvlgate(Qm,n,line);
		}
		else if(ch1=='R') // Rotation gate
		{
			div *= 2;
			if(ch2=='X')
			{
				Rm[0][0]=Rm[1][1]=Cmake(QMDDcos(m,div),Dzero);
				Rm[0][1]=Rm[1][0]=Cmake(Dzero,QMDDsin(m,-div));
			} else if(ch2=='Y')
			{
				Rm[0][0]=Rm[1][1]=Cmake(QMDDcos(m,div),Dzero);
				Rm[0][1]=Cmake(QMDDsin(m,-div),Dzero);
				Rm[1][0]=Cmake(QMDDsin(m,div),Dzero);
			} else if(ch2=='Z')
			{
				Rm[0][0]=Cmake(QMDDcos(m,div),QMDDsin(m,-div));
				Rm[0][1]=Rm[1][0]=Cmake(Dzero,Dzero);
				Rm[1][1]=Cmake(QMDDcos(m,div),QMDDsin(m,div));
			} else {
				printf("invalid rotation type  '%c'\n",ch2);
				throwException("",0);
			}
			f=QMDDmvlgate(Rm,n,line); // do rotation gate
		} else {
			printf("invalid gate type  '%c'\n",ch1);
			throwException("",0);
		}
	}
return(f);
}


QMDDedge QMDDreadGate(FILE *infile,QMDDrevlibDescription *circ)
{

  char lineFromInfile[256]; // read one line from the FILE
  char ch1;

  ch1=getch(infile);
  // skip characters until the first meaningful is read
  while(ch1==' '||ch1=='\n'||ch1=='#') {
	  if (ch1=='#')
		  skip2eol(infile);
	  ch1=getch(infile);
  }
  lineFromInfile[0] = ch1;
  getline(infile,&lineFromInfile[1]);

  //printf("readGateFromString: %s\n",lineFromInfile);
  
  return QMDDreadGateFromString(lineFromInfile,circ);

}

/*******************************************************************************/

QMDDrevlibDescription QMDDrevlibHeader(FILE *infile)
{
  int cont,header,n,p,i,j,k;
  char cmd[MAXSTRLEN],tLabel[MAXSTRLEN],tInput[MAXSTRLEN];
  char ch,ch1;
  CircuitLine temp;
  QMDDrevlibDescription circ;
  
  circ.nancillary=circ.ngarbage=0;
  header=1;
  if(VERBOSE) printf("Reading header");
  while(header)
  {
    ch=getch(infile);
    if(VERBOSE) printf("%c",ch);
    if(ch=='#') skip2eol(infile);
    else {
      while(ch==' '||ch=='\n') ch=getch(infile);
      if(ch!='.')
      {
        printf("invalid file\n");
	    circ.n=0;
	    return(circ);
      }
      getstr(infile,cmd);
      if(0==strcmp(cmd,"BEGIN")) header=0;    // end of header information
      else if(0==strcmp(cmd,"VERSION"))
      {
        ch1=getch(infile);
        while(ch1==' ') ch1=getch(infile);
        p=0;
        while(ch1!='\n')
        {
          circ.version[p++]=ch1;
          ch1=getch(infile);
        }
        circ.version[p]=0;
      } else if(0==strcmp(cmd,"NUMVARS"))
      {
        circ.n=n=getint(infile);
        if(VERBOSE) printf("\nnumber of variables %d\n",n);
      } else if(0==strcmp(cmd,"VARIABLES"))
      {
        for(p=n-1;p>=0;p--)
        {
          getstr(infile,circ.line[p].variable);
          strcpy(circ.line[p].input,circ.line[p].variable);				// set inputs to variable names
          strcpy(circ.line[p].input,circ.line[p].variable);				// set output to variable names
          circ.line[p].ancillary=circ.line[p].garbage='-';  // by default no lines are ancillaries or garbage
        }
        Nlabel=n;
      } else if(0==strcmp(cmd,"INPUTS"))
      {
        for(p=n-1;p>=0;p--)
        {
          getstr(infile,circ.line[p].input);
        }
      } else if(0==strcmp(cmd,"OUTPUTS"))
      {
        for(p=n-1;p>=0;p--)
        {
          getstr(infile,circ.line[p].output);
        }
      } else if(0==strcmp(cmd,"CONSTANTS"))
      {
        for(p=n-1;p>=0;p--)
        {
          circ.line[p].ancillary=getnbch(infile);
          if(circ.line[p].ancillary!='-') circ.nancillary++;
        }
        skip2eol(infile);
      } else if(0==strcmp(cmd,"GARBAGE"))
      {
        for(p=n-1;p>=0;p--)
        {
          circ.line[p].garbage=getnbch(infile);
          if(circ.line[p].garbage!='-') circ.ngarbage++;
        }
        skip2eol(infile);
      } else if(0==strcmp(cmd,"DEFINE"))
      {
        while(strcmp(cmd,"ENDDEFINE"))
        {
          skip2eol(infile);
          ch1=getch(infile);
          getstr(infile,cmd);
        }
        //skip2eol(infile);
      }
    }
  }

  if(VERBOSE) printf("completed.\n");
  for(i=0;i<circ.n;i++) circ.inperm[i]=i;

  return(circ);
}

//////////////////////////////////////////////////////////////////////////////////////////////

QMDDrevlibDescription QMDDcircuitRevlib(char *fname,QMDDrevlibDescription firstCirc, int match)
// reads a circuit in Revlib format: http://www.revlib.org/documentation.php 
{
	FILE *infile;

	QMDDrevlibDescription circ;

	int first,i,j,k,p,t,line[MAXN],order[MAXN],invorder[MAXN];
	char ch,tOutput[MAXSTRLEN],flag[MAXN],*cp,*cq,*cr;
	CircuitLine tline;

	QMDDedge e,f,olde;


	// get name of input file, open it and attach it to file (a global)
	infile=openTextFile(fname,'r');
	if(infile == NULL) {
		firstCirc.n = 0;
		return firstCirc; // failed to open infile.
	}

	//Read header from infile
	circ=QMDDrevlibHeader(infile);


	circ.ngate=circ.cgate=circ.tgate=circ.fgate=circ.pgate=circ.vgate=0;
	circ.qcost=circ.ngates=0;


	if(match)   //  match input order of this circuit to the one in firstCirc
	{
		for(i=0;i<firstCirc.n;i++)
			if(firstCirc.line[i].ancillary=='-'&&strcmp(firstCirc.line[i].input,circ.line[i].input))
			{
				for(j=i+1;j<circ.n;j++)
					if(0==strcmp(firstCirc.line[i].input,circ.line[j].input)) break;
				if(j==circ.n) printf("error in line match\n");
				tline=circ.line[i];
				circ.line[i]=circ.line[j];
				circ.line[j]=tline;
			}
	}


	first=1;
	e = QMDDident(0,circ.n);

	while(1) // read gates
	{

		f=QMDDreadGate(infile,&circ);
		if(f.p==NULL) break;

		if(first) // first gate in circuit
		{
			first = 0;
			e=f;
			QMDDincref(e);
		}
		else // second and subsequent gates
		{
			olde=e;
			e=QMDDmultiply(f,e); // multiply QMDD for gate * QMDD for circuit to date

			QMDDincref(e);
			QMDDdecref(olde);
		}
		if(GCswitch) QMDDgarbageCollect();
		
	}

	for(i=0;i<circ.n;i++) circ.outperm[i]=i;

	skip2eof(infile); // skip rest of input file

	circ.e=e;

	i=0;
	if(circ.ngate) circ.kind[i++]='N';
	if(circ.cgate) circ.kind[i++]='C';
	if(circ.tgate) circ.kind[i++]='T';
	if(circ.fgate) circ.kind[i++]='F';
	if(circ.pgate) circ.kind[i++]='P';
	if(circ.vgate) circ.kind[i++]='V';
	circ.kind[i]=0;

	i=0;
	if(circ.nancillary>0) circ.dc[i++]='C';
	if(circ.ngarbage>0) circ.dc[i++]='G';
	circ.dc[i]=0;

	fclose(infile);
	return(circ);
}


