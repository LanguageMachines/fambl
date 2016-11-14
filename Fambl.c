/* ==========================================================================

   Fambl.c

   Fambl - Family-Based Learning
   -----------------------------

   Copyright (c) 1997-2010 Antal van den Bosch
   ILK Research Group / Tilburg University
   Written by Antal van den Bosch / Antal.vdnBosch@uvt.nl

   Fambl is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   paramsearch is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, see <http://www.gnu.org/licenses/>.

   kickstart:

   [1] > make
   [2] > Fambl -h

   ====================================================================== */

#include<stdio.h>
#include<string.h>
#include<stddef.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<limits.h>
#include<getopt.h>

#define MAIN 1
#include "Common.h"
#undef MAIN

#include "Classify.h"
#include "Read.h"
#include "Metrics.h"
#include "Family.h"

/* main - carries Fambl through 1 complete experiment. A concededly
   unelegant procedure. Will be cleaned up, cut up, and generalised
   some day.
*/
int main(int argc, char **argv)
{
  /* all the functions that get called from main are here */
  void intro(void);
  void empty_memory(void);
  void bye(void);
  void timer(void);
  void set_default_values(void);
  void display_usage(void);

  /* some variables that are used for option checking */
  time_t totalbegintime,totalendtime;
  char   *rootname=NULL;
  char   willtest=0;
  char   readstring[1024];
  FILE   *test;
  int    optcount=0,getoptint=0,readint=0,learntime,testtime;

  char famblname[NAMELEN];
  char testname[NAMELEN];
  char dataname[NAMELEN];
  char cpsname[NAMELEN];

  setbuf(stdout,NULL);

  /* a welcome to the user; the current time; and a time stamp
     to measure total running time */
  intro();
  VERB=1;
  timer();
  time(&totalbegintime);

  /* set all changable global things to their default value */
  set_default_values();

  /* read all selected options and their arguments from the
     commandline */

  while ((getoptint=getopt(argc,argv,"f:t:u:bg:G:d:oneL:DXC:am:w:Wk:K:v:pijcr:hx:y:"))!=EOF)
  {
    switch (getoptint)
    {
      /* f: training file */
      case 'f':
	rootname=optarg;
      break;

      /* t: test file */
      case 't':
	strcpy(testname,optarg);
	willtest=1;
      break;

      /* o: write binary fambl file */
      case 'o':
	WRITEFAM=1;
	if (VERB)
	  fprintf(stderr,"     will write binary FAMBL file\n");
      break;

      /* e: read binary fambl file */
      case 'e':
	READFAM=1;
	if (VERB)
	  fprintf(stderr,"     will read binary FAMBL file\n");
      break;

      /* a: weight adaptation */
      case 'a':
	ALGORITHM=1;
	if (VERB)
	  fprintf(stderr,"     unpacking multi-valued into binary features\n");
      break;

      /* w: exemplar weighting */
      case 'w':
	sscanf(optarg,"%d",&readint);
	if ((readint<0)||(readint>2))
	  {
	    fprintf(stderr," <!> illegal value of -w parameter - must be in [0-2]\n\n");
	    exit(1);
	  }
	EWEIGHT=readint;
      break;

      /* W: write cps file */
      case 'W':
	PRINTCPS=1;
	if (VERB)
	  fprintf(stderr,"     CPS values will be printed to file\n");
      break;

      /* g: feature weighting */
      case 'g':
	sscanf(optarg,"%s",readstring);
	readint=atoi(readstring);
	if (strlen(readstring)==1)
	  {
	    if ((readint<0)||(readint>5))
	      { fprintf(stderr," <!> illegal value of -g parameter - must be in [0-5]\n\n");
	      exit(1);
	      }
	    FWEIGHT=readint;
	  }
	else
	  {
	    test=fopen(readstring,"r");
	    if (test==NULL)
	      {
		fprintf(stderr," <!> %s: no such weight file.\n\n",
			readstring);
		exit(1);
	      }
	    fclose(test);
	    FWEIGHT=8;
	  }
      break;

      /* G: read in binary feature weights */
      case 'G':
	sscanf(optarg,"%s",readstring);
	test=fopen(readstring,"r");
	if (test==NULL)
	  {
	    fprintf(stderr," <!> %s: no such weight file.\n\n",
		    readstring);
	    exit(1);
	  }
	fclose(test);
	FWEIGHT=9;
      break;

      /* d: distance weighting */
      case 'd':
	sscanf(optarg,"%d",&readint);
	if ((readint<0)||(readint>3))
	{ fprintf(stderr," <!> illegal value of -d parameter - must be in [0-3]\n\n");
	  exit(1);
        }
	DISTANCE=readint;
      break;

      /* m: value distance metric */
      case 'm':
	sscanf(optarg,"%d",&readint);
	if ((readint<0)||(readint>2))
	  {
	    fprintf(stderr," <!> illegal value of -m parameter - must be in [0-2]\n\n");
	    exit(1);
	  }
	METRIC=readint;
      break;

      /* L: MVDM threshold */
      case 'L':
	sscanf(optarg,"%d",&readint);
	if (readint<=0)
	{ fprintf(stderr," <!> illegal value of -L parameter - must be greater than 0\n\n");
	  exit(1);
        }
	MVDML=readint;
      break;

      /* x: collapse-disjunct-to-wildcard threshold */
      case 'x':
	sscanf(optarg,"%d",&readint);
	if (readint<=1)
	{ fprintf(stderr," <!> illegal value of -x parameter - must be greater than 1\n\n");
	  exit(1);
        }
	WILDCARDTHRESHOLD=readint;
	if (VERB)
	  fprintf(stderr,"     threshold of # values to collapse disjuncts into wildcards set to %d\n",
		 WILDCARDTHRESHOLD);
      break;

      /* u: feature type */
      case 'u':
	strcpy(feattype,optarg);
	if (strstr(feattype,"n")!=NULL)
	  numerics=1;
      break;

      /* n: numerics trigger */
      case 'n':
	numerics=1;
	if (VERB)
	  fprintf(stderr,"     all feature values are numeric\n");
      break;

      /* b: use plain IB1 */
      case 'b':
	FAMKA=0;
	REPORTRATE=10000;
	if (VERB)
	  fprintf(stderr,"     performing plain MBL (IB1); no family extraction\n");
      break;

      /* D: print NN distance along with classification output */
      case 'D':
	distprint=1;
	if (VERB)
	  fprintf(stderr,"     printing NN distance along with classification output\n");
      break;

      /* X: allow instances to be in several families */
      case 'X':
	inclusive=1;
	if (VERB)
	  fprintf(stderr,"     instances are allowed to be in several families\n");
      break;

      /* p: write classification accuracy to .% file*/
      case 'p':
	WRITEPERC=1;
	if (VERB)
	  fprintf(stderr,"     writing classification accuracy to separate .%% file\n");
      break;

      /* c: use centroids instead of hyperrectangles */
      case 'c':
	centroids=1;
	if (VERB)
	  fprintf(stderr,"     creating centroids instead of hyperrectangles\n");
      break;

      /* k: give explicit k for the k-NN engine */
      case 'k':
	sscanf(optarg,"%d",&readint);
	if ((readint<1)||(readint>MAXPATWIDTH))
	{ fprintf(stderr," <!> illegal value of -k parameter - must be in [1-%d]\n\n",MAXPATWIDTH);
	  exit(1);
        }
	FAMBLKA=readint;
	if (VERB)
	  fprintf(stderr,"     k in k-NN set to %d\n",
		 FAMBLKA);
      break;

      /* K: give explicit k for family extraction */
      case 'K':
	sscanf(optarg,"%d",&readint);
	if ((readint<0)||(readint>MAXPATWIDTH))
	{ fprintf(stderr," <!> illegal value of -K parameter - must be in [0-%d]\n\n",MAXPATWIDTH);
	  exit(1);
        }
	FAMKA=readint;
	if (VERB)
	  fprintf(stderr,"     K in family extraction set to %d\n",
		 FAMKA);
      break;

      /* C: give explicit cardinality ceiling for feature combinations */
      case 'C':
	sscanf(optarg,"%d",&readint);
	if ((readint<1)||(readint>MAXCOMBINATIONS))
	{ fprintf(stderr," <!> illegal value of -k parameter - must be in [1-%d]\n\n",MAXCOMBINATIONS);
	  exit(1);
        }
	NRCOMBINATIONS=readint;
	if (VERB)
	  fprintf(stderr,"     feature combination cardinality ceiling set to %d\n",
		 NRCOMBINATIONS);
      break;

      /* v: set verbosity level */
      case 'v':
	sscanf(optarg,"%d",&readint);
	if ((readint<0)||(readint>4))
	{ fprintf(stderr," <!> illegal value of verbosity level - must be in [0-4]\n\n");
	  exit(1);
        }
	else
	{ /* set verbosity level in an inefficient and silly way */
	  if (!readint)
	    VERB=0;
	  if (readint==1)
	    VERB=1;
	  if (readint==2)
	  { VERB=1;
	    VERB2=1;
	  }
	  if (readint==3)
	  { VERB=1;
	    VERB2=1;
	    VERB3=1;
	  }
	  if (readint==4)
	  { VERB=1;
	    VERB2=1;
	    VERB3=1;
	    VERB4=1;
	  }
	}
      break;

      /* i: unused option
      case 'i':
	if (VERB)
	  fprintf(stderr,"     unused option\n");
      break; */

      /* j: use joker (wildcard) values instead of value disjunctions */
      case 'j':
	WILDCARD=1;
	if (VERB)
	  fprintf(stderr,"     use joker (wildcard) values instead of value disjunctions\n");
      break;

      /* r: report rate */
      case 'r':
	sscanf(optarg,"%d",&readint);
	if (readint<1)
	{ fprintf(stderr," <!> illegal value of -r parameter - must be larger than zero\n\n");
	  exit(1);
        }
	REPORTRATE=readint;
      break;

      /* h: display usage */
      case 'h':
        display_usage();
	exit(1);
      break;

      /* unknown option */
      case '?':
	fprintf(stderr," <!> run Fambl -h for usage info.\n\n");
        exit(1);
    }
    optcount++;
  }

  /* now a lot of options need to be re-checked. If something seems
     truly wrong, bail out. */

  /* with automatic feature selection, set the sweeping functions
     accordingly */
  if (ALGORITHM)
    {
      sweep=sweep_nomvdm_fs;
      sweep_c=sweep_c_nomvdm_fs;
    }

  /* in case no options are found before the end of the commandline,
     or a strange string is found, complain and exit. */
  if (!optcount)
    {
      fprintf(stderr," <!> missing information on the commandline.\n");
      fprintf(stderr,"     run Fambl -h for usage info.\n\n");
      exit(1);
    }

  /* obligatory information: the -f statement. */
  if (rootname==NULL)
    {
      fprintf(stderr," <!> compulsory option -f missing.\n");
      fprintf(stderr,"     run Fambl -h for more usage info.\n\n");
      exit(1);
    }

  /* check for the datafile. If not there, complain and exit. */
  strcpy(dataname,rootname);
  if (fopen(dataname,"r")==NULL)
    {
      fprintf(stderr," <!> no data file %s\n\n",
	      dataname);
      exit(1);
    }

  /* initialize the name of the CPS file */
  if (PRINTCPS)
    {
      sprintf(cpsname,"%s.cps",
	      dataname);
    }

  /* if the -t option was selected, check for the test file. If it's not
     there, complain and exit. When it was not selected, the data file
     is also used for testing. */
  if (willtest)
    {
      if (fopen(testname,"r")==NULL)
	{
	  fprintf(stderr," <!> no test file %s\n",
		  testname);
	  exit(1);
	}
    }
  else
    {
      strcpy(testname,rootname);
      fprintf(stderr,"     using data file %s for testing\n",
	     testname);
  }

  /* prepare the name of the Fambl file. If an existing file should be
     read, check for its existence. If it's not there, complain and
     exit. */
  strcpy(famblname,rootname);
  strcat(famblname,".fambl");
  if (READFAM)
    {
      if (fopen(famblname,"r")==NULL)
	{
	  fprintf(stderr," <!> no Fambl file %s\n",
		  famblname);
	  exit(1);
	}
      else
	if (VERB2)
	  fprintf(stderr,"     found Fambl file %s\n",
		 famblname);
    }

  /* if the seed is not user-set, just take the time as seed. Otherwise,
     seed with the user's choice. */
  if (!SEED)
    srand48((unsigned long int)time(NULL));
  else
    srand48(SEED);

  /* tell about the chosen metric. */
  if (VERB)
  {
    fprintf(stderr,"     metric scheme set to ");

    if (FWEIGHT==0)
      fprintf(stderr,"no feature wgts, ");
    if (FWEIGHT==1)
      fprintf(stderr,"IG feature wgts, ");
    if (FWEIGHT==2)
      fprintf(stderr,"GR feature wgts, ");
    if (FWEIGHT==3)
      fprintf(stderr,"chi square feature wgts, ");
    if (FWEIGHT==4)
      fprintf(stderr,"shared variance feature wgts, ");
    if (FWEIGHT==5)
      fprintf(stderr,"log-likelihood feature wgts, ");
    if (FWEIGHT==6)
      fprintf(stderr,"GR feature wgts with IGTree-ish boosting, ");
    if (FWEIGHT==7)
      fprintf(stderr,"GR feature wgts with TRIBL-ish boosting, ");
    if (FWEIGHT==8)
      fprintf(stderr,"user-specified wghts from file, ");

    if (DISTANCE==0)
      fprintf(stderr,"no distance wgts, ");
    if (DISTANCE==1)
      fprintf(stderr,"linear-inverse distance wgts, ");
    if (DISTANCE==2)
      fprintf(stderr,"inverse distance wgts, ");
    if (DISTANCE==3)
      fprintf(stderr,"exponential-decay distance wgts, ");

    if (METRIC==0)
      fprintf(stderr,"no MVDM\n");
    if (METRIC==1)
      fprintf(stderr,"MVDM\n");
    if (METRIC==2)
      fprintf(stderr,"Jeffrey divergence\n");

    if (MVDML>1)
      fprintf(stderr,"     MVDM frequency threshold set to non-default %d\n",
	     MVDML);

  }

  /* until now only i/o checks have been made. Let's start things
     up, slowly. */

  /* mark the time; Fambl kicks off. */
  time(&begintime);

  /* read in the data. It will be used to count types, compute
     names and weights, and will form the
     basis for family probing and extraction. */
  readin_data( dataname );

  /* count value-class occurrences. These will form the basis
     for computing gain ratio values and vdm matrices later on.
     moreover, if low-frequent values are to be suppressed,
     this is also done here. */
  count_vcocc();

  /* if no gain ratio feature weighting is used, just set the
     weights to 1.0 and perform some precomputations
     that will optimize the k-NN part later on. If gain ratio
     is used, either read in the weights or compute them.
     Then, if a weight adaptation method was selected (IGTree
     or TRIBL), adapt the weight and do the precomputation all
     over again. */
  if ((FWEIGHT==0)&&(ALGORITHM==0))
    {
      equalize_gr();
      optimize_with_gr();
      compute_grsum();
    }
  else
    { /* if (ALGORITHM==0) */
      {
	if ((FWEIGHT==3)||(FWEIGHT==4))
	  compute_chisquare();
	if (FWEIGHT==5)
	  compute_logl();
	if ((FWEIGHT==1)||
	    (FWEIGHT==2)||
	    (FWEIGHT==6)||
	    (FWEIGHT==7))
	  compute_gr();
	if (FWEIGHT==8)
	  readin_weights(readstring);
	optimize_with_gr();
	compute_grsum();
	if ((FWEIGHT==6)||(FWEIGHT==7))
	  {
	    adapt_weights_to_algorithm();
	    optimize_with_gr();
	    compute_grsum();
	  }
      }
      if (ALGORITHM)
	{
	  if (FWEIGHT==9)
	    readin_binary_weights(readstring);
	  else
	    atomic_rooster();
	}
    }

  /* when using feature weighting, but even without using it, it is wise
     do a sort of the data on the most important feature, and compute an
     index on that feature. only when numeric values are used, forget about
     presorting. */
  if (!((READFAM)||(numerics)))
  { presort_data( dataname );
    count_vcocc();
  }

  /* if MVDM is used, compute the information needed: either the
     value-class cooccurrence tables or the full MVDM matrices. */
  if (METRIC>0)
  {
    if (!ALGORITHM)
      {
	compute_mvdm();
	sweep=sweep_mvdm;
	sweep_c=sweep_c_mvdm;
	optimize_with_gr();
	compute_grsum();
      }
  }

  time(&endtime);
  learntime=(int) endtime - (int) begintime;
  if (VERB)
    {
      fprintf(stderr,"     took %d second",
	      learntime);
      if (learntime==1)
	fprintf(stderr,"\n");
      else
	fprintf(stderr,"s\n");
    }

  /* if there is no Fambl file to read, start building one. If there
     is one, read it. */
  if (!READFAM)
  {
    /* perform the family extraction stage. */
    time(&begintime);
    make_families(famblname);
    time(&endtime);
    learntime+=(int) endtime - (int) begintime;
    if (VERB)
      {
	fprintf(stderr,"     took %d second",
		(int) endtime - (int) begintime);
	  if (((int) endtime - (int) begintime)==1)
	    fprintf(stderr,"\n");
	  else
	    fprintf(stderr,"s\n");
      }
    empty_memory();
  }
  else
  {
    empty_memory();
    readin_fambl( famblname );
  }

  /* if exemplar weights are to be used, compute them here */
  if (EWEIGHT>0)
    {
      time(&begintime);
      compute_exemplar_weights( dataname, cpsname );
      time(&endtime);
      learntime+=(int) endtime - (int) begintime;
      if (VERB)
	{
	  fprintf(stderr,"     took %d second",
		  (int) endtime - (int) begintime);
	  if (((int) endtime - (int) begintime)==1)
	    fprintf(stderr,"\n");
	  else
	    fprintf(stderr,"s\n");
	}
  }

  /* classify all instances in the file designated as test file */
  time(&begintime);
  perform_test(testname);
  time(&endtime);
  testtime=(int) endtime - (int) begintime;
  if (VERB)
    {
      fprintf(stderr,"     took %d second",
	      testtime);
      if (testtime==1)
	fprintf(stderr,"\n");
      else
	fprintf(stderr,"s\n");
      if (testtime>0)
	fprintf(stderr,"     (%.2f instances per second)\n",
		((1.*totpat)/
		 (1.*testtime)));
      else
	fprintf(stderr,"     (all %d instances within one second)\n",
		totpat);
    }

  /* what's the time? */
  if (VERB)
    {
      fprintf(stderr," <*> ");
      timer();
    }

  /* time stamp: we're ready. report how much we took, and quit. */
  time(&totalendtime);
  if (VERB)
    {
      fprintf(stderr,"     Fambl spent a total of %d seconds running;\n",
	      learntime+testtime);
      fprintf(stderr,"     %d on learning, %d on testing.\n",
	      learntime,testtime);
    }
  bye();

  return 0;
}

/* set_default_values - values in the absence of commandline overrules,
   and some initialisations of booleans and dynamic functions
*/
void set_default_values(void)
{
  /* default values of user-changable parameters */

  VERB=1;
  VERB2=VERB3=VERB4=0;
  FAMBLKA=1;
  DISTANCE=0;
  FAMKA=1000;
  READFAM=0;
  WRITEFAM=0;
  METRIC=1;
  MVDML=1;
  WILDCARD=0;
  WILDCARDVALUE=-2;
  WILDCARDTHRESHOLD=1000;
  REPORTRATE=1000;
  ALGORITHM=0;           /* no value weighting    */
  EWEIGHT=0;             /* no exemplar weighting */
  PRINTCPS=0;
  FWEIGHT=2;             /* shared variance       */
  NRCOMBINATIONS=1;
  WRITEPERC=0;
  SEED=0;
  mostfreqclass=0;

  /* booleans */
  notshowntwice=1;
  superfeat=-1;
  thereisasuperfeat=0;
  alreadymadeone=0;
  numerics=0;
  nrnumerics=0;
  centroids=0;
  distprint=0;
  inclusive=0;

  /* floats */
  returndist=0.0;

  /* dynamic functions */
  sweep=sweep_nomvdm;
  sweep_c=sweep_c_nomvdm;
  getthisrank=getthisrank_normal;

  /* special inits */
  DISCRETIZER=(((int) pow(2,(8*sizeof(value))))/2)-1;
  strcpy(feattype,"");
}

/* empty_memory - clear up pattern and class buffers when family
   creation is over
*/
void empty_memory(void)
{
  int i;

  for (i=0; i<PATWIDTH; i++)
    free(feature[i]);
  free(klass);
}

/* timer - understandable local time
*/
void timer(void)
{
  struct tm *curtime;
  time_t    bintime;

  time(&bintime);
  curtime = localtime(&bintime);
  if (VERB)
    fprintf(stderr,"current time: %s", asctime(curtime));
}

/* bye!
 */
void bye(void){
  for ( int i=0; i < PATWIDTH; ++i )
    {
      for ( int j=0; j<NRCLASSES; j++)
	{
	  free(valocc[i][j]);
	}
      free(valocc[i]);
    }
  for ( int i=0; i < PATWIDTH; ++i )
    {
      free( totalocc[i] );
    }
  free( totalocc );
  if (VERB)
    fprintf(stderr,"     Fambl ready.\n\n");
}

/* display_usage - show all commandline options with default values.
   also in file USAGE.
*/
void display_usage(void)
{
  fprintf(stderr,"usage: Fambl -f training file [-t test file] [-g metric] [-m metric] [-a]\n");
  fprintf(stderr,"             [-w function] [-b] [-c] [-D] [-k knn] [-t type] [-n] [-j]\n");
  fprintf(stderr,"             [-v verb] [-o] [-e] [-z] [-r reportr.] [-x seed] [-C cardin]\n\n");

  fprintf(stderr," obligatory:\n");
  fprintf(stderr,"  -f name      training file name\n");


  fprintf(stderr," custom (off when not selected, or default when not specified):\n");
  fprintf(stderr,"  -t name      test file name\n");
  fprintf(stderr,"  -g metric    set feature (value) weighting metric:\n");
  fprintf(stderr,"               0: no feature weighting\n");
  fprintf(stderr,"               1: information gain\n");
  fprintf(stderr,"               2: information gain ratio (default)\n");
  fprintf(stderr,"               3: chi square\n");
  fprintf(stderr,"               4: shared variance\n");
  fprintf(stderr,"               5: log-likelihood\n");
  /* hidden options :-)
  fprintf(stderr,"               <file>: read weights from file\n");
  fprintf(stderr,"               6: GR+IGTree-ish weight boosting\n");
  fprintf(stderr,"               7: GR+TRIBL-ish weight boosting\n");
  fprintf(stderr,"  -G <file>    read binary weights from <file>\n"); */
  fprintf(stderr,"  -m metric    set value difference metric:\n");
  fprintf(stderr,"               0: overlap (default)\n");
  fprintf(stderr,"               1: MVDM (modified value difference metric)\n");
  fprintf(stderr,"               2: JD (Jeffrey divergence)\n");
  fprintf(stderr,"  -d metric    set distance weighting metric:\n");
  fprintf(stderr,"               0: no distance weighting (default)\n");
  fprintf(stderr,"               1: linear-inversed distance weighting\n");
  fprintf(stderr,"               2: inversed distance weighting\n");
  fprintf(stderr,"               3: exponential-decay distance weighting\n");
  fprintf(stderr,"  -a           unpack multi-valued into binary (atomic) features\n");
  fprintf(stderr,"  -C cardin.   feature combination cardinality ceiling (default 1)\n");
  fprintf(stderr,"  -w function  set family weighting metric:\n");
  fprintf(stderr,"               0: no family weighting (default)\n");
  fprintf(stderr,"               1: class prediction strength\n");
  fprintf(stderr,"               2: class prediction strength with Laplace correction\n");
  fprintf(stderr,"  -W           write <datafile>.cps file (only with -w)\n");
  fprintf(stderr,"  -b           perform basic MBL (IB1), i.e., create only single-instance\n");
  fprintf(stderr,"               families\n");
  fprintf(stderr,"  -k knn       set the k in k-NN classification (default 1)\n");
  fprintf(stderr,"  -K fam-knn   set the k in family extraction (default 1000)\n");
  fprintf(stderr,"  -L thresh.   set the MVDM frequency threshold (default 1)\n");
  fprintf(stderr,"  -u type      sets feature type according to string \"type\", e.g. nnsss sets\n");
  fprintf(stderr,"               1 and 2 numeric, and 3, 4, 5 symbolic (default all s)\n");
  fprintf(stderr,"  -n           declare all feature values as numeric (shortcut for -t nnnn...)\n");
  fprintf(stderr,"  -j           use joker (wildcard) values instead of value disjunctions\n");
  fprintf(stderr,"  -x thresh.   collapse disjuncts to wildcards above number of values\n");
  fprintf(stderr,"  -X           allow instances to be in more than one family\n");
  fprintf(stderr,"  -c           collapse hyperrectangles to into centroids (numeric features)\n");
  fprintf(stderr,"  -D           print NN distances along with classification output\n");
  fprintf(stderr,"  -v verb.     set verbosity level [0-4] (default 1)\n");
  fprintf(stderr,"  -p           write a .%% file containing classification accuracy\n");
  fprintf(stderr,"  -o           write binary Fambl file filestem.fambl\n");
  fprintf(stderr,"  -e           read binary Fambl file filestem.fambl\n");
  fprintf(stderr,"  -r reportr.  rate at which to report on progress in processing families or\n");
  fprintf(stderr,"               test instances (default 1000)\n\n");

}

/* intro - name dropping plus version number (which gets hand-coded
   each time a version is updated, just as it is updated in the README
   file: in principle, version and date mentioned here should match the
   ones in README) */
void intro(void)
{
  fprintf(stderr,"\n");
  fprintf(stderr,"Fambl (Family-based learning), version 2.3.4, 20 June 2010\n");
  fprintf(stderr,"(c) 1997-2010 ILK Research Group, Tilburg University\n");
  fprintf(stderr,"http://ilk.uvt.nl / antalb@uvt.nl\n");
}
