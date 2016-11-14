/* Common.h

   Copyright 1997-2004 ILK / Tilburg University
   Written by Antal van den Bosch
*/

#ifndef COMMON_H
#define COMMON_H
/* rerouting declaration of global variables */

#ifdef MAIN
#define EXTERN
#else
#define EXTERN extern
#endif

/* (max/def) quantity defines: */

#define MAXNREXP     1500000  /* maximal # expressions in FAMBL file        */
#define MAXEXPRESSION 131072  /* maximal expression length                  */
#define MAXPATWIDTH     2048  /* maximal pattern width (# features)         */
#define MAXVALLEN        256  /* maximal string length of values            */
#define MAXNRCLASSES    8192  /* maximal number of different classes        */
#define MAXNRVALS      32768  /* maximal number of values per feature       */
#define TOPMVDM          100  /* nr of values per feat with prestored mvdm  */
#define MAXBONUSES    131072  /* top nr of atomic features to remember      */
#define BINS              20  /* # discretization bins for numeric features */
#define MAXCOMBINATIONS   32  /* maximal number of value combinations       */
#define ALMOSTNOTHING 0.00001 /* used in computing float inequality         */
#define NAMELEN          100  /* length for filenames                       */
#define MAXFAMKA        1024  /* maximal family K                           */

/* handy defines: */
#define min(x, y) ((x)<=(y) ? (x) : (y))
#define max(x, y) ((x)>=(y) ? (x) : (y))

/* "value" as a shortname for int. */
typedef int (value);

EXTERN value  *feature[MAXPATWIDTH];
EXTERN value  *klass;
EXTERN value  *expression[MAXNREXP];
EXTERN value  *expressionclass;
EXTERN char   hit[MAXPATWIDTH];
EXTERN float  featgrvals[MAXPATWIDTH];
EXTERN float  restvals[MAXPATWIDTH];
EXTERN float  bonusrestvals[MAXPATWIDTH];
EXTERN float  *perclasssum;
EXTERN value  featorder[MAXPATWIDTH];
EXTERN char   feattype[MAXPATWIDTH];
EXTERN int    *expressionoccurrence;
EXTERN int    PATWIDTH,NRPAT,NRFAM,FAMBLKA,VERB,VERB2,VERB3,VERB4,FAMKA,
  nrtypes,nrexpressions,NRCLASSES,maxexpfound,READFAM,DEFPROBE,METRIC,
  famtop,superfeat,thisoccurrence,REPORTRATE,ALGORITHM,totpat,
  nrbonuses,EWEIGHT,WILDCARD,WILDCARDVALUE,WRITEFAM,DISCRETIZER,
  nrnumerics,FWEIGHT,nrapplicables,inst,mostfreqclass,
  NRCOMBINATIONS,DISTANCE,MVDML,WILDCARDTHRESHOLD,WRITEPERC,
  PRINTCPS;
EXTERN float  returndist,mindist;
EXTERN unsigned long int SEED;
EXTERN char   *classes[MAXNRCLASSES];
EXTERN char   *values[MAXPATWIDTH][MAXNRVALS];
EXTERN char   *selected;
EXTERN int    *group;
EXTERN int    *firstgroupmembers;
EXTERN value  *bonusselections[MAXBONUSES];
EXTERN float  *classweight[MAXBONUSES];
EXTERN float  bonuses[MAXBONUSES];
EXTERN int    orderedbonus[MAXBONUSES];
EXTERN int    *applicable_bonuses;
EXTERN float  ffmin[MAXPATWIDTH];
EXTERN float  ffmax[MAXPATWIDTH];
EXTERN int    *focalstart,*focalend;
EXTERN float  ***mvdm;
EXTERN value  **mvdmindex;
EXTERN value  **invertedmvdmindex;
EXTERN int    *nrmvdm;
EXTERN int    *mvdmbotfreq;
EXTERN int    **valocc[MAXPATWIDTH];
EXTERN int    **totalocc;
EXTERN value  thispattern[MAXPATWIDTH];
EXTERN int    nrvalues[MAXPATWIDTH];
EXTERN int    oldnrvalues[MAXPATWIDTH];
EXTERN int    combinationlimit[MAXPATWIDTH];
EXTERN int    overallclassocc[MAXNRCLASSES];
EXTERN float  maxgrsum,identical;
EXTERN char   notshowntwice,thereisasuperfeat,
              alreadymadeone,numerics,skipprobe,
              centroids,distprint,inclusive;
EXTERN int    circles;
EXTERN time_t begintime,endtime;

float  numscale(value jut);

/* dynamic functions */

EXTERN void   (*sweep)(int sweeper);
void   sweep_mvdm(int sweeper);
void   sweep_nomvdm(int sweeper);
void   sweep_nomvdm_fs(int sweeper);

EXTERN void   (*sweep_c)(int sweeper);
void   sweep_c_mvdm(int sweeper);
void   sweep_c_nomvdm(int sweeper);
void   sweep_c_nomvdm_fs(int sweeper);

EXTERN void   (*getthisrank)(void);
void   getthisrank_normal(void);
void   getthisrank_exact(void);

/* other functions in Common.c */

void   alloc_error(void);
void   my_system (char *command);
void   select_applicable_bonuses_family(int inst);
void   select_applicable_bonuses();
void   sort_bonuses(void);

#endif
