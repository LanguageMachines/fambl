/* Family.c

   Copyright 1997-2004 ILK / Tilburg University
   Written by Antal van den Bosch
*/

#include<vector>
#include<stdio.h>
#include<string.h>
#include<stddef.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<limits.h>

#include "Common.h"
#include "Family.h"

int    *currentfamily,*thisrankid;
int    nrcurrentfamily,totalfamka;
float  *thisrankdist;

char   *checked;
int    clusteredness[MAXNRCLASSES];
std::vector<int> order;
float  *ordernr;

/** make_families
*/
void make_families( const char famblname[NAMELEN] )
{
  FILE   *famblfile = 0;
  FILE   *tempsort;
  value  *tempexp,*optionalvalues;
  char   in;
  register int i,j,k;
  char   tempsortname[50];
  char   commandline[1024];
  int    currentinst=0,nroptionalvalues,totalexpressionlen,
    totalchecked,totalmembers,totalfamilies,nextfield=0;
  float  clusterednesssum;
  value  fammin,fammax;
  void   compute_instance_weights(void);
  int    familysize[1024];
  int res;
  if (VERB)
    fprintf(stderr," <*> family extraction stage started\n");

  if (WRITEFAM)
  {
    if (VERB3)
      fprintf(stderr,"     dumping expressions in file %s\n",famblname);
    famblfile=fopen(famblname,"w");
  }

  for (i=0; i<1024; i++)
    familysize[i]=0;

  allocfamily();

  if (NRPAT>MAXNREXP)
    {
      expressionclass=(value*)malloc(MAXNREXP*sizeof(value));
      expressionoccurrence=(int*)malloc(MAXNREXP*sizeof(int));
    }
  else
    {
      expressionclass=(value*)malloc(NRPAT*sizeof(value));
      expressionoccurrence=(int*)malloc(NRPAT*sizeof(int));
    }
  tempexp=(value*)malloc(MAXEXPRESSION*sizeof(value));
  optionalvalues=(value*)malloc(MAXNRVALS*sizeof(value));

  maxexpfound=0;

  for (i=0; i<NRPAT; i++)
    checked[i]=0;

  for (i=0; i<NRCLASSES; i++)
    clusteredness[i]=0;

  /* first generate an ordering of all instances; highest CPS first */
  ordernr=(float*)malloc(NRPAT*sizeof(float));
  if (ordernr==NULL)
    alloc_error();
  order.resize(NRPAT,0);

  /* when distance or circle-k are 0, make sure that family extraction
     is done the easy way. else, also sort the randomised ordering to
     ensure random selection of instances at family extraction time */
  if (FAMKA==0)
  {
    getthisrank=getthisrank_exact;
    for (i=0; i<NRPAT; i++)
      order[i]=i;
  }
  else
  {
    if (VERB2)
      {
	fprintf(stderr,"     ordering instances randomly\n");
      }
    /* randomize the training set */
    compute_instance_weights();
    if (VERB2)
    fprintf(stderr,"     sorting instances\n");
    sprintf(tempsortname,"fambls%d.tmp",
            (int) (lrand48()*((unsigned int)time(NULL))%UINT_MAX));
    tempsort=fopen(tempsortname,"w");
    for (i=0; i<NRPAT; i++)
      fprintf(tempsort,"%f %d\n",ordernr[i],i);
    fclose(tempsort);
    sprintf(commandline,"sort -rn %s -o %s -T .\n",
            tempsortname,tempsortname);
    if (VERB4)
      fprintf(stderr,"     executing commandline:\n     %s",commandline);
    my_system(commandline);
    if (VERB3)
      fprintf(stderr,"     reading numbers back from sorted file: %s\n",
	      tempsortname );
    tempsort=fopen(tempsortname,"r");
    for (i=0; i<NRPAT; i++){
       res = fscanf(tempsort,"%f %d\n",&ordernr[i],&order[i]);
       if ( res == EOF ){
	 fprintf( stderr, "premature end-of-file on %s\n", tempsortname );
	 exit(1);
       }
       if ( res != 2 ){
	 fprintf( stderr, "invalid entry in %s\n", tempsortname );
	 exit(1);
       }
    }
    fclose(tempsort);
    sprintf(commandline,"rm %s\n",tempsortname);
    my_system(commandline);
    if (VERB2)
      fprintf(stderr,"\r     sorting ready.          \n");
  }

  nrcurrentfamily=0;
  totalmembers=totalfamka=totalfamilies=totalexpressionlen=currentinst=0;

  while ((currentinst<NRPAT)&&(totalfamilies<MAXNREXP))
  {
    if (VERB3)
      fprintf(stderr,"     family # %d is built from instance %d\n",
	     totalfamilies,order[currentinst]);

    if (VERB3)
      fprintf(stderr,"\n     tracing family of instance %7d\n",
	     order[currentinst]);

    /* find the immediate family of that instance
       by looking for relatives */

    /* first, initialise */
    for (i=0; i<famtop; i++)
      thisrankid[i]=-1;
    currentfamily[0]=order[currentinst];
    checked[order[currentinst]]=1;

    /* if working with bonuses, determine the subset of applicable bonus
       selections */

    if (ALGORITHM)
      select_applicable_bonuses_family(currentfamily[0]);

    /* get a ranking of friendly neighbours qualifying for
       family membership */
    inst=currentfamily[0];
    getthisrank();

    nrcurrentfamily=1;
    k=0;
    while (thisrankid[k]!=-1)
      {
	currentfamily[nrcurrentfamily]=thisrankid[k];
	nrcurrentfamily++;
	k++;
      }

    for (i=0; i<nrcurrentfamily; i++)
      checked[currentfamily[i]]=1;

    if ((nrcurrentfamily-2)<1024)
      familysize[nrcurrentfamily-2]++;

    /* now convert it all to a family expression */
    nextfield=0;
    for (i=0; i<PATWIDTH; i++)
    {
      /* symbolic feature: make a disjunction of feature values */
      if (feattype[i]=='s')
	{
	  nroptionalvalues=0;
	  for (j=0; j<nrcurrentfamily; j++)
	    {
	      in=0;
	      for (k=0; ((!in)&&(k<nroptionalvalues)); k++)
		if (optionalvalues[k]==feature[i][currentfamily[j]])
		  in=1;
	      if (!in)
		{
		    {
		      optionalvalues[nroptionalvalues]=
			feature[i][currentfamily[j]];
		      if (VERB4)
			fprintf(stderr,"        found: %s\n",
			       values[i][optionalvalues[nroptionalvalues]]);
		      nroptionalvalues++;
		    }

		}
	    }
	  if (VERB4)
	    {
	      fprintf(stderr,"     feature %d: %d values found\n",i+1,nroptionalvalues);
	      for (j=0; j<nroptionalvalues; j++)
		{
		  fprintf(stderr,"        value %d: %s\n",
			 j+1,values[i][optionalvalues[j]]);
		}
	    }
	}
      /* numeric values: determine min and max, and define the point (single
         value) or range (multiple values) */
      else
	{
	  fammin=DISCRETIZER;
	  fammax=0;
	  for (j=0; j<nrcurrentfamily; j++)
	    {
	      fammax=max(feature[i][currentfamily[j]],fammax);
	      fammin=min(feature[i][currentfamily[j]],fammin);
	    }
	  if (fammin==fammax)
	    {
	      nroptionalvalues=1;
	      optionalvalues[0]=fammin;
	      if (VERB4)
		{
		  fprintf(stderr,"     feature %d: single value %6.4f\n",
			 i,
			 numscale(optionalvalues[0]));
		}
	    }
	  else
	    {
	      if (centroids)
		{
		  nroptionalvalues=1;
		  optionalvalues[0]=((fammin+fammax)/2);
		}
	      else
		{
		  nroptionalvalues=2;
		  optionalvalues[0]=fammin;
		  optionalvalues[1]=fammax;
		  if (VERB4)
		    {
		      fprintf(stderr,"     feature %d: value range %6.4f - %6.4f\n",
			     i,
			     numscale(optionalvalues[0]),
			     numscale(optionalvalues[1]));
		    }
		}
	    }
	}

      /* filling value and occurrence into expression */
      if (((WILDCARD)&&(nroptionalvalues>1))||
	  (nroptionalvalues>WILDCARDTHRESHOLD))
	{
	  tempexp[PATWIDTH+1+nextfield]=WILDCARDVALUE;
	  nroptionalvalues=1;
	}
      else
	{
	  for (j=0; j<nroptionalvalues; j++)
	    tempexp[PATWIDTH+1+nextfield+j]=optionalvalues[j];
	}

      /* updating index field */
      tempexp[i]=PATWIDTH+1+nextfield;
      nextfield+=nroptionalvalues;
    }
    tempexp[PATWIDTH]=PATWIDTH+1+nextfield;

    expression[totalfamilies]=
      (value*)malloc((PATWIDTH+1+nextfield)*sizeof(value));
    for (i=0; i<(PATWIDTH+1+nextfield); i++)
      expression[totalfamilies][i]=tempexp[i];
    expressionclass[totalfamilies]=klass[currentfamily[0]];
    expressionoccurrence[totalfamilies]=nrcurrentfamily;

    /* determine amount of memory used.
       first count feature index when used, and values */
    if (!WILDCARD)
      totalexpressionlen+=((PATWIDTH+1+nextfield)*sizeof(value));
    else
      totalexpressionlen+=((2*PATWIDTH)*sizeof(value));
    /* then class and occurrence counter. */
    totalexpressionlen+=(sizeof(value)+sizeof(int));

    if ((PATWIDTH+1+nextfield)>maxexpfound)
      maxexpfound=PATWIDTH+1+nextfield;

    clusteredness[klass[currentfamily[0]]]++;

    if (VERB3)
    {
      fprintf(stderr,"     found a family with class %s and %d members                      \n",classes[klass[currentfamily[0]]],nrcurrentfamily);
      for (i=0; i<nrcurrentfamily; i++)
	{
	  fprintf(stderr,"     %3d. ",i+1);
	  for (j=0; j<PATWIDTH; j++)
	    {
	      fprintf(stderr,"%s,",
		      values[j][feature[j][currentfamily[i]]]);
	    }
	  fprintf(stderr," class %d/%s\n",
		 klass[currentfamily[i]],
		 classes[klass[currentfamily[0]]]);
	}

      fprintf(stderr,"     converted into the following expression:\n     ");
      for (i=0; i<PATWIDTH; i++)
	{
	  if (feattype[i]=='s')
	    {
	      for (j=expression[totalfamilies][i];
		   j<expression[totalfamilies][i+1];
		   j++)
		{
		  if (expression[totalfamilies][j]==WILDCARDVALUE)
		    fprintf(stderr,"*,");
		  else
		    {
		      fprintf(stderr,"%s,",values[i][expression[totalfamilies][j]]);
		    }
		}
	      fprintf(stderr," ");
	    }
	  else
	    {
	      if (expression[totalfamilies][i]==expression[totalfamilies][i+1]+1)
		fprintf(stderr,"%6.4f,",
		       numscale(expression[totalfamilies][expression[totalfamilies][i]]));
	      else
		fprintf(stderr,"%6.4f-%6.4f,",
		       numscale(expression[totalfamilies][expression[totalfamilies][i]]),
		       numscale(expression[totalfamilies][expression[totalfamilies][i]+1]));
	    }
	}
      fprintf(stderr,"\n     class %d/%s, occ. %d\n\n",
	     expressionclass[totalfamilies],
	     classes[expressionclass[totalfamilies]],
	     expressionoccurrence[totalfamilies]);
    }

    if (WRITEFAM)
      {
	fwrite(expression[totalfamilies],
	       sizeof(value),
	       PATWIDTH+1+nextfield,
	       famblfile);
	fwrite(&expressionclass[totalfamilies],
	       sizeof(value),
	       1,
	       famblfile);
	fwrite(&expressionoccurrence[totalfamilies],
	       sizeof(int),
	       1,
	       famblfile);
      }

    totalmembers+=nrcurrentfamily;
    totalfamilies++;

    if ((VERB)&&((totalfamilies%REPORTRATE)==0))
      {
	totalchecked=0;
	for (i=0; i<NRPAT; i++)
	  if (checked[i])
	    totalchecked++;
	fprintf(stderr,"\r     %6d families (covering%6.2f%%), %7.2f members av.\n",
	       totalfamilies,
	       ((100.*totalchecked)/(1.*NRPAT)),
	       (1.*totalmembers)/(1.*totalfamilies));
      }

    if (currentinst<NRPAT)
      while ((currentinst<NRPAT)&&(checked[order[currentinst]]))
	currentinst++;
  }

  expressionclass=(value*)realloc(expressionclass,totalfamilies*sizeof(value));
  expressionoccurrence=(int*)realloc(expressionoccurrence,totalfamilies*sizeof(int));

  nrexpressions=totalfamilies;
  clusterednesssum=0;
  for (i=0; i<NRCLASSES; i++)
    clusterednesssum+=
      clusteredness[i]*
      ((1.*overallclassocc[i])/(1.*NRPAT));

  /* summarize everything */
  if (VERB)
    {
      if (VERB2)
	fprintf(stderr,"     family extraction stage is finished.                                   \n");
      fprintf(stderr,"     # families                        : %9d\n",totalfamilies);
      fprintf(stderr,"     average # members                 : %9.4f\n",
	      (1.*totalmembers)/(1.*totalfamilies));
      fprintf(stderr,"     average k distances               : %9.4f\n",
	      (1.*totalfamka)/(1.*totalfamilies));
      fprintf(stderr,"     average description length (bytes): %9.4f\n",
	     (1.*totalexpressionlen)/(1.*totalfamilies));
      fprintf(stderr,"     compression (raw memory)          : %9.4f %%\n",
	     100.-((100.*totalexpressionlen)/
		   ((NRPAT*PATWIDTH*sizeof(value))+
		    (NRPAT*sizeof(value)))));
      fprintf(stderr,"     compression (vs instance types)   : %9.4f %%\n",
	     100.-((100.*totalexpressionlen)/
		   ((nrtypes*PATWIDTH*sizeof(value))+
		    (nrtypes*sizeof(value))+
		    (nrtypes*sizeof(int)))));
      fprintf(stderr,"     #type vs. #fam reduction          : %9.4f %%\n",
	     100.-((100.*totalfamilies)/(1.*nrtypes)));
      fprintf(stderr,"     clusteredness                     : %9.2f\n",
	     clusterednesssum);
    }

  // family stats
  if (VERB3)
    {
      for (i=0; i<=50; i++)
	fprintf(stderr,"     family size %2d - %7d - %8.5f\n",
		i,familysize[i],(100.*familysize[i])/(1.*totalfamilies));
    }

  if (WRITEFAM)
    fclose(famblfile);

  freefamily();
  free(tempexp);
  free(optionalvalues);
}

/* allocfamily - allocate the temporary buffers needed for tracing
   families
*/
void allocfamily(void)
{
  checked=(char*)malloc(NRPAT*sizeof(char));
  if (checked==NULL)
    alloc_error();
  currentfamily=(int*)malloc(famtop*sizeof(int));
  if (currentfamily==NULL)
    alloc_error();
  thisrankdist=(float*)malloc(famtop*sizeof(float));
  if (thisrankdist==NULL)
    alloc_error();
  thisrankid=(int*)malloc(famtop*sizeof(int));
  if (thisrankid==NULL)
    alloc_error();

}

/* freefamily - free the temporary buffers needed for tracing families
*/
void freefamily(void)
{
  free(ordernr);
  free(checked);
  free(currentfamily);
  free(thisrankdist);
  free(thisrankid);
}

/* getthisrank_normal - for instance inst, get a ranking of all friendly
   neighbours within the set thresholds
*/
void getthisrank_normal(void)
{
  char   samecat,useless,kachanged=1;
  register int i,k;
  int    nrrank,thisrank,kachange,glitch=0,thisstart,thisend,max;
  float  previousdist,boundary;
  char   went_through_focal=0;

  mindist=0.0;
  identical=0.0;
  boundary=maxgrsum-ALMOSTNOTHING;

  if (VERB3)
  {
    fprintf(stderr,"\r     computing neighbour ranking of instance %7d, ",
	   inst);
    for (i=0; i<PATWIDTH; i++)
    {
      if (feattype[i]=='s')
        fprintf(stderr,"%s,",values[i][feature[i][inst]]);
      else
	fprintf(stderr,"%6.4f,",numscale(feature[i][inst]));
    }
    fprintf(stderr,"%s\n",classes[klass[inst]]);
  }

  nrrank=0;

  /* start with matching all instances with the same feature value at
     the most important feature when using feature weighting. */
  if (!numerics)
    {
      thisstart=focalstart[feature[featorder[0]][inst]];
    thisend=focalend[feature[featorder[0]][inst]];
  }
  else
  {
    thisstart=0;
    thisend=NRPAT;
    went_through_focal=1;
  }
  max=thisstart;

  while (max<NRPAT)
  {
    if ((klass[inst]!=klass[max])||
	(!checked[max])||
	(inclusive))
    {
      /* compute similarity, sets global thingies */
      sweep(max);

      if ((identical>mindist-ALMOSTNOTHING)&&
	  (circles<FAMKA))
      {
	/* insert the number in the rank */
	thisrank=0;

	samecat=1;
	if (klass[inst]!=klass[max])
	  samecat=0;

	useless=0;
	kachanged=0;
	if (nrrank>0)
	  {
	    while ((identical<thisrankdist[thisrank])&&
		   (thisrankid[thisrank]!=-1)&&
		   (thisrank<nrrank))
	      thisrank++;
	    if ((thisrankid[thisrank]==-1)&&
		(identical<thisrankdist[thisrank]))
	      useless=1;
	    if (!useless)
	      {
		if (!samecat)
		  while ((identical==thisrankdist[thisrank])&&
			 (thisrankid[thisrank]!=-1)&&
			 (thisrank<nrrank))
		    thisrank++;
		if (identical>thisrankdist[thisrank])
		  kachanged=1;
	      }
	  }
	if ((!useless)&&(thisrank<NRPAT))
	  {
	    if (thisrank<nrrank)
	      {
		for (k=nrrank; k>thisrank; k--)
		  {
		    thisrankid[k]   = thisrankid[k-1];
		    thisrankdist[k] = thisrankdist[k-1];
		  }
	      }
	    thisrankdist[thisrank]=identical;
	    if (samecat)
	      thisrankid[thisrank]=max;
	    else
	      {
		thisrankid[thisrank]=-1;
		nrrank=thisrank;
	      }
	    if (nrrank<NRPAT)
	      nrrank++;

	    if (kachanged)
	      {
		kachange=0;
		glitch=nrrank;
		previousdist=boundary;
		for (i=0; ((kachange<FAMKA)&&(i<nrrank)); i++)
		  {
		    if (thisrankdist[i]<previousdist)
		      {
			previousdist=thisrankdist[i];
			kachange++;
			if (kachange==FAMKA)
			  glitch=i;
		      }
		  }
		if (glitch<nrrank)
		  nrrank=glitch;
	      }
	    kachanged=0;
	    mindist=thisrankdist[nrrank-1];
	  }
      }
    }
    /* when the focal feature is looked through, start from the beginning.
       when it's already looked through, skip it when encountering it. */
    max++;
    if (!went_through_focal)
    {
      if (max==thisend)
      {
	max=0;
        if (max==thisstart)
	  max=thisend;
        went_through_focal=1;
      }
    }
    else
    {
      if (max==thisstart)
	max=thisend;
    }
  }

  if (VERB3)
    fprintf(stderr,"     ready ranking, %d ranked\n",
	   nrrank);
  thisrankid[nrrank]=-1;

  kachange=0;
  previousdist=thisrankdist[0];
  if (nrrank>1)
    {
      for (i=1; ((kachange<FAMKA)&&(i<nrrank)); i++)
	{
	  if (thisrankdist[i]<previousdist)
	    {
	      previousdist=thisrankdist[i];
	      kachange++;
	    }
	}
    }
  totalfamka+=kachange+1;
}

/* getthisrank_exact - for instance inst, get all clones of the same type.
   used instead of getthisrank_normal when family-similarity thresholds are
   maximally strict.
*/
void getthisrank_exact(void)
{ register int i;
  int    nrrank,max;
  char   same;

  if (VERB3)
  {
    fprintf(stderr,"\r     fetching clones of instance %7d, ",
	   inst);
    for (i=0; i<PATWIDTH; i++)
      fprintf(stderr,"%s,",values[i][feature[i][inst]]);
    fprintf(stderr,"%s\n",classes[klass[inst]]);
  }

  nrrank=0;

  /* go back and forth from the current instance; since the instances
     are sorted the clones will be surrounding it */

  /* first go back
  same=1;
  max=inst-1;
  while ((max>=0)&&(klass[max]==klass[inst])&&(same))
  {
    same=1;
    for (i=0; ((same)&&(i<PATWIDTH)); i++)
      if (feature[i][inst]!=feature[i][max])
	same=0;
    if ((same)&&(nrrank<NRPAT))
    {
      thisrankid[nrrank]=max;
      thisrankdist[nrrank]=0.0;
      nrrank++;
    }
    max--;
    } */

  /* then go forth */
  same=1;
  max=inst+1;
  while ((max<NRPAT)&&(klass[max]==klass[inst])&&(same))
  {
    same=1;
    for (i=0; ((same)&&(i<PATWIDTH)); i++)
      if (feature[i][inst]!=feature[i][max])
	same=0;
    if ((same)&&(nrrank<NRPAT))
    {
      thisrankid[nrrank]=max;
      thisrankdist[nrrank]=0.0;
      nrrank++;
    }
    max++;
  }

  if (VERB3)
    fprintf(stderr,"     ready ranking, %d ranked\n",nrrank);
  thisrankid[nrrank]=-1;
}

/* compute_instance_weights -
 */
void compute_instance_weights(void)
{
  register int i;

  /* now determine the order of things: CPS with Laplace */
  for (i=0; i<NRPAT; i++)
  {
    ordernr[i]=
      1.*(unsigned int) (lrand48()*((unsigned int)time(NULL))%UINT_MAX);
    order[i]=i;
  }
}
