/** Classify.c

Copyright 1997-2004 ILK / Tilburg University
Written by Antal van den Bosch */

#include<stdio.h>
#include<string.h>
#include<stddef.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<limits.h>

#include "Common.h"

int   *candidateclass;
float *candidateoccur;
float *candidatematch;
int   *candidateid;
int   *nrnn;
int   *nrcornn;
int   *candid;
float  *eweights;

#include "Classify.h"

void alloccandidates(void)
{
  int nrcand;

  if (nrexpressions>MAXNREXP)
    nrcand=MAXNREXP;
  else
    nrcand=nrexpressions+1;
  candidateclass=(int*)malloc(nrcand*sizeof(int));
  if (candidateclass==NULL)
    alloc_error();
  candidateoccur=(float*)malloc(nrcand*sizeof(float));
  if (candidateoccur==NULL)
    alloc_error();
  candidatematch=(float*)malloc(nrcand*sizeof(float));
  if (candidatematch==NULL)
    alloc_error();
  candidateid=(int*)malloc(nrcand*sizeof(int));
  if (candidateid==NULL)
    alloc_error();
}

void freecandidates(void)
{
  free(candidateclass);
  free(candidateoccur);
  free(candidatematch);
  free(candidateid);
}

value klassify(void)
{
  int   i,j,k;
  int   nrcandidates,rank,
        nrclasslist,nrsuperclasslist,
        max,filledk=0;
  char  in,same;
  value classlist[MAXNRCLASSES];
  value superclasslist[MAXNRCLASSES];
  value returnclass=0,thisclass;
  float classoccs[MAXNRCLASSES];
  float nearestn,furthestn,maxocc;

  if (VERB2)
  {
    fprintf(stderr,"     classifying ");
    for (i=0; i<PATWIDTH; i++)
    {
      if (thispattern[i]==-1)
	fprintf(stderr,"UNKNOWN_VALUE,");
      else
	fprintf(stderr,"%s,",values[i][thispattern[i]]);
    }
    fprintf(stderr,"\n");
  }

  nrcandidates=0;
  mindist=0.0;
  identical=0.0;
  candidatematch[0]=0.0;
  filledk=0;

  max=0;

  while (max<nrexpressions)
  {
    sweep_c(max);

    if (EWEIGHT)
      identical*=eweights[max];

    if //((identical>0.0)&&
      ((identical>mindist-ALMOSTNOTHING)||(filledk<FAMBLKA))
      {
	thisclass=expressionclass[max];

	rank=0;
	while ((identical<candidatematch[rank])&&
	       (rank<nrcandidates)&&
	       (rank<MAXNREXP))
	  {
	    rank++;
	  }

	if (rank<nrcandidates)
	  {
	    for (k=nrcandidates; k>rank; k--)
	      {
		candidateclass[k]=candidateclass[k-1];
		candidatematch[k]=candidatematch[k-1];
		candidateoccur[k]=candidateoccur[k-1];
		candidateid[k]=candidateid[k-1];
	      }
	  }

	candidateclass[rank]=thisclass;
	candidatematch[rank]=identical;
	candidateoccur[rank]=thisoccurrence;
	candidateid[rank]=max;
	nrcandidates++;

	if (nrcandidates>1)
	  {
	    filledk=1;
	    for (i=1; i<nrcandidates; i++)
	      if (candidatematch[i]!=candidatematch[i-1])
		filledk++;

	    if (filledk>FAMBLKA)
	      {
		k=nrcandidates-1;
		while ((k>1)&&
		       (candidatematch[k]==candidatematch[k-1]))
		  k--;
		k--;
		//mindist=candidatematch[k];
		nrcandidates=k+1;
		filledk--;
	      }
	  }

	mindist=candidatematch[nrcandidates-1];

	if (VERB4)
	  {
	    fprintf(stderr,"     %d distances (lowest %f), ranked as follows:\n",
		   filledk,mindist);
	    for (i=0; i<nrcandidates; i++)
	      {
		fprintf(stderr,"     %3d. identical %f\n",
		       i,candidatematch[i]);
	      }
	  }

      }
    max++;
  }

  // distance weighting
  if (DISTANCE>0)
    {
      for (k=0; k<nrcandidates; k++)
      {
	if (DISTANCE==1)
	  {
	    nearestn=maxgrsum-candidatematch[0];
	    furthestn=maxgrsum-candidatematch[nrcandidates-1];
	    if (nearestn!=furthestn)
	      candidateoccur[k]=
		candidateoccur[k]*
		((furthestn-(maxgrsum-candidatematch[k]))/
		 (furthestn-nearestn));
	  }
	if (DISTANCE==2)
	  {
	    candidateoccur[k]=
	      candidateoccur[k]*
	      (1/((maxgrsum-candidatematch[k])+ALMOSTNOTHING));
	  }
	if (DISTANCE==3)
	  {
	    candidateoccur[k]=
	      candidateoccur[k]*
	      (exp((double) -1.*(maxgrsum-candidatematch[k])));
	  }
      }
    }

  if (VERB3)
    {
      fprintf(stderr,"     %d candidates found\n",
             nrcandidates);
      if (nrcandidates>0)
	for (k=0; k<nrcandidates; k++)
	  {
	    fprintf(stderr,"     %3d. class %s, occ. %.4f, ident %.4f = dist. %.4f\n",
		   k+1,
		   classes[candidateclass[k]],
		   candidateoccur[k],
		   candidatematch[k],
		   maxgrsum-candidatematch[k]);
	  }
    }

  if (nrcandidates==0)
    {
      returnclass=mostfreqclass;
      returndist=maxgrsum;
    }
  else
    {

      if (distprint)
	returndist=maxgrsum-candidatematch[0];


      /* we are going for regular k-nn */
      same=1;
      for (i=1; ((same)&&(i<nrcandidates)); i++)
	if (candidateclass[i-1]!=candidateclass[i])
	  same=0;

      if (same)
	{
	  returnclass=candidateclass[0];
	  if (VERB4)
	    fprintf(stderr,"     unambiguous expression set: output class is %s\n",
		   classes[returnclass]);
	}
      else
	{
	  if (VERB4)
	    fprintf(stderr,"     class ambiguity in expression set: %d classes\n",
		   nrcandidates);

	  nrclasslist=0;
	  for (i=0; i<nrcandidates; i++)
	    {
	      in=0;
	      for (j=0; ((!in)&&(j<nrclasslist)); j++)
		if (classlist[j]==candidateclass[i])
		  {
		    in=1;
		    classoccs[j]+=candidateoccur[i];
		  }
	      if (!in)
		{
		  classlist[nrclasslist]=candidateclass[i];
		  classoccs[nrclasslist]=candidateoccur[i];
		  nrclasslist++;
		}
	    }
	  maxocc=0.0;
	  nrsuperclasslist=0;
	  for (i=0; i<nrclasslist; i++)
	    {
	      if (classoccs[i]==maxocc)
		{
		  superclasslist[nrsuperclasslist]=classlist[i];
		  nrsuperclasslist++;
		}
	      if (classoccs[i]>maxocc)
		{
		  superclasslist[0]=classlist[i];
		  nrsuperclasslist=1;
		  maxocc=classoccs[i];
		}
	    }

	  if (nrsuperclasslist==1)
	    {
	      /* if (VERB3)
		 fprintf(stderr,"     one class, %s, occurs most often\n",
		 classes[superclasslist[0]]); */
	      returnclass=superclasslist[0];
	    }
	  else
	    {
	      /* if (VERB3)
		 fprintf(stderr,"     WE HAVE A TIE, between\n"); */
	      maxocc=0.0;
	      for (j=0; j<nrsuperclasslist; j++)
		{
		  k=0;
		  while (superclasslist[j]!=k)
		    k++;
		  if ((1.*overallclassocc[k])>maxocc)
		    {
		      maxocc=(1.*overallclassocc[k]);
		      returnclass=superclasslist[j];
		    }
		  /* if (VERB3)
		     fprintf(stderr,"         %s, overall %d\n",
		     classes[k],overallclassocc[k]); */
		}
	      /* if (VERB3)
		 fprintf(stderr,"     utter tie won by %s\n",
		 classes[returnclass]); */
	    }
	}

      if (VERB3)
	fprintf(stderr,"     => class %s\n",
	       classes[returnclass]);
    }
  return returnclass;
}

void klassify_eweight(value thiscat)
{
  int   i,k;
  int   nrcandidates,rank,limit,max,
        EWEIGHTKA,afterexact;
  value thisclass;

  if (VERB2)
    {
      fprintf(stderr,"     classifying ");
      for (i=0; i<PATWIDTH; i++)
	{
	  if (thispattern[i]==-1)
	    fprintf(stderr,"UNKNOWN_VALUE,");
	  else
	    fprintf(stderr,"%s,",values[i][thispattern[i]]);
	}
      fprintf(stderr,"\n");
    }

  EWEIGHTKA=FAMBLKA+1;

  nrcandidates=0;
  mindist=0.0;
  identical=0.0;

  max=0;

  while (max<nrexpressions)
  {
    sweep_c(max);

    if (identical>mindist-ALMOSTNOTHING)
    {
      thisclass=expressionclass[max];

      rank=0;
      while ((identical<candidatematch[rank])&&
	     (rank<nrcandidates)&&
	     (rank<MAXNREXP))
	rank++;

      if (rank<nrcandidates)
	{
	  for (k=nrcandidates; k>rank; k--)
	    {
	      candidateclass[k]=candidateclass[k-1];
	      candidatematch[k]=candidatematch[k-1];
	      candidateoccur[k]=candidateoccur[k-1];
	      candid[k]=candid[k-1];
	    }
	}
      candidateclass[rank]=thisclass;
      candidatematch[rank]=identical;
      candidateoccur[rank]=thisoccurrence;
      candid[rank]=max;
      nrcandidates++;

      limit=1;
      k=0;
      while ((k<EWEIGHTKA)&&(limit<nrcandidates))
	{
	  if (candidatematch[limit-1]>candidatematch[limit])
	    k++;
	  if (k<EWEIGHTKA)
	    limit++;
	}

      nrcandidates=limit;
      mindist=candidatematch[nrcandidates-1];
    }
    max++;
  }

  /* skip the exact matches */
  afterexact=0;
  while (candidatematch[afterexact]>maxgrsum-ALMOSTNOTHING)
     afterexact++;

  if (VERB3)
    {
      fprintf(stderr,"     %d candidates found\n",
             nrcandidates-afterexact);
      if (nrcandidates-afterexact>0)
	for (k=afterexact; k<nrcandidates; k++)
	  {
	    fprintf(stderr,"     %3d. class %s, occ. %.2f, match %f\n",
		   k-afterexact,
		   classes[candidateclass[k]],
		   candidateoccur[k],
		   candidatematch[k]);
	  }
    }

  /* update the CPS counters */
  for (i=afterexact; i<nrcandidates; i++)
    {
      nrnn[candid[i]]++;
      if (candidateclass[i]==thiscat)
	nrcornn[candid[i]]++;
    }


  if (EWEIGHT==1)
    eweights[i]=(1.*nrcornn[i])/(1.*nrnn[i]);
  if (EWEIGHT==2)
    eweights[i]=((1.*nrcornn[i])+1.)/((1.*nrnn[i])+(1.*NRCLASSES));

}

void perform_test( const char testname[NAMELEN] )
{
  FILE  *testfile,*outfile,*percfile=NULL;
  char  *readvalue;
  char  outname[1024];
  char  percname[1024];
  char  kpart[10];
  int   i,j,corpat,returnclass,readcat=0;
  char  *line;
  char  *linecopy;
  char  delimiter;
  char  deltok[10];
  char  thisclass[MAXVALLEN+1];
  float numericreadvalue;

  line=(char*)malloc(500000*sizeof(char));
  linecopy=(char*)malloc(500000*sizeof(char));

  strcpy(outname,testname);
  if (FAMKA==0)
    strcat(outname,".IB1");
  else
    strcat(outname,".Fambl");

  if (METRIC==0)
    strcat(outname,".O");
  if (METRIC==1)
    strcat(outname,".M");
  if (METRIC==2)
    strcat(outname,".J");

  if (MVDML!=1)
    {
      sprintf(kpart,".L%d",
	      MVDML);
      strcat(outname,kpart);
    }

  if (FWEIGHT==0)
    strcat(outname,".nw");
  if (FWEIGHT==1)
    strcat(outname,".ig");
  if (FWEIGHT==2)
    strcat(outname,".gr");
  if (FWEIGHT==3)
    strcat(outname,".x2");
  if (FWEIGHT==4)
    strcat(outname,".sv");
  if (FWEIGHT==5)
    strcat(outname,".ll");
  if (FWEIGHT==6)
    strcat(outname,".tree");
  if (FWEIGHT==7)
    strcat(outname,".tribl");
  if ((FWEIGHT==8)||(FWEIGHT==9))
    strcat(outname,".ud");

  sprintf(kpart,".k%d",
	  FAMBLKA);
  strcat(outname,kpart);

  sprintf(kpart,".K%d",
	  FAMKA);
  strcat(outname,kpart);


  if (DISTANCE)
    {
      if (DISTANCE==1)
	strcat(outname,".IL");
      if (DISTANCE==2)
	strcat(outname,".ID");
      if (DISTANCE==3)
	strcat(outname,".ED");
    }

  if (EWEIGHT==1)
    strcat(outname,".cps");
  if (EWEIGHT==2)
    strcat(outname,".cpsl");
  if (WILDCARD)
    strcat(outname,".wild");
  if (WILDCARDTHRESHOLD!=1000)
    {
      sprintf(kpart,".x%d",
	      WILDCARDTHRESHOLD);
      strcat(outname,kpart);
    }

  if (inclusive)
    strcat(outname,".X");

  if (ALGORITHM)
    {
      strcat(outname,".bin");
      if (NRCOMBINATIONS>1)
	{
	  sprintf(kpart,".C%d",
		  NRCOMBINATIONS);
	  strcat(outname,kpart);
	}
    }

  if (WRITEPERC)
    {
      strcpy(percname,outname);
      strcat(percname,".%");
    }

  strcat(outname,".out");

  if (VERB)
    {
      fprintf(stderr," <*> starting test with k=%d\n     writing output to %s\n",
	     FAMBLKA,outname);
      if (WRITEPERC)
	fprintf(stderr,"     also writing classification accuracy to %s\n",
	       percname);
    }

  totpat=corpat=0;

  // determine delimiter!
  testfile=fopen(testname,"r");
  if (testfile==NULL)
    {
      fprintf(stderr," <!> test file %s not found\n",
	      testname);
      exit(1);
    }
  if ( fgets(line,MAXEXPRESSION,testfile) == NULL ){
    fprintf( stderr, "read failure on file %s\n", testname );
    exit(1);
  }

  if (strstr(line," "))
    delimiter=' ';
  else
    {
      if (strstr(line,"\t"))
        delimiter='\t';
      else
        delimiter=',';
    }
  sprintf(deltok,"%c\n",
          delimiter);
  fclose(testfile);
  if (VERB3)
    fprintf(stderr,"     delimiter: [%c]\n",
            delimiter);

  outfile=fopen(outname,"w");
  if (WRITEPERC)
    percfile=fopen(percname,"w");

  setbuf(outfile,NULL);

  alloccandidates();

  /* loop through the testfile */
  testfile=fopen(testname,"r");
  if ( fgets(line,MAXEXPRESSION,testfile) == NULL ){
    fprintf( stderr, "read failure on file %s\n", testname );
    exit(1);
  }
  while (!feof(testfile))
  {
    if ((int)strlen(line)>((2*PATWIDTH)+1))
      {
	strcpy(linecopy,line);

	readvalue=strtok(line,deltok);
	for (i=0; i<PATWIDTH; i++)
	  {
	    /* symbolic? */
	    if (feattype[i]=='s')
	      {
		j=0;
		while ((j<nrvalues[i])&&(strcmp(values[i][j],readvalue)!=0))
		  j++;
		if (j<nrvalues[i])
		  thispattern[i]=j;
		else
		  thispattern[i]=-1;
	      }
	    /* numeric? */
	    else
	      {
		sscanf(readvalue," %f",&numericreadvalue);
		thispattern[i]=(value)
		  (((numericreadvalue-ffmin[i])*(1.*DISCRETIZER))/
		   (ffmax[i]-ffmin[i]));
	      }
	    readvalue=strtok(NULL,deltok);

	  }
	strcpy(thisclass,readvalue);
	j=0;
	while ((j<NRCLASSES)&&(strcmp(classes[j],thisclass)!=0))
	  j++;
	readcat=j;

	/* if working with bonuses, determine the subset of applicable bonus
	   selections */
	if (ALGORITHM)
	  select_applicable_bonuses();

	/* do the actual classification */
	if ((NRCLASSES==0)||
	    ((ALGORITHM)&&(nrbonuses==0)))
	  returnclass=0;
	else
	  returnclass=klassify();

	/* (note: it happens that a data set contains just one class or
	   no features. */

	if (VERB2)
	  fprintf(stderr,"     target: %s, predicted: %s\n",
		  classes[readcat],classes[returnclass]);

	totpat++;

	/* score */
	if (readcat==returnclass)
	  corpat++;
	else
	  if (VERB3)
	    {
	      fprintf(stderr,"     mismatch! ");
	      for (i=0; i<(int)strlen(linecopy)-1; i++)
		fprintf(stderr,"%c",linecopy[i]);
	      fprintf(stderr," %s\n",classes[returnclass]);
	    }

	/* send classification to output file (timbl style) */
	j=(int)strlen(linecopy)-2;
	if (linecopy[j]=='.')
	  j--;
	for (i=0; i<=j; i++)
	  fprintf(outfile,"%c",linecopy[i]);
	fprintf(outfile,",%s",classes[returnclass]);
	if (distprint)
	  fprintf(outfile," %f",
		  returndist);
	fprintf(outfile,"\n");
      }

    if ((VERB)&&(!(totpat%REPORTRATE))&&(totpat>0))
      fprintf(stderr,"     %7d instances processed, %8.4f %% correct\n",
	      totpat,(100.*corpat)/(1.*totpat));
    if ( fgets(line,MAXEXPRESSION,testfile) == NULL && !feof(testfile)){
      fprintf( stderr, "read failure on file %s\n", testname );
      exit(1);
    }
  }
  fclose(testfile);
  fclose(outfile);

  /* report on number and percentage of correct classifications */
  if (VERB)
    {
      fprintf(stderr,"     %7d instances out of %7d classified correctly\n",
	     corpat,totpat);
      fprintf(stderr,"     Fambl score: %8.4f %% correct instances\n",
	     (100.*corpat)/(1.*totpat));
    }
  else
    fprintf(stderr,"%.4f\n",
	   (100.*corpat)/(1.*totpat));

  if (WRITEPERC)
    {
      fprintf(percfile,"%f\n",
	      (100.*corpat)/(1.*totpat));
      fclose(percfile);
    }

  /* candidate buffer is ready to be liberated */
  freecandidates();
  free(line);
  free(linecopy);
}

void compute_exemplar_weights( const char dataname[NAMELEN],
			       const char cpsname[NAMELEN] )
{
  FILE  *testfile;
  FILE  *cpsfile=NULL;
  value readcat=0;
  char  *readvalue;
  int   i,j,k;
  char  delimiter;
  char  deltok[10];
  char  *line;
  char  *linecopy;
  char  thisclass[MAXVALLEN+1];

  line=(char*)malloc(500000*sizeof(char));
  linecopy=(char*)malloc(500000*sizeof(char));

  if (VERB)
    fprintf(stderr," <*> starting computation of family weights\n");

  nrnn=(int*)malloc(nrexpressions*sizeof(int));
  if (nrnn==NULL)
    alloc_error();
  nrcornn=(int*)malloc(nrexpressions*sizeof(int));
  if (nrcornn==NULL)
    alloc_error();
  eweights=(float*)malloc(nrexpressions*sizeof(float));
  if (eweights==NULL)
    alloc_error();
  candid=(int*)malloc(nrexpressions*sizeof(int));
  if (candid==NULL)
    alloc_error();
  alloccandidates();

  totpat=0;
  for (i=0; i<nrexpressions; i++)
    {
      nrnn[i]=0;
      nrcornn[i]=0;
    }

  // determine delimiter
  testfile=fopen(dataname,"r");
  if ( fgets(line,MAXEXPRESSION,testfile) == NULL ){
    fprintf( stderr, "read failure on file %s\n", dataname );
    exit(1);
  }
  if (strstr(line," "))
    delimiter=' ';
  else
    {
      if (strstr(line,"\t"))
        delimiter='\t';
      else
        delimiter=',';
    }
  sprintf(deltok,"%c\n",
          delimiter);
  fclose(testfile);
  if (VERB3)
    fprintf(stderr,"     delimiter: [%c]\n",
            delimiter);


  testfile=fopen(dataname,"r");
  if ( fgets(line,MAXEXPRESSION,testfile) == NULL ){
    fprintf( stderr, "read failure on file %s\n", dataname );
    exit(1);
  }
  while (!feof(testfile))
    {
      if ((int)strlen(line)>((2*PATWIDTH)+1))
	{
	  strcpy(linecopy,line);

	  /*
	  for (i=0; i<((int)strlen(line)-2); i++)
	    if (line[i]=='.')
	      line[i]='P';
	  */

	  readvalue=strtok(line,deltok);
	  for (i=0; i<PATWIDTH; i++)
	    {
	      j=0;
	      while ((j<nrvalues[i])&&(strcmp(values[i][j],readvalue)!=0))
		j++;
	      thispattern[i]=j;
	      readvalue=strtok(NULL,deltok);
	    }
	  strcpy(thisclass,readvalue);
	  j=0;
	  while ((j<NRCLASSES)&&(strcmp(classes[j],thisclass)!=0))
	    j++;
	  readcat=j;

	  klassify_eweight(readcat);
	  totpat++;
	}

      if ((VERB)&&(!(totpat%REPORTRATE))&&(totpat>0))
	fprintf(stderr,"     %7d training instances processed\n",
	       totpat);
      if ( fgets(line,MAXEXPRESSION,testfile) == NULL && !feof(testfile)){
	fprintf( stderr, "read failure on file %s\n", dataname );
	exit(1);
      }
    }
  fclose(testfile);

  /* the relevant bit: computation of CPS values with or without Laplacian
     correction */

  if (PRINTCPS)
    cpsfile=fopen(cpsname,"w");

  for (i=0; i<nrexpressions; i++)
    {
      if (nrnn[i]>0)
	{
	  if (EWEIGHT==1)
	    eweights[i]=(1.*nrcornn[i])/(1.*nrnn[i]);
	  if (EWEIGHT==2)
	    eweights[i]=((1.*nrcornn[i])+1.)/((1.*nrnn[i])+(1.*NRCLASSES));
	}
      else
	eweights[i]=0.0;
      if (VERB3)
	fprintf(stderr,"     family %6d: %3d NN, %3d correct NN, CPS %6.4f\n",
	       i,nrnn[i],nrcornn[i],eweights[i]);

      if (PRINTCPS)
	{
	  for (j=0; j<PATWIDTH; j++)
	    {
	      if (feattype[j]=='s')
		{
		  for (k=expression[i][j];
		       k<expression[i][j+1];
		       k++)
		    {
		      if (expression[i][k]==WILDCARDVALUE)
			fprintf(cpsfile,"*,");
		      else
			fprintf(cpsfile,"%s,",values[j][expression[i][k]]);
		    }
		}
	      else
		{
		  if (expression[i][j]==expression[i][j+1]+1)
		    fprintf(cpsfile,"%6.4f,",
			   numscale(expression[i][expression[i][j]]));
		  else
		    fprintf(cpsfile,"%6.4f-%6.4f,",
			   numscale(expression[i][expression[i][j]]),
			   numscale(expression[i][expression[i][j]+1]));
		}
	    }
	  fprintf(cpsfile,"%s\t%f\n",
		 classes[expressionclass[i]],
		 eweights[i]);
	}
    }

  if (PRINTCPS)
    fclose(cpsfile);

  /* get rid of all unnecessary allocated buffers */
  free(line);
  free(linecopy);
  free(nrnn);
  free(nrcornn);
  free(candid);
  freecandidates();
}
