/* Read.c

   Copyright 1997-2004 ILK / Tilburg University
   Written by Antal van den Bosch
*/

#include<stdio.h>
#include<string.h>
#include<stddef.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<limits.h>

#include "Common.h"
#include "Read.h"

void readin_fambl( const char famblname[NAMELEN] )
{
  FILE  *famblfile;
  int   i,j,count=0;
  size_t res;
  value readbuf[MAXPATWIDTH+1];

  if (VERB)
    fprintf(stderr,"     reading Fambl file %s\n",
	   famblname);

  nrexpressions=maxexpfound=0;

  for (i=0; i<NRCLASSES; i++)
    overallclassocc[i]=0;

  expressionclass=(value*)malloc(MAXNREXP*sizeof(value));
  expressionoccurrence=(int*)malloc(MAXNREXP*sizeof(int));

  famblfile=fopen(famblname,"r");
  while (!feof(famblfile))
    {
      for (i=0; i<=PATWIDTH; i++){
	res = fread(&readbuf[i],
		    sizeof(value),
		    1,
		    famblfile);
	if ( res == 0 ){
	  fprintf( stderr, "read error on file %s\n", famblname );
	}
      }
      expression[count]=(value*)malloc(readbuf[PATWIDTH]*sizeof(value));
      maxexpfound=max(readbuf[PATWIDTH],maxexpfound);
      for (i=0; i<=PATWIDTH; i++)
	expression[count][i]=readbuf[i];
      for (i=PATWIDTH+1; i<readbuf[PATWIDTH]; i++)
	{
	  res = fread(&expression[count][i],
		      sizeof(value),
		      1,
		      famblfile);
	  if ( res == 0 ){
	    fprintf( stderr, "read error on file %s\n", famblname );
	  }
	}
      res = fread(&expressionclass[count],
		  sizeof(value),
		  1,
		  famblfile);
      if ( res == 0 ){
	fprintf( stderr, "read error on file %s\n", famblname );
      }
      res = fread(&expressionoccurrence[count],
		  sizeof(int),
		  1,
		  famblfile);
      if ( res == 0 ){
	fprintf( stderr, "read error on file %s\n", famblname );
      }
      overallclassocc[expressionclass[count]]+=expressionoccurrence[count];

      if (VERB4)
	{
	  fprintf(stderr,"     expression %d:\n     ",count);
	  for (i=0; i<PATWIDTH; i++)
	    {
	      for (j=expression[count][i];
		   j<expression[count][i+1];
		   j++)
		fprintf(stderr,"%s,",
			values[i][expression[count][j]]);
	      fprintf(stderr," ");
	    }
	  fprintf(stderr,"\n     class %d/%s, occurrence %d\n\n",
		  expressionclass[count],classes[expressionclass[count]],
		  expressionoccurrence[count]);
	}

      count++;
    }
  nrexpressions=count-1;
  if (VERB)
    fprintf(stderr,"     read %d Fambl expressions\n",
	   nrexpressions);
  fclose(famblfile);
  expressionclass=(value*)realloc(expressionclass,nrexpressions*sizeof(value));
  expressionoccurrence=(int*)realloc(expressionoccurrence,nrexpressions*sizeof(int));
}

void readin_binary_weights( const char filename[NAMELEN])
{
  FILE  *data;
  int   i,j,featnumber;
  value valuenumber;
  char  readvalue[1024];
  char  copyvalue[1024];
  char  *part;
  char  match;
  char  featnumberstring[32];
  char  valuestring[1024];
  float readweight;
  size_t res;
  // set all atomic weights to 0.0
  for (i=0; i<MAXBONUSES; i++)
    {
      bonuses[i]=0.0;
    }
  nrbonuses=0;

  if (VERB)
    fprintf(stderr,"     reading binary feature weights from file %s\n",
	    filename);

  data=fopen(filename,"r");
  if (data==NULL)
    {
      fprintf(stderr," <!> %s: file not found\n\n",
	      filename);
      exit(1);
    }
  while (!feof(data))
    {
      res = fscanf(data,"%s %f ",
		   readvalue,&readweight);
      if ( res == 0 )
	{
	  fprintf( stderr, "read error from %s\n", filename );
	}
      strcpy(copyvalue,readvalue);
      part=strtok(copyvalue,"|\0");
      strcpy(featnumberstring,part);
      sscanf(featnumberstring,"%d",&featnumber);
      part=strtok(NULL,"|\0");
      strcpy(valuestring,part);
      if (VERB4)
	fprintf(stderr,"     read %s (weight %f),\n       containing feature number %d, value %s\n",
		readvalue,readweight,featnumber,valuestring);
      valuenumber=-1;
      j=0;
      match=0;
      while ((j<nrvalues[featnumber])&&(!match))
	{
	  if (values[featnumber][j][0]==readvalue[0])
	    if (strcmp(values[featnumber][j],readvalue)==0)
	      {
		match=1;
		valuenumber=j;
	      }
	  j++;
	}
      if (VERB4)
	fprintf(stderr,"     found value %d of feature %d\n",
		valuenumber,featnumber);

      bonusselections[nrbonuses]=
	(value*)malloc(PATWIDTH*sizeof(value));
      if (bonusselections[nrbonuses]==NULL)
	alloc_error();

      for (i=0; i<PATWIDTH; i++)
	{
	  if (i==featnumber)
	    bonusselections[nrbonuses][i]=valuenumber;
	  else
	    bonusselections[nrbonuses][i]=-1;
	}

      bonuses[nrbonuses]=readweight;
      if (VERB4)
	{
	  fprintf(stderr,"     generated the following bonusselection with weight %f:\n",
		  bonuses[nrbonuses]);
	  for (i=0; i<PATWIDTH; i++)
	    {
	      fprintf(stderr,"     %2d. %d\n",
		      i,bonusselections[nrbonuses][i]);
	    }
	}
      nrbonuses++;
    }

  /* allocate the applicable bonuses buffer */
  applicable_bonuses=(int*)malloc(nrbonuses*sizeof(int));
  if (applicable_bonuses==NULL)
    {
      alloc_error();
      exit(1);
    }

  if (VERB)
  {
    fprintf(stderr,"\r     %d binary features found\n",
	    nrbonuses);
  }
}

void readin_weights( const char filename[NAMELEN] )
{
  FILE  *data;
  int   i,readint,linenr=1;
  float readweight;
  size_t res;
  // set all weights to 0.0
  for (i=0; i<PATWIDTH; i++)
    featgrvals[i]=0.0;

  if (VERB)
    fprintf(stderr,"     reading feature weights from file %s\n",
	    filename);

  data=fopen(filename,"r");
  if (data==NULL)
    {
      fprintf(stderr," <!> %s: file not found\n\n",
	      filename);
      exit(1);
    }
  while (!feof(data))
    {
      res = fscanf(data,"%d %f ",
		   &readint,&readweight);
      if ( (res == 0) ||
	   (readint<1)||
	   (readint>PATWIDTH))
	{
	  fprintf(stderr," <!> erroneous feature number found at line %d of %s\n\n",
		  linenr,filename);
	  exit(1);
	}
      if (readweight<0.0)
	{
	  fprintf(stderr," <!> no negative weights allowed; found one at line %d of %s\n\n",
		  linenr,filename);
	  exit(1);
	}

      // set the weight
      featgrvals[readint-1]=readweight;
      if (VERB3)
	fprintf(stderr,"     feature %3d: weight %f\n",
		readint-1,readweight);

      linenr++;
    }
  fclose(data);

}

void readin_data( const char dataname[NAMELEN] )
{
  FILE  *data;
  int   i,j;
  value nrreadclasses=0,code=0;
  int   counter=0,hit=0,ffcounter,topocc;
  char  *line;
  char  *previousline;
  char  *dummy;
  char  in;
  char  *thisvalue;
  char  thisclass[MAXVALLEN+1];
  char  commandline[1024];
  char  tmpdataname[50];
  char  delimiter;
  char  deltok[10];
  int   ffindex[MAXPATWIDTH];
  float *floatfeature[MAXPATWIDTH];

  if (VERB)
    fprintf(stderr,"     sorting and reading data base %s\n",
           dataname);

  // initialize read buffers
  line=(char*)malloc(MAXEXPRESSION*sizeof(char));
  previousline=(char*)malloc(MAXEXPRESSION*sizeof(char));
  strcpy(previousline,"#");

  // determine delimiter
  data=fopen(dataname,"r");
  if ( fgets(line,MAXEXPRESSION,data) == NULL ){
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
  fclose(data);
  if (VERB3)
    fprintf(stderr,"     delimiter: [%c]\n",
	    delimiter);

  /* first, create a shadow data file that is sorted.
     it's given a name with six random digits: e.g. fambl123456.tmp. */

  sprintf(tmpdataname,"fambl%d.tmp",
	  (int) (lrand48()*((unsigned int)time(NULL))%UINT_MAX));

  /* do a simple sort to make looking for clones easier */
  sprintf(commandline,"sort %s -o %s -T .\n",
          dataname,tmpdataname);
  if (VERB4)
    fprintf(stderr,"     executing commandline:\n     %s",commandline);
  my_system(commandline);

  /* readin the data file, while counting types */

  NRPAT=nrtypes=0;
  data=fopen(tmpdataname,"r");
  if ( fgets(line,MAXEXPRESSION,data) == NULL ){
    fprintf( stderr, "read failure on file %s\n", tmpdataname );
    exit(1);
  }
  while (!feof(data))
    {
      NRPAT++;
      if (strcmp(previousline,line)!=0)
	{
	  nrtypes++;
	  strcpy(previousline,line);
	}
      if ( fgets(line,MAXEXPRESSION,data) == NULL && !feof(data) ){
	fprintf( stderr, "read failure on file %s\n", tmpdataname );
	exit(1);
      }
    }
  fclose(data);

  PATWIDTH=0;
  dummy=strtok(line,deltok);
  while (dummy!=NULL)
    {
      dummy=strtok(NULL,deltok);
      if (dummy!=NULL)
	PATWIDTH++;
    }

  /* when PATWIDTH is not in accordance with a user-specified
     feature type string, object and quit. if the user did not
     specify anything, set the feature type string to "all symbolic"
     (default). */
  if (strlen(feattype)>0)
    {
      if ((int)strlen(feattype)!=PATWIDTH)
	{
	  fprintf(stderr," <!> prespecified number of features (%d) conflicts with data (%d).\n\n",(int)strlen(feattype),PATWIDTH);
	  exit(0);
	}
      for (i=0; i<PATWIDTH; i++)
	{
	  if (feattype[i]=='n')
	    {
	      floatfeature[nrnumerics]=(float*)malloc(NRPAT*sizeof(float));
	      if (floatfeature[nrnumerics]==NULL)
		alloc_error();
	      ffindex[nrnumerics]=i;
	      ffmax[nrnumerics]=(float) INT_MIN;
	      ffmin[nrnumerics]=(float) INT_MAX;
	      nrnumerics++;
	    }
	}
    }
  else
    {
      for (i=0; i<PATWIDTH; i++)
	{
	  strcat(feattype," ");
	  if (numerics)
	    {
	      feattype[i]='n';
	      floatfeature[nrnumerics]=(float*)malloc(NRPAT*sizeof(float));
	      if (floatfeature[nrnumerics]==NULL)
		alloc_error();
	      ffindex[nrnumerics]=i;
	      ffmax[nrnumerics]=(float) INT_MIN;
	      ffmin[nrnumerics]=(float) INT_MAX;
	      nrnumerics++;
	    }
	  else
	    feattype[i]='s';
	}
    }

  /* allocate feature value and class dynamic buffers */

  for (i=0; i<PATWIDTH; i++)
    {
      feature[i]=(value*)malloc((NRPAT+1)*sizeof(value));
      if (feature[i]==NULL)
	alloc_error();
    }
  klass=(value*)malloc((NRPAT+1)*sizeof(value));
  if (klass==NULL)
    alloc_error();

  if (VERB2)
    fprintf(stderr,"     allocation ready; accessing & reading file\n");

  // start reading data
  data=fopen(tmpdataname,"r");
  if (VERB2)
    fprintf(stderr,"     file %s opened\n",tmpdataname);
  for (i=0; i<PATWIDTH; i++)
    nrvalues[i]=0;
  nrreadclasses=0;

  if ( fgets(line,MAXEXPRESSION,data) == NULL ){
    fprintf( stderr, "read failure on file %s\n", tmpdataname );
    exit(1);
  }
  while (!feof(data))
  {
    if (VERB2)
      if ((!(counter%(REPORTRATE*100)))&&(counter>0))
      { fprintf(stderr,"     %7d instances read\n",counter);
        fflush(stdout);
      }

    thisvalue=strtok(line,deltok);
    ffcounter=0;
    for (i=0; i<PATWIDTH; i++)
      {
	/* symbolic? get the code */
	if (feattype[i]=='s')
	  {
	    in=0;
	    for (j=0; ((j<nrvalues[i])&&(!in)); j++)
	      if (thisvalue[0]==values[i][j][0])
		if (strcmp(thisvalue,values[i][j])==0)
		  {
		    in=1;
		    code=j;
		  }
	    if (!in)
	      {
		values[i][nrvalues[i]]=
		  (char*)malloc((strlen(thisvalue)+1)*sizeof(char));
		strcpy(values[i][nrvalues[i]],thisvalue);
		code=nrvalues[i];
		nrvalues[i]++;
	      }
	    feature[i][counter]=code;
	  }

	/* numeric? get value and determine min & max */
	else
	  {
	    sscanf(thisvalue," %f",&floatfeature[ffcounter][counter]);
	    ffmax[ffcounter]=max(floatfeature[ffcounter][counter],ffmax[ffcounter]);
	    ffmin[ffcounter]=min(floatfeature[ffcounter][counter],ffmin[ffcounter]);
	    ffcounter++;
	  }
	thisvalue=strtok(NULL,deltok);
      }
    strcpy(thisclass,thisvalue);

    in=0;
    if (nrreadclasses>0)
      for (j=0; ((j<nrreadclasses)&&(!in)); j++)
	if (thisclass[0]==classes[j][0])
	  if (strcmp(thisclass,classes[j])==0)
	    {
	      in=1;
	      code=j;
	      hit=j;
	    }
    if (in)
      overallclassocc[hit]++;
    else
      {
	classes[nrreadclasses]=
	  (char*)malloc((strlen(thisclass)+1)*sizeof(char));
	strcpy(classes[nrreadclasses],thisclass);
	overallclassocc[nrreadclasses]=1;
	code=nrreadclasses;
	nrreadclasses++;
      }

    klass[counter]=code;
    counter++;

    if ( fgets(line,MAXEXPRESSION,data) == NULL && !feof(data) ){
      fprintf( stderr, "read failure on file %s\n", tmpdataname );
      exit(1);
    }
  }

  NRCLASSES=nrreadclasses;
  fclose(data);

  if ((VERB)&&(notshowntwice))
    {
      fprintf(stderr,"     data base has %d instances\n",
	      NRPAT);
      fprintf(stderr,"                   %d instance types\n",
	      nrtypes);
      fprintf(stderr,"                   %d features\n",
	      PATWIDTH);
      fprintf(stderr,"                   %d classes\n",
	      NRCLASSES);

      notshowntwice=0;
    }

  /* remove shadow data file */

  sprintf(commandline,"rm %s\n",tmpdataname);
  my_system(commandline);

  /* convert 4-bit floating point values to slightly discretised values of type
     "value" (see Common.h) */
  for (i=0; i<nrnumerics; i++)
    {
      for (j=0; j<NRPAT; j++)
	{
	  feature[ffindex[i]][j]=(value)
	    (((floatfeature[i][j]-ffmin[i])*(1.*DISCRETIZER))/
	     (ffmax[i]-ffmin[i]));
	}
      free(floatfeature[i]);
    }

  /* determine the maximal size that a family will take up */

  famtop=0;
  for (i=0; i<NRCLASSES; i++)
    if (overallclassocc[i]>famtop)
      famtop=overallclassocc[i];
  famtop++;
  famtop=min(famtop,MAXNREXP);

  // compute most frequent class
  topocc=0;
  for (i=0; i<NRCLASSES; i++)
    if (overallclassocc[i]>topocc)
      {
	mostfreqclass=i;
	topocc=overallclassocc[i];
      }

  free(line);
  free(previousline);

}

void presort_data( const char *dataname )
{
  FILE *data;
  int  i,j;
  int  counter,nrfocal=0,nrsuppr=0;
  value code=0,thisfocal=-1;
  char *line;
  char *previousline;
  char in;
  char *thisvalue;
  char delimiter;
  char deltok[10];
  char thisclass[MAXVALLEN+1];
  char commandline[1024];
  char tmpdataname[50];

  // initialize read buffers
  line=(char*)malloc(MAXEXPRESSION*sizeof(char));
  previousline=(char*)malloc(MAXEXPRESSION*sizeof(char));
  strcpy(previousline,"#");

  // determine delimiter
  data=fopen(dataname,"r");
  if ( fgets(line,MAXEXPRESSION,data) == NULL ){
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
  fclose(data);

  if (VERB)
    fprintf(stderr,"     presorting and rereading instance base %s\n",
           dataname);

  /* first, create a shadow data file that is sorted */

  sprintf(tmpdataname,"fambl%d.tmp",(unsigned int) lrand48()%UINT_MAX);

  /* do a sort on the most important feature */
  sprintf(commandline,"nice sort -t'%c' %s -o %s -k %d,%d -T .\n",
	  delimiter,dataname,tmpdataname,featorder[0]+1,featorder[0]+2);
  if (VERB4)
    fprintf(stderr,"     executing commandline:\n     %s",commandline);
  my_system(commandline);

  /* allocate the focal index arrays */
  focalstart=(int*)malloc((nrvalues[featorder[0]]+1)*sizeof(int));
  focalend=(int*)malloc((nrvalues[featorder[0]]+1)*sizeof(int));

  /* readin the data file, while counting types */

  data=fopen(tmpdataname,"r");
  if (VERB2)
    fprintf(stderr,"     file %s opened\n",
	    tmpdataname);

  counter=0;
  if ( fgets(line,MAXEXPRESSION,data) == NULL ){
    fprintf( stderr, "read failure on file %s\n", tmpdataname );
    exit(1);
  }
  while (!feof(data))
    {
      if (VERB)
	if ((counter>0)&&
	    (counter%(REPORTRATE*100)==0))
	  fprintf(stderr,"     %7d instances read\n",
		  counter);

      thisvalue=strtok(line,deltok);
      for (i=0; i<PATWIDTH; i++)
	{

	  code=0;
	  in=0;
	  for (j=0; ((j<nrvalues[i])&&(!in)); j++)
	    if (thisvalue[0]==values[i][j][0])
	      if (strcmp(thisvalue,values[i][j])==0)
		{
		  in=1;
		  code=j;
		}

	  feature[i][counter]=code;

	  if ((i==featorder[0])&&(code!=thisfocal))
	    {
	      if (nrfocal>0)
		focalend[thisfocal]=counter;
	      focalstart[code]=counter;
	      nrfocal++;
	      thisfocal=code;
	    }
	  thisvalue=strtok(NULL,deltok);

	}
      strcpy(thisclass,thisvalue);

      in=0;
      for (j=0; ((j<NRCLASSES)&&(!in)); j++)
	if (thisclass[0]==classes[j][0])
	  if (strcmp(thisclass,classes[j])==0)
	    {
	      in=1;
	      code=j;
	    }

      klass[counter]=code;
      counter++;
      if ( fgets(line,MAXEXPRESSION,data) == NULL && !feof(data )){
	fprintf( stderr, "read failure on file %s\n", tmpdataname );
	exit(1);
      }
    }
  fclose(data);
  focalend[thisfocal]=counter;

  if (VERB4)
    {
      fprintf(stderr,"     index on focal feature:\n");
      for (i=0; i<nrvalues[featorder[0]]; i++)
	fprintf(stderr,"     %2d (%s): %d - %d\n",
		i,values[featorder[0]][i],focalstart[i],focalend[i]);
    }

  if ((VERB)&&(nrsuppr>0))
    fprintf(stderr,"     %d feature value tokens suppressed\n",
	   nrsuppr);

  /* remove shadow data file */

  sprintf(commandline,"rm %s\n",tmpdataname);
  my_system(commandline);

  free(line);
  free(previousline);

}
