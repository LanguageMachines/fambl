/* Metrics.c

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

/* compute_grsum - precompute some useful gain ratio sums, and check
   whether there is a superfeature heavier than the sum of the rest
*/
void compute_grsum(void)
{ int   i,j;
  float othersum;

  maxgrsum=0.0;
  for (i=0; i<PATWIDTH; i++)
    {
      maxgrsum+=featgrvals[i];

      othersum=0.0;
      for (j=0; j<PATWIDTH; j++)
	if (j!=i)
	  othersum+=featgrvals[j];
      if (othersum<featgrvals[i])
	{
	  superfeat=i;
	  thereisasuperfeat=1;
	  if (VERB2)
	    fprintf(stderr,"     super heavy feature (# %d) detected\n",i);
	}
    }

  if (METRIC>0)
    maxgrsum*=2.;

  if (VERB3)
    fprintf(stderr,"     sum of feature weights: %f\n",maxgrsum);
}

/* equalize_gr - flatten feature weights when only MVDM is used
*/
void equalize_gr(void)
{
  int i;

  for (i=0; i<PATWIDTH; i++)
    {
      featgrvals[i]=1.0;
      featorder[i]=i;
    }
  if (ALGORITHM==0)
    maxgrsum=2.*PATWIDTH;
}

/* adapt_weights_to_algorithm - adapt the computed GR weights such that
   they emulate the feature ordering assumed in IGTree (complete
   ordering) or TRIBL (only ordering on the most important features).
*/
void adapt_weights_to_algorithm(void)
{
  int i,j;
  int cutofffeat=0;
  float sum=0.0,summ,mean,vari,cutoff;

  /* in case of IGTree, just start all over and create weights that
     are larger than the sum of those of all less important features */
  if (FWEIGHT==6)
    {
      for (i=PATWIDTH-1; i>=0; i--)
	{
	  sum=1.0;
	  if (i<(PATWIDTH-1))
	    {
	      for (j=i+1; j<PATWIDTH; j++)
		sum+=featgrvals[featorder[j]];
	    }
	  featgrvals[featorder[i]]=sum;
	}
    }
  /* in case of TRIBL, first the cutoff point needs to be found. The
     default heuristic, implemented here, is to take this to be the
     average GR plus one standard deviation. This seems to work quite
     reasonably. Then proceed as above. */
  if (FWEIGHT==7)
    {
      summ=0.0;
      for (i=0; i<PATWIDTH; i++)
	summ+=featgrvals[i];
      mean=summ/(1.*PATWIDTH);
      vari=0.0;
      for (i=0; i<PATWIDTH; i++)
	vari+=((featgrvals[i]-mean)*(featgrvals[i]-mean));
      vari/=(1.*PATWIDTH);
      cutoff=mean+sqrt(vari);
      if (VERB3)
	fprintf(stderr,"     cutoff value found by TRIBL: %f\n",
	       cutoff);
      if (featgrvals[featorder[cutofffeat]]>cutoff)
	while (featgrvals[featorder[cutofffeat]]>cutoff)
	  cutofffeat++;
    if (VERB2)
      fprintf(stderr,"     TRIBL sets cutoff point after %d features\n",
	     cutofffeat);
    if (cutofffeat>0)
      {
	for (i=cutofffeat-1; i>=0; i--)
	  {
	    sum=0.0;
	    if (i<(PATWIDTH-1))
	      {
		for (j=i+1; j<PATWIDTH; j++)
		  sum+=featgrvals[featorder[j]];
	      }
	    featgrvals[featorder[i]]=sum;
	  }
      }
    }
  if (VERB)
    {
      fprintf(stderr,"     adapted GR feature weights:\n");
      for (i=0; i<PATWIDTH; i++)
	fprintf(stderr,"     feature %3d : GR %10.6f\n",i,featgrvals[i]);
    }
}

/* compute_gr - a terrible lot of hassle for a strange information-theory
   based feature weighting heuristic (cf. Quinlan, 1993). Uses the usual method
   for symbolic features, and uses a simple heuristic for numerics: cut the
   data into 20 bins and regard these as the values.
   also compute the top IGs to be obtained.
*/
void compute_gr(void)
{
  int   i,j,k,m;
  value classnr;
  int   nrsubs;
  int   suboccs[MAXNRCLASSES];
  float sum,pocc,dbentropy,splitinfo,average;

  if (VERB)
    fprintf(stderr,"     computing feature gain ratio values\n");

  /* first, compute dbentropy */

  sum=0.0;
  for (classnr=0; classnr<NRCLASSES; classnr++)
  { pocc=((1.*overallclassocc[classnr])/(1.*NRPAT));
    sum+=(pocc*(log((double) pocc)/log(2.0)));
  }
  dbentropy=sum*-1.;
  if (VERB2)
    fprintf(stderr,"     entropy: %7.4f (redundancy %7.4f)\n",
	   dbentropy,1-(dbentropy/((log((double) NRCLASSES)/log(2.0)))));

  if (VERB3)
    fprintf(stderr,"     computing gain ratio values\n");

  average=0.0;
  for (i=0; i<PATWIDTH; i++)
  {
    featgrvals[i]=0.0;
    splitinfo=0.0;
    if ((nrvalues[i]>1)||(feattype[i]=='n'))
      {
	/* symbolic features */
	if (feattype[i]=='s')
	  {
	    for (j=0; j<nrvalues[i]; j++)
	      {
		nrsubs=0;
		for (k=0; k<NRCLASSES; k++)
		  {
		    suboccs[k]=valocc[i][k][j];
		    nrsubs+=valocc[i][k][j];
		  }
		sum=0.0;
		for (m=0; m<NRCLASSES; m++)
		  {
		    if (nrsubs>0)
		      pocc=(1.*suboccs[m])/(1.*nrsubs);
		    else
		      pocc=0.0;
		    if (pocc>0.0)
		      sum+=(pocc*(log((double) pocc)/log(2.0)));
		  }
		if (nrsubs>0)
		  {
		    featgrvals[i]+=
		      ((-1.*sum)*((1.*nrsubs)/(1.*NRPAT)));
		    splitinfo+=(((1.*nrsubs)/NRPAT)*
				(log((double) ((1.*nrsubs)/(1.*NRPAT)))/
				 log(2.0)));
		  }
	      }
	  }
	/* numeric features: discretize numeric continuum in evenly-sized bins
	   (awkward but simple) */
	else
	  {
	    for (j=0; j<BINS; j++)
	      {
		nrsubs=0;
		for (k=0; k<NRCLASSES; k++)
		  suboccs[k]=0;
		nrsubs=0;
		for (k=0; k<NRPAT; k++)
		  {
		    if (((feature[i][k]>=((int)((1.*j*DISCRETIZER)/(1.*BINS))))&&
			 (feature[i][k]<((int)((1.*(j+1)*DISCRETIZER)/(1.*BINS)))))||
			((j==BINS-1)&&
			 (feature[i][k]==DISCRETIZER)))
		      {
			nrsubs++;
			suboccs[klass[k]]++;
		      }
		  }
		sum=0.0;
		for (m=0; m<NRCLASSES; m++)
		  {
		    if (nrsubs>0)
		      pocc=(1.*suboccs[m])/(1.*nrsubs);
		    else
		      pocc=0.0;
		    if (pocc>0.0)
		      sum+=(pocc*(log((double) pocc)/log(2.0)));
		  }
		if (nrsubs>0)
		  {
		    featgrvals[i]+=
		      ((-1.*sum)*((1.*nrsubs)/(1.*NRPAT)));
		    splitinfo+=(((1.*nrsubs)/NRPAT)*
				(log((double) ((1.*nrsubs)/(1.*NRPAT)))/
				 log(2.0)));
		  }
	      }
	  }

	splitinfo*=-1.;
	/* when a feature has one value, or a very major value,
	   splitinfo is (or approaches) 0. In that case, GR
	   can be thought of as 1.0. This isn't math - it's heuristics. */
	if (FWEIGHT>=2)
	  {
	    if (splitinfo<ALMOSTNOTHING)
	      featgrvals[i]=1.0;
	    else
	      featgrvals[i]=
		(dbentropy-featgrvals[i])/splitinfo;
	  }
	else
	  featgrvals[i]=dbentropy-featgrvals[i];
      }
    else
      featgrvals[i]=1.0;

    /* report */
    if (VERB)
      {
	if (feattype[i]=='s')
	  fprintf(stderr,"     feature %3d (%5d values): %10.6f\n",
		  i,
		  nrvalues[i],
		  featgrvals[i]);
	else
	  fprintf(stderr,"     feature %3d (   numeric  ): %10.6f\n",
		  i,
		  featgrvals[i]);
      }
    average+=featgrvals[i];
  }

  /* compute an average */
  average/=(1.*PATWIDTH);
  if (VERB)
    {
      if (FWEIGHT>=2)
	fprintf(stderr,"     average gain ratio: %f\n",
		average);
      else
	fprintf(stderr,"     average information gain: %f\n",
		average);
    }
}

/* motherfeat - determines the mother feature of a value number */
int motherfeat(int valuenr)
{ int i;

  i=0;
  while ((i<PATWIDTH)&&(valuenr>=0))
    {
      valuenr-=nrvalues[i];
      i++;
    }
  return i-1;
}

/* minoffset - determines the mother feature offset of a value number */
int minoffset(int valuenr)
{ int i,prev;

  i=0;
  prev=valuenr;
  while ((i<PATWIDTH)&&(valuenr>=0))
    {
      prev=valuenr;
      valuenr-=nrvalues[i];
      i++;
    }
  return prev;
}

float largel(float p, int k, int n)
{
  float result;

  if ((p<0.00000001)||(p>0.99999999))
    result=0.0;
  else
    result= (float) ( (1.*k)*(log((float) p)) + ((1.*(n-k))*(log(1.-(float) p))) );

  return result;
}

/* atomic_rooster - automatic unpacking of multi-valued features into
   binary (atomic) ones. Additional option: combine atomic features
   into new features up to a combination-cardinality ceiling.
 */
void atomic_rooster(void)
{
  value classnr;
  int   **thisvalocc;
  int   *thistotalocc;
  value **cvalue;
  int   i,j,k,l,m,n=0,nrsubs,sublength,thisnrvalues,thismaxcombination,
    thisbonuses=0,selpointcount,k1,k2,n1,n2;
  float sum=0.0,pocc,dbentropy,chisquare,logl,expected=0.0,splitinfo,
    littlefeatgrval,competeval=0.0,thiscompeteval,p,p1,p2,thislogl,
    hf,kf,thismvdm,mf;
  value thisvalue[MAXPATWIDTH];
  int   selpoint[MAXPATWIDTH];
  char  *suspect[MAXPATWIDTH];
  char  in,ready,pass,superpass,match,used,allsuspect,thissuspect;
  void  increaseselect(int sublength);

  if ((PATWIDTH-1)>NRCOMBINATIONS)
    thismaxcombination=NRCOMBINATIONS;
  else
    thismaxcombination=(PATWIDTH-1);

  thisbonuses=MAXBONUSES;

  /* initialise list of suspect feature values */
  for (i=0; i<PATWIDTH; i++)
    {
      suspect[i]=(char*)malloc(nrvalues[i]*sizeof(char));
      if (suspect[i]==NULL)
	{
	  alloc_error();
	  exit(1);
	}
      for (j=0; j<nrvalues[i]; j++)
	suspect[i][j]=0;
    }

  nrbonuses=0;
  for (i=0; i<thismaxcombination; i++)
    combinationlimit[i]=0;

  /* allocate the temporary buffers */

  thisvalocc=(int**)malloc(NRCLASSES*sizeof(int*));
  if (thisvalocc==NULL)
    alloc_error();
  for (j=0; j<NRCLASSES; j++)
    {
      thisvalocc[j]=(int*)malloc(nrtypes*sizeof(int));
      if (thisvalocc[j]==NULL)
	alloc_error();
    }
  thistotalocc=(int*)malloc(nrtypes*sizeof(int));
  if (thistotalocc==NULL)
    alloc_error();

  selected=(char*)malloc(PATWIDTH*sizeof(char));
  if (selected==NULL)
    alloc_error();

  cvalue=(value**)malloc(nrtypes*sizeof(value*));
  if (cvalue==NULL)
    alloc_error();
  for (i=0; i<nrtypes; i++)
    {
      cvalue[i]=(value*)malloc(thismaxcombination*sizeof(value));
      if (cvalue[i]==NULL)
	alloc_error();
    }

  /* compute dbentropy */
  sum=0.0;
  for (classnr=0; classnr<NRCLASSES; classnr++)
    {
      pocc=((1.*overallclassocc[classnr])/(1.*NRPAT));
      sum+=(pocc*(log((double) pocc)/log(2.0)));
    }
  dbentropy=sum*-1.;
  if (VERB2)
    fprintf(stderr,"     entropy: %7.4f bits (redundancy %7.4f)\n",
	   dbentropy,1-(dbentropy/((log((double) NRCLASSES)/log(2.0)))));

  if (VERB)
    fprintf(stderr,"     automatic unpacking of multi-valued to atomic values\n");

  allsuspect=0;

  for (sublength=1;
       ((!allsuspect)&&(sublength<=thismaxcombination));
       sublength++)
    {
      /* remember that when going to each new sublength, all
	 features with sublength-1 compounding are below the current
	 value of "nrbonuses" */
      combinationlimit[sublength-1]=nrbonuses;

      for (i=0; i<PATWIDTH; i++)
	{
	  if (i<sublength)
	    selected[i]=1;
	  else
	    selected[i]=0;
	}

      ready=0;
      while (!ready)
	{
	  selpointcount=0;
	  for (m=0; m<PATWIDTH; m++)
	    if (selected[m])
	      {
		selpoint[selpointcount]=m;
		selpointcount++;
	      }

	  superpass=1;

	  /* if sublength>2, check if any of the subpart frames of
	     the current selection has held an informative compound.
	     if not, then save the trouble. */
	  if ((superpass)&&(sublength>2))
	    {
	      superpass=0;

	      /* now find out whether subpart frames of the current compound
		 contain informative compounds to begin with */

	      for (l=0; ((!superpass)&&(l<sublength)); l++)
		{
		  for (m=combinationlimit[sublength-2];
		       ((!superpass)&&(m<combinationlimit[sublength-1]));
		       m++)
		    {
		      match=1;
		      for (n=0; ((match)&&(n<sublength)); n++)
			{
			  if (n!=l)
			    if (bonusselections[m][selpoint[n]]==-1)
			      match=0;
			}
		      if (match)
			superpass=1;
		    }
		}
	    }

	  if (superpass)
	    {
	      /* create the combined feature; count everything there is
		 to count */

	      thisnrvalues=0;
	      for (k=0; k<NRPAT; k++)
		{
		  thissuspect=0;
		  for (m=0; m<sublength; m++)
		    {
		      thisvalue[m]=feature[selpoint[m]][k];
		      if (suspect[selpoint[m]][thisvalue[m]])
			thissuspect=1;
		    }
		  if (!thissuspect)
		    {
		      in=0;
		      m=0;
		      while ((!in)&&(m<thisnrvalues))
			{
			  if (cvalue[m][0]==thisvalue[0])
			    {
			      in=1;
			      for (l=1; ((in)&&(l<sublength)); l++)
				if (cvalue[m][l]!=thisvalue[l])
				  in=0;
			    }
			  m++;
			}

		      if (!in)
			{
			  for (l=0; l<sublength; l++)
			    cvalue[thisnrvalues][l]=thisvalue[l];

			  for (j=0; j<NRCLASSES; j++)
			    thisvalocc[j][thisnrvalues]=0;
			  thisvalocc[klass[k]][thisnrvalues]=1;
			  thistotalocc[thisnrvalues]=1;

			  thisnrvalues++;
			}
		      else
			{
			  if (!thissuspect)
			    {
			      thisvalocc[klass[k]][m-1]++;
			      thistotalocc[m-1]++;
			    }
			}
		    }
		}

	      if (VERB2)
		{
		  fprintf(stderr,"\r     fts %4d ",
			 selpoint[0]);
		  if (selpoint[1]%100==0)
		    for (l=1; l<sublength; l++)
		      fprintf(stderr,"%4d ",selpoint[l]);
		}

	      /* start selecting the features. at the atomic and double level,
		 this is straightforward. at the third and further levels,
		 things must be focused to what has a chance of pay-off. */

	      for (k=0; k<thisnrvalues; k++)
		{
		  superpass=1;

		  /* if one of the feature values involved is suspect, skip
		     the whole thing immediately. */

		  if (superpass)
		    for (l=0; ((superpass)&&(l<sublength)); l++)
		      if (suspect[selpoint[l]][cvalue[k][l]])
			superpass=0;

		  if ((superpass)&&(sublength>2))
		    {
		      superpass=0;

		      /* now find out whether subparts of the current compound
			 value have been seen as atomics and compounds just
			 under this level. if not, there's not going to be
			 any pay-off. if there are informative subparts,
			 compute their summed information.  */

		      competeval=0.0;
		      for (l=0; l<sublength; l++)
			{
			  thiscompeteval=0.0;
			  for (m=combinationlimit[sublength-2];
			       m<combinationlimit[sublength-1];
			       m++)
			    {
			      match=1;
			      for (n=0; ((match)&&(n<sublength)); n++)
				if (n!=l)
				  if (bonusselections[m][selpoint[n]]!=cvalue[k][n])
				    match=0;
			      if (match)
				{
				  superpass=1;
				  thiscompeteval=bonuses[m];
				  n=0;
				  while ((bonusselections[n][selpoint[l]]!=
					  cvalue[k][l])&&
					 (n<combinationlimit[1]))
				    n++;

				  if (n<combinationlimit[1])
				    thiscompeteval+=bonuses[n];
				  if (thiscompeteval>competeval)
				    competeval=thiscompeteval;
				}
			    }
			}
		      /* competeval/=(1.*sublength); */
		    }

		  if (superpass)
		    {
		      /* first, compute littlefeatgrval, i.e., the weight
			 of the individual constructed feature value.
			 again, do either IG, GR, chi square, shared variance,
			 or log likelihood, and on top of that, MVDM / JD
                         which are also weights, effectively, with binary
		         features */

		      littlefeatgrval=0.0;

		      if (thistotalocc[k]>0)
			{
			  if ((FWEIGHT==1)||(FWEIGHT==2)||
			      (FWEIGHT==6)||(FWEIGHT==7))
			    {
			      /* the value case */
			      nrsubs=thistotalocc[k];
			      sum=0.0;
			      for (m=0; m<NRCLASSES; m++)
				{
				  if (thisvalocc[m][k]>0)
				    {
				      pocc=0.0;
				      if (nrsubs>0)
					pocc=(1.*thisvalocc[m][k])/(1.*nrsubs);
				      else
					pocc=0.0;
				      sum+=(pocc*(log((double) pocc)/log(2.0)));
				    }
				}
			      if (nrsubs>0)
				littlefeatgrval+=
				  ((-1.*sum)*((1.*nrsubs)/(1.*NRPAT)));

			      /* the non-value case */
			      sum=0.0;
			      nrsubs=NRPAT-thistotalocc[k];
			      for (m=0; m<NRCLASSES; m++)
				{
				  if (overallclassocc[m]>thisvalocc[m][k])
				    {
				      pocc=0.0;
				      if (nrsubs>0)
					pocc=(1.*(overallclassocc[m]-thisvalocc[m][k]))/
					  (1.*nrsubs);
				      else
					pocc=0.0;
				      sum+=(pocc*(log((double) pocc)/log(2.0)));
				    }
				}
			      if (nrsubs>0)
				littlefeatgrval+=
				  ((-1.*sum)*((1.*nrsubs)/(1.*NRPAT)));

			      /* GR selected? Then do split info!! */

			      splitinfo=0.0;
			      if (FWEIGHT>=2)
				{

				  /* the value case */
				  splitinfo+=(((1.*thistotalocc[k])/
					       (1.*NRPAT))*
					      (log((double)
						   ((1.*thistotalocc[k])/
						    (1.*NRPAT)))/
					       log(2.0)));

				  /* the non-value case */
				  splitinfo+=(((1.*(NRPAT-thistotalocc[k]))/
					       (1.*NRPAT))*
					      (log((double)
						   ((1.*(NRPAT-thistotalocc[k]))/
						    (1.*NRPAT)))/
					       log(2.0)));

				  splitinfo*=-1.;
				  littlefeatgrval=
				    (dbentropy-littlefeatgrval)/splitinfo;
				}
			      else
				littlefeatgrval=dbentropy-littlefeatgrval;

			      /* post-correct possible 1-value cases
			      if (NRPAT==thistotalocc[k])
			      littlefeatgrval=1.0; */

			    }

			  /* shared variance (chi square) */
			  if ((FWEIGHT==3)||
			      (FWEIGHT==4))
			    {
			      chisquare=0.0;

			      for (l=0; l<NRCLASSES; l++)
				{
				  /* the value case */
				  expected=
				    ((1.*overallclassocc[l])*(1.*thistotalocc[k]))/(1.*NRPAT);
				  if (expected>0.0)
				    chisquare+=
				      ((((1.*thisvalocc[l][k])-expected)*
					((1.*thisvalocc[l][k])-expected))/
				       expected);

				  /* the non-value case */
				  expected=
				    ((1.*overallclassocc[l])*(1.*(NRPAT-thistotalocc[k])))/
				    (1.*NRPAT);
				  if (expected>0.0)
				    chisquare+=
				      ((((1.*(overallclassocc[l]-thisvalocc[l][k]))-expected)*
					((1.*(overallclassocc[l]-thisvalocc[l][k]))-expected))/
				       expected);
				}

			      if (FWEIGHT==4)
				{ /* compute Cramer's phi - note that degree of freedom
				     is necessarily 1 */
				  chisquare=chisquare/(1.*NRPAT);
				}
			      littlefeatgrval=chisquare;
			    }

			  /* log-likelihood */
			  if (FWEIGHT==5)
			    {
			      logl=0.0;

			      for (l=0; l<NRCLASSES; l++)
				{
				  k1=thisvalocc[l][k];
				  n1=overallclassocc[l];
				  k2=thistotalocc[k]-thisvalocc[l][k];
				  n2=NRPAT-overallclassocc[l];

				  p1=(1.*k1)/(1.*n1);
				  p2=(1.*k2)/(1.*n2);

				  p=(1.*(k1+k2))/(1.*(n1+n2));

				  thislogl=2.*(largel(p1,k1,n1)+
					       largel(p2,k2,n2)-
					       largel(p,k1,n1)-
					       largel(p,k2,n2));

				  logl+=thislogl;
				}

			      littlefeatgrval=logl;
			    }

			  // for MVDM and JD, things multiply with
			  // the weight already computed

			  /* MVDM and JD, with smoothing */
			  if (METRIC)
			    {

			      if (FWEIGHT==0)
				littlefeatgrval=1.0;

			      //fprintf(stderr,"littlefeatgrval: %f\n",
			      //   littlefeatgrval);

			      thismvdm=0.0;
			      if (thistotalocc[k]>=MVDML)
				{

				  for (l=0; l<NRCLASSES; l++)
				    {
				      hf=(1.*thisvalocc[l][k])/(1.*thistotalocc[k]);
				      kf=(1.*(overallclassocc[l]-thisvalocc[l][k]))/(1.*(NRPAT-thistotalocc[k]));
				      if (METRIC==1)
					thismvdm+=
					  fabs(hf-kf);
				      else
					{
					  mf=(hf+kf)/2.;
					  if (mf>0.0)
					    {
					      if (hf>0.0)
						thismvdm+=
						  (hf*(log(hf/mf)/log(2.)));
					      if (kf>0.0)
						thismvdm+=
						  (kf*(log(kf/mf)/log(2.)));
					    }
					}
				    }
				  thismvdm=2.-thismvdm;
				  //if (thismvdm<0.0)
				  //thismvdm=0.0;
				}
			      else
				{
				  thismvdm=0.0;
				}

			      littlefeatgrval*=thismvdm;

			      //fprintf(stderr,"times thismvdm %f gives %f\n",
			      //     thismvdm,littlefeatgrval);
			    }
			}
		      /* then if the feature outweighs the occurrence threshold,
			 go add it to the list of bonuses */

		      pass=1;

		      /* the "2" case is just in between. we'll have
			 to check the competition here of the sum of its
			 atomic parts */
		      if (sublength==2)
			{
			  competeval=0.0;
			  for (l=0; l<2; l++)
			    {
			      n=0;
			      while ((bonusselections[n][selpoint[l]]!=
				      cvalue[k][l])&&
				     (n<combinationlimit[1]))
				n++;
			      if (n<combinationlimit[1])
				competeval+=bonuses[n];
			    }
			  /* competeval/=2.0; */
			}

		      /* here's the critical comparison that will accept
			 only extra-informative new compound features! */

		      if (sublength>=2)
			{
			  pass=0;
			  if (littlefeatgrval>(competeval+ALMOSTNOTHING))
			    pass=1;
			}

		      if ((pass)&&
			  (nrbonuses<thisbonuses))
			{
			  if (sublength>=2)
			    littlefeatgrval-=competeval;

			  if (VERB2)
			    {
			      fprintf(stderr,"\r     feature %6d, occ%6d, wgt %12.9f:",
				      nrbonuses,
				      thistotalocc[k],
				      littlefeatgrval);
			      m=0;
			      for (l=0; l<PATWIDTH; l++)
				if (selected[l])
				  {
				    fprintf(stderr," %s (%d) ",
					    values[l][cvalue[k][m]],l);
				    m++;
				  }
			      fprintf(stderr,"\n");
			    }

			  /* store the newly found bonus */
			  bonusselections[nrbonuses]=
			    (value*)malloc(PATWIDTH*sizeof(value));
			  if (bonusselections[nrbonuses]==NULL)
			    alloc_error();

			  m=0;
			  for (l=0; l<PATWIDTH; l++)
			    {
			      if (selected[l])
				{
				  bonusselections[nrbonuses][l]=cvalue[k][m];
				  m++;
				}
			      else
				bonusselections[nrbonuses][l]=-1;
			    }

			  bonuses[nrbonuses]=littlefeatgrval;

			  if (nrbonuses<thisbonuses)
			    nrbonuses++;
			}
		    }
		}
	    }
	  /* are we done yet with searching through all combinations? */
	  ready=1;
	  for (i=PATWIDTH-1; ((ready)&&(i>(PATWIDTH-1-sublength))); i--)
	    if (selected[i]!=1)
	      ready=0;
	  if (!ready)
	    increaseselect(sublength);
	}

      /* update the suspect list */
      for (m=0; m<PATWIDTH; m++)
	for (n=0; n<nrvalues[m]; n++)
	  {
	    if (!suspect[m][n])
	      {
		used=0;
		for (l=combinationlimit[sublength-1];
		     ((!used)&&(l<nrbonuses));
		     l++)
		  {
		    if (bonusselections[l][m]==n)
		      used=1;
		  }
		if (!used)
		  {
		    suspect[m][n]=1;
		    if (VERB3)
		      fprintf(stderr,"\r     discarded val %5d of feat %3d (%s)\n",
			     n,m,values[m][n]);
		  }
	      }
	  }
      /* is everything suspect? bail out! */
      allsuspect=1;
      for (m=0; ((allsuspect)&&(m<PATWIDTH)); m++)
	for (n=0; ((allsuspect)&&(n<nrvalues[m])); n++)
	  if (suspect[m][n]==0)
	    allsuspect=0;
      if ((VERB2)&&(allsuspect))
	fprintf(stderr,"\r     no more extra information in more compound features to be expected\n");

      if (VERB)
	{ if ((sublength>1)&&
	      (nrbonuses-combinationlimit[sublength-1]>0))
	  fprintf(stderr,"\r     %6d compound features made up of %d values\n",
		 nrbonuses-combinationlimit[sublength-1],sublength);
	else
	  if ((sublength==1)&&
	      (VERB2))
	    fprintf(stderr,"\r     %6d atomic-value features\n",
		   nrbonuses);
	}

    }

  if (VERB)
  {
    if (combinationlimit[1]>0)
      fprintf(stderr,"\r     %d features found, of which %d compound\n",
	     nrbonuses,nrbonuses-combinationlimit[1]);
    else
      fprintf(stderr,"\r     %d binary features found\n",
	     nrbonuses);
  }

  /* free the temporary buffers */

  for (i=0; i<nrtypes; i++)
    free(cvalue[i]);
  free(cvalue);
  free(selected);
  free(thistotalocc);
  free(thisvalocc);
  for (i=0; i<PATWIDTH; i++)
    free(suspect[i]);

  /* allocate the applicable bonuses buffer */
  applicable_bonuses=(int*)malloc(nrbonuses*sizeof(int));
  if (applicable_bonuses==NULL)
    {
      alloc_error();
      exit(1);
    }

}

/* increaseselect - stoopid function for increasing a bit array */
void increaseselect(int sublength)
{
  int i,j;
  int  nrones;

  /* if the end of the line is filled with 1s but not all, shift
     the whole bunch back */
  nrones=0;
  i=PATWIDTH-1;
  while (selected[i])
  { if (selected[i])
      nrones++;
    i--;
  }
  if ((nrones>0)&&(nrones<sublength))
  { for (i=PATWIDTH-1-nrones; i<PATWIDTH; i++)
      selected[i]=0;
    i=PATWIDTH-1;
    while ((!selected[i])&&(i>0))
      i--;
    selected[i]=0;
    i++;
    for (j=i; j<i+nrones+1; j++)
      selected[j]=1;
  }
  else
  { /* go for the next zero */
    i=PATWIDTH-1;
    while (selected[i])
      i--;

    /* go for the next one */
    while (!selected[i])
      i--;

    /* if this one is followed by a zero, shift it */
    if (i<PATWIDTH-1)
      if ((selected[i])&&(!selected[i+1]))
      { selected[i+1]=1;
        selected[i]=0;
      }
  }

  /* fprintf(stderr,"selected: ");
  for (i=0; i<PATWIDTH; i++)
    fprintf(stderr,"%d",selected[i]);
  fprintf(stderr,"\n"); */
}

/* optimize_with_gr: when GR is used, compute the ordering of features
   and restvalues (sums of GRs of less important features) for use
   later on in k-NN to speed up things */
void optimize_with_gr(void)
{
  int   i,j;
  int   maxnr=0;
  float max;
  int   checked[MAXPATWIDTH];

  if (METRIC==0)
    {
      /* compute feature order for optimised k-nn later on */
      for (i=0; i<PATWIDTH; i++)
	checked[i]=0;
      for (i=0; i<PATWIDTH; i++)
	{
	  max=0.0;
	  for (j=0; j<PATWIDTH; j++)
	    {
	      if ((!checked[j])&&(featgrvals[j]>max))
		{
		  max=featgrvals[j];
		  maxnr=j;
		}
	    }
	  featorder[i]=maxnr;
	  checked[maxnr]=1;
	}
      /* also compute rest values for optimised k-nn later on */
      for (i=0; i<PATWIDTH; i++)
	{
	  restvals[i]=0.0;
	  for (j=i; j<PATWIDTH; j++)
	    restvals[i]+=featgrvals[featorder[j]];
	}
      restvals[PATWIDTH]=0.0;
    }
  else
    {
      /* compute feature order */
      for (i=0; i<PATWIDTH; i++)
	checked[i]=0;
      for (i=0; i<PATWIDTH; i++)
	{
	  max=0.0;
	  for (j=0; j<PATWIDTH; j++)
	    {
	      if ((!checked[j])&&((2.*featgrvals[j])>max))
		{
		  max=(2.*featgrvals[j]);
		  maxnr=j;
		}
	    }
	  featorder[i]=maxnr;
	  checked[maxnr]=1;
	}
      /* compute remaining GR sums per feature-ordered position */
      for (i=0; i<PATWIDTH; i++)
	{
	  restvals[i]=0.0;
	  for (j=i; j<PATWIDTH; j++)
	    restvals[i]+=(2.*featgrvals[featorder[j]]);
	}
      restvals[PATWIDTH]=0.0;
    }
}

/* count_vcocc - count value-class coocurrences, for use
   in gain ratio and vdm computations later on. also
   suppresses low-frequent values, if selected (-z option) -
   thanks to Jakub Zavrel */
void count_vcocc(void)
{
  int i,j,k;

  if (VERB2)
    fprintf(stderr,"     counting value-class coocurrences\n");

  /* build matrices feature by feature */

  if (totalocc==NULL)
    totalocc=(int**)calloc(PATWIDTH,sizeof(int*));
  if (totalocc==NULL)
    alloc_error();
  for (i=0; i<PATWIDTH; i++)
    if (feattype[i]=='s')
    {
      if (VERB2)
        fprintf(stderr,"     creating value-class matrix of feature %d\n",
	       i);

      /* allocate for each class a new value occurrence buffer */
      if (valocc[i]==NULL)
	valocc[i]=(int**)calloc(NRCLASSES,sizeof(int*));
      if (valocc[i]==NULL)
	alloc_error();

      /* first, count raw occurrences (values per class and values total) */

      for (j=0; j<NRCLASSES; j++)
	{
	  if (valocc[i][j]==NULL)
	    valocc[i][j]=(int*)calloc((nrvalues[i]+1),sizeof(int));
	  if (valocc[i][j]==NULL)
	    alloc_error();
	  for (k=0; k<nrvalues[i]; k++)
	    valocc[i][j][k]=0;
	}

      if (totalocc[i]==NULL)
	totalocc[i]=(int*)calloc((nrvalues[i]+1),sizeof(int));
      if (totalocc[i]==NULL)
	alloc_error();
      for (j=0; j<nrvalues[i]; j++)
	totalocc[i][j]=0;

      for (j=0; j<NRPAT; j++)
	{
	  valocc[i][klass[j]][feature[i][j]]++;
	  totalocc[i][feature[i][j]]++;
	}
    }
}

/* compute_mvdm - compute Cost and Salzberg's (1993) modified value
   difference metric, based on Stanfill and Waltz' (1986) value
   difference metric (and not to be confused with Domingos' (1995)
   simplified value difference metric)
*/
void compute_mvdm(void)
{
  int   i=0,j,k,l,nrtop,rank,indexj;
  float hf,kf,m;
  float thismvdm;
  int   topfreq[TOPMVDM+1];
  value topvalue[TOPMVDM+1];

  if (VERB)
    {
      if (METRIC==1)
	fprintf(stderr,"     computing MVDM matrices\n");
      if (METRIC==2)
	fprintf(stderr,"     computing Jeffrey divergence matrices\n");
    }

  // first, initialize top-100 (max) indices of most frequent values
  // per feature for which MVDM is going to be prestored.

  nrmvdm=(int*)malloc((PATWIDTH+1)*sizeof(int));
  if (nrmvdm==NULL)
    alloc_error();
  mvdmindex=(value**)malloc((PATWIDTH+1)*sizeof(value*));
  if (mvdmindex==NULL)
    alloc_error();
  invertedmvdmindex=(value**)malloc((PATWIDTH+1)*sizeof(value*));
  if (invertedmvdmindex==NULL)
    alloc_error();
  mvdmbotfreq=(int*)malloc((PATWIDTH+1)*sizeof(int));
  if (mvdmbotfreq==NULL)
    alloc_error();
  mvdm=(float***)malloc((PATWIDTH+1)*sizeof(float**));
  if (mvdm==NULL)
    alloc_error();

  for (i=0; i<PATWIDTH; i++)
    {
      if (VERB2)
	fprintf(stderr,"     finding top frequent values of feature %d\n",
		i);

      if (nrvalues[i]<TOPMVDM)
	nrmvdm[i]=nrvalues[i];
      else
	nrmvdm[i]=TOPMVDM;
      mvdmindex[i]=(value*)malloc((nrmvdm[i]+1)*sizeof(value));
      invertedmvdmindex[i]=(value*)malloc((nrvalues[i]+1)*sizeof(value));

      nrtop=1;
      topfreq[0]=totalocc[i][0];
      topvalue[0]=0;
      for (j=1; j<nrvalues[i]; j++)
	{
	  rank=0;
	  while ((rank<nrtop)&&
		 (totalocc[i][j]<topfreq[rank]))
	    rank++;
	  if (rank<nrmvdm[i])
	    {
	      if (nrtop>rank)
		for (k=nrtop; k>rank; k--)
		  {
		    topfreq[k]=topfreq[k-1];
		    topvalue[k]=topvalue[k-1];
		  }
	      topfreq[rank]=totalocc[i][j];
	      topvalue[rank]=j;
	      if (nrtop<nrmvdm[i])
		nrtop++;
	    }
	}
      if (VERB3)
	{
	  fprintf(stderr,"     most frequent values:\n");
	  for (j=0; j<nrmvdm[i]; j++)
	    fprintf(stderr,"     %3d. %9d - %s\n",
		    j,topfreq[j],values[i][topvalue[j]]);
	}
      for (j=0; j<nrmvdm[i]; j++)
	{
	  mvdmindex[i][j]=topvalue[j];
	  invertedmvdmindex[i][topvalue[j]]=j;
	}
      mvdmbotfreq[i]=topfreq[nrmvdm[i]-1];
    }

  /* patch: featgrvals must be initialised with automatic feature
     selection, since they are not used anywhere else with that
     option */
  if (ALGORITHM)
    for (i=0; i<PATWIDTH; i++)
      featgrvals[i]=1.0;

  /* build matrices feature by feature */

  for (i=0; i<PATWIDTH; i++)
    if (feattype[i]=='s')
      {
	if (VERB2)
	  fprintf(stderr,"     computing top-%d matrix for feature %d\n",
		  nrmvdm[i],i);

	mvdm[i]=(float**)malloc((nrmvdm[i]+1)*sizeof(float*));
	if (mvdm[i]==NULL)
	  alloc_error();

	for (j=0; j<nrmvdm[i]; j++)
	  {
	    mvdm[i][j]=(float*)malloc((nrvalues[i]+1)*sizeof(float));
	    indexj=mvdmindex[i][j];
	    for (k=0; k<nrvalues[i]; k++)
	      {
		mvdm[i][j][k]=0.0;
		if ((totalocc[i][indexj]>=MVDML)&&
		    (totalocc[i][k]>=MVDML))
		  {
		    thismvdm=0.0;
		    for (l=0; l<NRCLASSES; l++)
		      {
			hf=(1.*valocc[i][l][indexj])/(1.*totalocc[i][indexj]);
			kf=(1.*valocc[i][l][k])/(1.*totalocc[i][k]);
			if (METRIC==1)
			  thismvdm+=
			    fabs(hf-kf);
			else
			  {
			    m=(hf+kf)/2.;
			    if (m>0.0)
			      {
				if (hf>0.0)
				  thismvdm+=
				    (hf*(log(hf/m)/log(2.)));
				if (kf>0.0)
				  thismvdm+=
				    (kf*(log(kf/m)/log(2.)));
			      }
			  }
		      }
		    mvdm[i][j][k]=thismvdm;
		    //mvdm[i][j][k]=(2.-thismvdm)*featgrvals[i];
		  }
		if (VERB4)
		  fprintf(stderr,"     similarity at feature %d between %s and %s: %6.4f\n",
			 i,values[i][indexj],values[i][k],mvdm[i][j][k]);
	      }
	  }
      }
}

/* compute_chisquare - compute either the chi square, or the
   shared variance from each feature's chi-square value, on the
   basis of its value-class table.
   */
void compute_chisquare(void)
{
  register int j,k;
  int i,freedom;
  float chisquare,expected;
  int   thistotalocc[BINS];
  int   thisvalocc[MAXNRCLASSES][BINS];

  if (VERB)
    {
      if (FWEIGHT==3)
	fprintf(stderr,"     computing chi-square feature weights\n");
      if (FWEIGHT==4)
	fprintf(stderr,"     computing shared variance feature weights\n");
    }

  for (i=0; i<PATWIDTH; i++)
    {
      if (feattype[i]=='s')
	{
	  chisquare=0.0;
	  for (j=0; j<NRCLASSES; j++)
	    {
	      for (k=0; k<nrvalues[i]; k++)
		{
		  expected=((1.*overallclassocc[j])*(1.*totalocc[i][k]))/(1.*NRPAT);
		  if (expected>0.0)
		    chisquare+=
		      ((((1.*valocc[i][j][k])-expected)*((1.*valocc[i][j][k])-expected))/
		       expected);
		}
	    }
	  if (FWEIGHT==4)
	    {
	      /* compute Cramer's phi */
	      if (NRCLASSES<nrvalues[i])
		freedom=NRCLASSES;
	      else
		freedom=nrvalues[i];
	      if (freedom>1)
		chisquare=chisquare/(1.*(NRPAT*(freedom-1)));
	      else
		chisquare=0.0;
	    }

	  featgrvals[i]=chisquare;
	  if (VERB)
	    {
	      if (FWEIGHT==3)
		fprintf(stderr,"     feature %3d (%5d values): chi-square %10.6f\n",
		       i,
		       nrvalues[i],
		       featgrvals[i]);
	      if (FWEIGHT==4)
		fprintf(stderr,"     feature %3d (%5d values): shared variance %10.6f\n",
		       i,
		       nrvalues[i],
		       featgrvals[i]);
	    }
	}
      /* numeric features: discretize numeric continuum in evenly-sized bins
	 (awkward but simple) */
      else
	{
	  for (j=0; j<BINS; j++)
	    {
	      thistotalocc[j]=0;
	      for (k=0; k<NRCLASSES; k++)
		thisvalocc[k][j]=0;
	      for (k=0; k<NRPAT; k++)
		{
		  if (((feature[i][k]>=((int)((1.*j*DISCRETIZER)/(1.*BINS))))&&
		       (feature[i][k]<((int)((1.*(j+1)*DISCRETIZER)/(1.*BINS)))))||
		      ((j==BINS-1)&&
		       (feature[i][k]==DISCRETIZER)))
		    {
		      thistotalocc[j]++;
		      thisvalocc[klass[k]][j]++;
		    }
		}
	    }
	  chisquare=0.0;
	  for (j=0; j<NRCLASSES; j++)
	    {
	      for (k=0; k<BINS; k++)
		{
		  expected=((1.*overallclassocc[j])*(1.*thistotalocc[k]))/(1.*NRPAT);
		  if (expected>0.0)
		    {
		      chisquare+=
			((((1.*thisvalocc[j][k])-expected)*((1.*thisvalocc[j][k])-expected))/
			 expected);
		    }
		}
	    }

	  if (FWEIGHT==4)
	    {
	      /* compute Cramer's phi */
	      if (NRCLASSES<BINS)
		freedom=NRCLASSES;
	      else
		freedom=BINS;
	      if (freedom>1)
		chisquare=chisquare/(1.*(NRPAT*(freedom-1)));
	      else
		chisquare=0.0;
	    }

	  featgrvals[i]=chisquare;
	  if (VERB)
	    {
	      if (FWEIGHT==3)
		fprintf(stderr,"     feature %3d (   numeric  ): chi-square %10.6f\n",
		       i,
		       featgrvals[i]);
	      else
		fprintf(stderr,"     feature %3d (   numeric  ): shared variance %10.6f\n",
		       i,
		       featgrvals[i]);
	    }
	}
    }
}


/* compute_logl - compute log likelihood of feature, which is the sum
                  of the sum of value-class log likelihoods
   */
void compute_logl(void)
{
  int   i,j,k,k1,k2,n1,n2;
  float p,p1,p2,thislogl,logl,max;
  int   thistotalocc[BINS];
  int   thisvalocc[MAXNRCLASSES][BINS];

  if (VERB)
    fprintf(stderr,"     computing log-likelihood feature weights\n");

  for (i=0; i<PATWIDTH; i++)
    {
      if (feattype[i]=='s')
	{
	  logl=0.0;
	  for (j=0; j<NRCLASSES; j++)
	    {
	      for (k=0; k<nrvalues[i]; k++)
		{
		  k1=valocc[i][j][k];
		  n1=overallclassocc[j];
		  k2=totalocc[i][k]-valocc[i][j][k];
		  n2=NRPAT-overallclassocc[j];

		  p1=(1.*k1)/(1.*n1);
		  p2=(1.*k2)/(1.*n2);

		  p=(1.*(k1+k2))/(1.*(n1+n2));

		  thislogl=2.*(largel(p1,k1,n1)+
			       largel(p2,k2,n2)-
			       largel(p,k1,n1)-
			       largel(p,k2,n2));

		  logl+=thislogl;
		}
	    }

	  featgrvals[i]=logl;
	  if (VERB2)
	    fprintf(stderr,"     feature %3d (%5d values): raw log-l %10.6f\n",
		   i,
		   nrvalues[i],
		   featgrvals[i]);
	}
      /* numeric features: discretize numeric continuum in evenly-sized bins
	 (awkward but simple) */
      else
	{
	  for (j=0; j<BINS; j++)
	    {
	      thistotalocc[j]=0;
	      for (k=0; k<NRCLASSES; k++)
		thisvalocc[k][j]=0;
	      for (k=0; k<NRPAT; k++)
		{
		  if (((feature[i][k]>=((int)((1.*j*DISCRETIZER)/(1.*BINS))))&&
		       (feature[i][k]<((int)((1.*(j+1)*DISCRETIZER)/(1.*BINS)))))||
		      ((j==BINS-1)&&
		       (feature[i][k]==DISCRETIZER)))
		    {
		      thistotalocc[j]++;
		      thisvalocc[klass[k]][j]++;
		    }
		}
	    }
	  logl=0.0;

	  for (j=0; j<NRCLASSES; j++)
	    {
	      for (k=0; k<BINS; k++)
		{
		  k1=thisvalocc[j][k];
		  n1=overallclassocc[j];
		  k2=thistotalocc[k]-thisvalocc[j][k];
		  n2=NRPAT-overallclassocc[j];

		  p1=(1.*k1)/(1.*n1);
		  p2=(1.*k2)/(1.*n2);

		  p=(1.*(k1+k2))/(1.*(n1+n2));

		  thislogl=2.*(largel(p1,k1,n1)+
			       largel(p2,k2,n2)-
			       largel(p,k1,n1)-
			       largel(p,k2,n2));

		  logl+=thislogl;
		}
	    }

	  featgrvals[i]=logl;
	  if (VERB2)
	    fprintf(stderr,"     feature %3d (   numeric  ): raw log-l %10.6f\n",
		   i,
		   featgrvals[i]);
	}
    }

  if (ALGORITHM==0)
    {
      /* normalize to 0.0 - 1.0 */
      max=0.0;
      for (i=0; i<PATWIDTH; i++)
	if (featgrvals[i]>max)
	  max=featgrvals[i];
      for (i=0; i<PATWIDTH; i++)
	featgrvals[i]/=max;
      if (VERB)
	for (i=0; i<PATWIDTH; i++)
	  {
	    if (feattype[i]=='s')
	      fprintf(stderr,"     feature %3d (%5d values): norm. log-likelihood %10.6f\n",
		     i,
		     nrvalues[i],
		     featgrvals[i]);
	    else
	      fprintf(stderr,"     feature %3d (   numeric  ): norm. log-likelihood %10.6f\n",
		     i,
		     featgrvals[i]);

	  }
    }
  else
    if (VERB)
      for (i=0; i<PATWIDTH; i++)
	{
	  if (feattype[i]=='s')
	    fprintf(stderr,"     feature %3d (%5d values): log-likelihood %10.6f\n",
		   i,
		   nrvalues[i],
		   featgrvals[i]);
	  else
	    fprintf(stderr,"     feature %3d (   numeric  ): log-likelihood %10.6f\n",
		   i,
		   featgrvals[i]);

	}

}
