/** Common.c

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
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <errno.h>

#include "Common.h"

/* alloc_error - give a general message than we've run out of memory
 */
void alloc_error(void)
{ fprintf(stderr," <!> allocation error: not enough memory.\n\n");
  exit(1);
}

/* my_system - run a csh command. adopted from the system() man page
 */
void my_system (char *command)
{
  int sysresult = system( command );
  if ((sysresult==127)||(sysresult==-1))
    {
      fprintf(stderr," <!> Fambl system(%s) failed- aborting.\n\n", command );
      exit(0);
    }
}


/* numdiff - return a floating point difference between two
   value-discretized values */
float  numdiff(value jut, value jul)
{
  return 1.-(fabs((double) ((jut-jul)/(1.*DISCRETIZER))));
}

/* numscale - return a floating point value of a discretized value */
float numscale(value jut)
{
  return (1.*jut)/(1.*DISCRETIZER);
}

/* sweep_mvdm - when using MVDM, just compute the total sum
   of feature-value similarities between the test instance (inst) and a
   memory instance (sweeper). Check whether it has a chance of
   surpassing the lowest distance in the k-NN ranking, otherwise stop
   summing immediately.
*/
void sweep_mvdm(int sweeper)
{
  int   i,j,l;
  value thiss,that;
  float thismvdm,hf,kf,m;

  identical=0.0;

  for (i=0;
       (((identical+restvals[i])>mindist-ALMOSTNOTHING)&&
	(i<PATWIDTH));
       i++)
    {
      j=featorder[i];
      if (feattype[j]=='s')
	{
	  if (feature[j][sweeper]==feature[j][inst])
	    identical+=2.*featgrvals[j];
	  else
	    {
	      thiss=feature[j][sweeper];
	      that=feature[j][inst];
	      thismvdm=0.0;
	      if ((totalocc[j][thiss]>=MVDML)&&
		  (totalocc[j][that]>=MVDML))
		{
		  if (totalocc[j][thiss]>=mvdmbotfreq[j])
		    thismvdm+=
		      mvdm[j][invertedmvdmindex[j][thiss]][that];
		  else
		    {
		      if (totalocc[j][that]>=mvdmbotfreq[j])
			thismvdm+=
			  mvdm[j][invertedmvdmindex[j][that]][thiss];
		      else
			{
			  for (l=0; l<NRCLASSES; l++)
			    {
			      hf=((1.*valocc[j][l][thiss])/
				  (1.*totalocc[j][thiss]));
			      kf=((1.*valocc[j][l][that])/
				  (1.*totalocc[j][that]));

			      if (METRIC==1)
				thismvdm+=fabs(hf-kf);
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
			}
		    }
		}
	      else
		thismvdm=2.;
	      thismvdm=(2.-thismvdm)*featgrvals[j];
	      identical+=thismvdm;
	    }
	}
      else
	identical+=
	  (numdiff(feature[j][sweeper],feature[j][inst])*
	   featgrvals[j]);
    }
}

/* sweep_nomvdm - compute similarity with next instance using only
   GR or flat weighting
*/
void sweep_nomvdm(int sweeper)
{
  register int i,j;
  float thr;

  identical=0.0;
  circles=0;
  thr=mindist-ALMOSTNOTHING;

  for (i=0;
       (((identical+restvals[i])>thr)&&
	(i<PATWIDTH));
       i++)
    {
      j=featorder[i];
      if (feattype[j]=='s')
	{
	  if (feature[j][sweeper]==feature[j][inst])
	    identical+=featgrvals[j];
	  else
	    circles++;
	}
      else
	identical+=
	  (numdiff(feature[j][sweeper],feature[j][inst])*
	   featgrvals[j]);
    }
}

/* sweep_nomvdm_fs - compute similarity via unpacked feature values */
void sweep_nomvdm_fs(int sweeper)
{
  register int   i,j;
  float thr;
  char  match;

  identical=0.0;
  circles=0;
  thr=mindist-ALMOSTNOTHING;

  for (i=0;
       (((identical+bonusrestvals[orderedbonus[i]])>thr)&&
	(i<nrapplicables));
       i++)
    {
      match=1;
      for (j=0; ((match)&&(j<PATWIDTH)); j++)
	if (bonusselections[applicable_bonuses[orderedbonus[i]]][j]!=-1)
	  if (bonusselections[applicable_bonuses[orderedbonus[i]]][j]!=feature[j][sweeper])
	    match=0;
      if (match)
	{
	  identical+=bonuses[applicable_bonuses[orderedbonus[i]]];
	  if (VERB4)
	    {
	      fprintf(stderr,"     applied bonus    %3d:",
		     i);
	      for (j=0; j<PATWIDTH; j++)
		{
		  if (bonusselections[applicable_bonuses[orderedbonus[i]]][j]==-1)
		    fprintf(stderr," *");
		  else
		    fprintf(stderr," %s",
			   values[j][bonusselections[applicable_bonuses[orderedbonus[i]]][j]]);
		}
	      fprintf(stderr,"\n");
	    }
	}
      else
	circles++;
      /* fprintf(stderr," %dth bonus, identical %f, restvals %f, mindist %f\n",
	 i,identical,bonusrestvals[i],mindist); */
    }

}

/* sweep_c_nomvdm_fs - compute similarity via automatically selected
   features */
void sweep_c_nomvdm_fs(int sweeper)
{
  register int  i,j,k;
  char match,among;
  float thr;

  identical=0.0;
  thr=mindist-ALMOSTNOTHING;

  for (i=0;
       (((identical+bonusrestvals[orderedbonus[i]])>thr)&&
	(i<nrapplicables));
	i++)
    {
      match=1;
      for (j=0; ((match)&&(j<PATWIDTH)); j++)
	{
	  if (bonusselections[applicable_bonuses[orderedbonus[i]]][j]!=-1)
	    {
	      if (FAMKA==0)
		{
		  among=0;
		  if (bonusselections[applicable_bonuses[orderedbonus[i]]][j]==
		      expression[sweeper][expression[sweeper][j]])
		    among=1;
		  if (!among)
		    match=0;
		}
	      else
		{
		  among=0;
		  for (k=expression[sweeper][j];
		       k<expression[sweeper][j+1];
		       k++)
		    {
		      if ((bonusselections[applicable_bonuses[orderedbonus[i]]][j]==
			   expression[sweeper][k])||
			  (expression[sweeper][k]==WILDCARDVALUE))
			among=1;
		    }
		  if (!among)
		    match=0;
		}
	    }
	}
      if (match)
	{
	  identical+=bonuses[applicable_bonuses[orderedbonus[i]]];
	}
    }
  thisoccurrence=expressionoccurrence[sweeper];
}

/* sweep_c_mvdm - compute similarity with next instance using MVDM plus
   GR or flat weighting.
 */
void sweep_c_mvdm(int sweeper)
{
  int   i,j,l,k;
  float thismvdm,hf,kf,m;
  float closest;
  value thiss,that;

  identical=0.0;

  for (i=0;
       ((i<PATWIDTH)&&
	((identical+restvals[i])>mindist-ALMOSTNOTHING));
       i++)
    {
      k=featorder[i];
      if (feattype[k]=='s')
	{
	  if (expression[sweeper][expression[sweeper][k]]==WILDCARDVALUE)
	    {
	      identical+=(2.*featgrvals[k]);
	    }
	  else
	    {
	      if (thispattern[k]>-1)
		{
		  thiss=thispattern[k];
		  closest=0.0;
		  for (j=expression[sweeper][k];
		       j<expression[sweeper][k+1];
		       j++)
		    {
		      that=expression[sweeper][j];
		      /* compute MVD values on the fly */
		      thismvdm=0.0;
		      if ((totalocc[k][thiss]>=MVDML)&&
			  (totalocc[k][that]>=MVDML))
			{
			  if (totalocc[k][thiss]>=mvdmbotfreq[k])
			    {
			      /* fprintf(stderr,"taking from mvdm matrix: %d is inverted %d, mvdm %f\n",
				      thiss,invertedmvdmindex[k][thiss],
				      mvdm[k][invertedmvdmindex[k][thiss]][that]); */
			    thismvdm+=
			      mvdm[k][invertedmvdmindex[k][thiss]][that];
			    }
			  else
			    {
			      if (totalocc[k][that]>=mvdmbotfreq[k])
				thismvdm+=
				  mvdm[k][invertedmvdmindex[k][that]][thiss];
			      else
				{
				  for (l=0; l<NRCLASSES; l++)
				    {
				      hf=((1.*valocc[k][l][thiss])/
					  (1.*totalocc[k][thiss]));
				      kf=((1.*valocc[k][l][that])/
					  (1.*totalocc[k][that]));
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
				}
			    }
			}
		      else
			{
			  thismvdm=2.;
			}
		      thismvdm=(2.-thismvdm)*featgrvals[k];
		      if (thismvdm>closest)
			closest=thismvdm;
		    }
		  identical+=closest;
		}
	    }
	}
      else
	{
	  if (expression[sweeper][k]==expression[sweeper][k+1]-1)
	    {
	      identical+=(featgrvals[k]*
			  (numdiff(expression[sweeper][expression[sweeper][k]],
				   thispattern[k])));
	    }
	  else
	    {
	      if ((thispattern[k]>=expression[sweeper][expression[sweeper][k]])&&
		  (thispattern[k]<=expression[sweeper][expression[sweeper][k]+1]))
		{
		  identical+=featgrvals[k];
		}
	      else
		{
		  if (thispattern[k]<expression[sweeper][expression[sweeper][k]])
		    {
		      identical+=(featgrvals[k]*
				  (numdiff(expression[sweeper][expression[sweeper][k]],
					   thispattern[k])));
		    }
		  else
		    {
		      identical+=(featgrvals[k]*
				  (numdiff(expression[sweeper][expression[sweeper][k]+1],
					   thispattern[k])));
		    }
		}
	    }
	}
    }
  thisoccurrence=expressionoccurrence[sweeper];
}

/* sweep_c_nomvdm - compute similarity with next instance using only
   GR or flat weighting
 */
void sweep_c_nomvdm(int sweeper)
{
  register int i,j,k;
  char in;

  identical=0.0;
  in=0;

  for (i=0;
       (((identical+restvals[i])>mindist-ALMOSTNOTHING)&&
	(i<PATWIDTH));
       i++)
    {
      k=featorder[i];
      if (feattype[k]=='s')
	{
	  if (FAMKA==0)
	    {
	      if (expression[sweeper][expression[sweeper][k]]==thispattern[k])
		{
		  identical+=featgrvals[k];
		  in=1;
		}
	    }
	  else
	    {
	      for (j=expression[sweeper][k];
		   j<expression[sweeper][k+1];
		   j++)
		{
		  if ((expression[sweeper][j]==thispattern[k])||
		      (expression[sweeper][j]==WILDCARDVALUE))
		    {
		      identical+=featgrvals[k];
		      in=1;
		    }
		}
	    }
	}
      else
	{
	  if (expression[sweeper][k]==expression[sweeper][k+1]-1)
	    {
	      identical+=(featgrvals[k]*
			  (numdiff(expression[sweeper][expression[sweeper][k]],
				   thispattern[k])));
	      in=1;
	    }
	  else
	    {
	      if ((thispattern[k]>=expression[sweeper][expression[sweeper][k]])&&
		  (thispattern[k]<=expression[sweeper][expression[sweeper][k]+1]))
		{
		  identical+=featgrvals[k];
		  in=1;
		}
	      else
		{
		  if (thispattern[k]<expression[sweeper][expression[sweeper][k]])
		    {
		      identical+=(featgrvals[k]*
				  (numdiff(expression[sweeper][expression[sweeper][k]],
					   thispattern[k])));
		      in=1;
		    }
		  else
		    {
		      identical+=(featgrvals[k]*
				  (numdiff(expression[sweeper][expression[sweeper][k]+1],
					   thispattern[k])));
		      in=1;
		    }
		}
	    }
	}
    }
  thisoccurrence=in*expressionoccurrence[sweeper];
}

/* select_applicable_bonuses - when feature value combinations are used,
   check for the currently applicable ones */
void select_applicable_bonuses(void)
{
  register int  i,j,k;
  char match;

  nrapplicables=0;

  for (i=nrbonuses-1; i>=0; i--)
  {
    match=1;
    k=0;
    for (j=0; ((j<PATWIDTH)&&(match)); j++)
      if (bonusselections[i][j]!=-1)
	{
	  if (bonusselections[i][j]!=thispattern[j])
	    match=0;
	  else
	    k++;
	}
    if ((match)&&(k>0))
      {
	applicable_bonuses[nrapplicables]=i;
	if (VERB4)
	  {
	    fprintf(stderr,"     applicable bonus %3d:",
		   nrapplicables);
	    for (j=0; j<PATWIDTH; j++)
	      {
		if (bonusselections[i][j]==-1)
		  fprintf(stderr," *");
		else
		  fprintf(stderr," %s",
			 values[j][bonusselections[i][j]]);
	      }
	    fprintf(stderr," - weight %f\n",
		   bonuses[i]);

	  }
	nrapplicables++;
      }
  }
  if (VERB4)
    fprintf(stderr,"     %d bonuses applicable\n",
	   nrapplicables);

  sort_bonuses();

  if (VERB4)
    {
      for (i=0; i<nrapplicables; i++)
	fprintf(stderr,"     sorted bonus %2d: %2d %f - restval %f\n",
	       i,orderedbonus[i],bonuses[applicable_bonuses[orderedbonus[i]]],
	       bonusrestvals[orderedbonus[i]]);
      fprintf(stderr,"     maxgrsum %f\n",
	     maxgrsum);
    }
}

/* select_applicable_bonuses_family - when feature value combinations are used,
   check for the currently applicable ones (family version) */
void select_applicable_bonuses_family(int inst)
{
  register int  i,j,k;
  char match;

  nrapplicables=0;

  for (i=nrbonuses-1; i>=0; i--)
    {
      match=1;
      k=0;
      for (j=0; ((j<PATWIDTH)&&(match)); j++)
	{
	  if (bonusselections[i][j]!=-1)
	    {
	      if (bonusselections[i][j]!=feature[j][inst])
		match=0;
	      else
		k++;
	    }
	}
      if ((match)&&(k>0))
	{
	  applicable_bonuses[nrapplicables]=i;
	  if (VERB4)
	    {
	      fprintf(stderr,"     applicable bonus %3d:",
		     nrapplicables);
	      for (j=0; j<PATWIDTH; j++)
		{
		  if (bonusselections[i][j]==-1)
		    fprintf(stderr," *");
		  else
		    fprintf(stderr," %s",
			   values[j][bonusselections[i][j]]);
		}
	      fprintf(stderr,"\n");
	    }
	  nrapplicables++;
	}
    }
  if (VERB4)
    fprintf(stderr,"     %d bonuses applicable\n",
	   nrapplicables);

  sort_bonuses();
}

/* sort_bonuses - figure out the ordering of current bonuses according
   to their weight, and precompute the summed weight to be gained when
   they are checked in that order */
void sort_bonuses(void)
{
  register int    i,j,hulp;
  char   ordered;

  /* bubble sort the bonus weights */
  for (i=0; i<nrapplicables; i++)
    orderedbonus[i]=i;
  ordered=0;
  while (!ordered)
    {
      for (i=0; i<nrapplicables-1; i++)
	if (bonuses[applicable_bonuses[orderedbonus[i]]]<
	    bonuses[applicable_bonuses[orderedbonus[i+1]]])
	  {
	    hulp=orderedbonus[i];
	    orderedbonus[i]=orderedbonus[i+1];
	    orderedbonus[i+1]=hulp;
	  }
      ordered=1;
      for (i=0; ((ordered)&&(i<nrapplicables-1)); i++)
	if (bonuses[applicable_bonuses[orderedbonus[i]]]<
	    bonuses[applicable_bonuses[orderedbonus[i+1]]])
	  ordered=0;
    }

  /* then compute rest values for optimised k-nn later on */
  for (i=0; i<nrapplicables; i++)
    {
      bonusrestvals[orderedbonus[i]]=0.0;
      for (j=i; j<nrapplicables; j++)
	{
	  if (METRIC==0)
	    bonusrestvals[orderedbonus[i]]+=bonuses[applicable_bonuses[orderedbonus[j]]];
	  else
	    bonusrestvals[orderedbonus[i]]+=(2.*bonuses[applicable_bonuses[orderedbonus[j]]]);
	}
    }
  bonusrestvals[nrapplicables]=0.0;
  maxgrsum=bonusrestvals[orderedbonus[0]];

}
