



/**************************************************************************************

                        ASSESSING DIRECTIONAL SPATIAL RELATIONSHIPS
				 BY IMPOSING PHYSICAL CONSIDERATIONS ON THE F-HISTOGRAMS

***************************************************************************************

									   Author
									   ------
		             P. Matsakis (pmatsakis@cecs.missouri.edu).

								       Version
								       -------
								 1.1, February 2001.

									   Based on
									   --------
     [1] P. Matsakis, Relations spatiales structurelles et interpretation d'images,
		 Ph.D. dissertation, IRIT, Universite Paul Sabatier, Toulouse, France, 1998.
	 [2] P. Matsakis, J. Keller, L. Wendling, J. Marjamaa, O. Sjahputera,
	     "Linguistic Description of Relative Positions in Images",
		 TSMC-B, vol. 31, no. 4, August 2001.

---------------------------------------------------------------------------------------

This module provides the functions degreeOfTruth() and degreeOfTruthSmu() described
below. An example of main() function is also presented at the very end of the file.






=======================================================================================
double degreeOfTruth (double *histogram, int numberDirections, double alpha);

degreeOfTruth(histogram,numberDirections,alpha) is equivalent to:
degreeOfTruthSmu(histogram,numberDirections,alpha,&STrapezoid,&muTriangular). 






=======================================================================================
double degreeOfTruthSmu (double *histogram, int numberDirections, double alpha,
					     double (*S)(double *,int), double (*mu)(double));

Assesses a proposition such as "A is in direction alpha of B".
---------------------------------------------------------------------------------------
in     | - histogram: A force histogram. The size of the array 'histogram'
       |   is 'numberDirections+1': 'histogram[0]' is the resultant of
	   |   forces in direction 0, 'histogram[numberDirections/4]' the
	   |   resultant of forces in direction PI/2, etc.
       | - numberDirections: The number of directions in which
       |   forces have been considered. It is supposed to be a
       |   positive multiple of 4 (e.g., 16, 32, 64, 120, 180).
       | - alpha: An angle measure that belongs to [0,2*PI]. 
       | - S: The aim of the function pointed by 'S' is to compute the
       |   threshold 'tau' (see [2]). The present module provides three
	   |   such functions: 'SMin', 'SMax', and 'STrapezoid'.
       | - mu: The function pointed by 'mu' defines the directional
	   |   spatial relations between points (see [2]). The present
	   |   module provides one such function: 'muTriangular'.
---------------------------------------------------------------------------------------
return | Assume that 'histogram' is a force histogram associated with two
       | objects A and B (see "FHistogramRaster.c" and "FHistogramVector.c").
	   | The value returned by 'degreeOfTruthSmu' is the degree of truth of the
	   | proposition "A is in direction alpha of B" (in the histogram's opinion).
	   | It belongs to the interval [0,1]. It is 0 if the proposition is totally
	   | false, 1 if it is totally true.


**************************************************************************************/








#include <stdlib.h>
#include <math.h>
#include <stdio.h>




/*==========================================================================
Constantes symboliques
----------------------------------------------------------------------------
Les deux premieres constantes sont utilisees par la fonction
'STrapezoid', les deux suivantes par la fonction 'muTriangular'.

Le seuil force est calcule a partir d'un trapeze symetrique par rapport a
l'angle etudie. SEUIL_PARAM_HT est compris entre 0.0 et 1.0 : il indique
la demi-longueur du noyau du trapeze, l'unite etant pi/2. SEUIL_PARAM_BS
est compris entre SEUIL_PARAM_HT et 1.0 : il indique la demi-longueur du
support.

La relation est representee par un trapeze lui-aussi symetrique par rapport
a l'angle etudie. RELAT_PARAM_HT et RELAT_PARAM_BS definissent ce trapeze.
==========================================================================*/


#define SEUIL_PARAM_HT	0.25
#define	SEUIL_PARAM_BS	0.5

#define RELAT_PARAM_HT	0.0
#define	RELAT_PARAM_BS	1.0

#define PI 3.141592653589 




/*==========================================================================
These 'S' functions may be used by 'degreeOfTruthSmu'.
----------------------------------------------------------------------------
in      | - histogram: A histogram of (effective) forces.
        |   The size of 'histogram' is 'numberDirections+1'.
	    |   'histogram[numberDirections/2]' is the resultant
	    |   force in the direction we are interested in.
        | - numberDirections: The number of directions in which 
        |   forces have been considered. It is supposed to be a
	    |   positive multiple of 4 (e.g., 16, 32, 64, 120, 180).
----------------------------------------------------------------------------
return  | The threshold 'tau' (see [2]).
----------------------------------------------------------------------------
WARNING | If you want to define your own 'S' function, keep in mind the
        | fact that for any 'i' such that 'i<=numberDirections/4' or
	    | 'i>=3*numberDirections/4', your 'S' function should not use
		| 'histogram[i]' (which is supposed to be zero anyway).
==========================================================================*/


double SMin (double *histogram, int numberDirections) {

	int i, n;
	double min;

	n=numberDirections/4;
	i=3*n-1; min=histogram[i];
	for(--i;i>n;--i) if(histogram[i]<min) min=histogram[i];
	return(min);
}


double SMax (double *histogram, int numberDirections) {

	int i, n;
	double max;

	n=numberDirections/4;
	i=3*n-1; max=histogram[i];
	for(--i;i>n;--i) if(histogram[i]>max) max=histogram[i];
	return(max);
}


double STrapezoid (double *histogram, int numberDirections) {

	int i, ht, bs, n, auxi;
	double increment, poids, somme;

	n=numberDirections/4;
	ht=n*SEUIL_PARAM_HT;
	bs=n*SEUIL_PARAM_BS;
	if(bs==ht) bs++;
	somme=0.0;
	n*=2;

	/* Partie laterale gauche.
	--------------------------*/
	auxi=n-ht;
	increment=1.0/(bs-ht);
	poids=increment;
	for(i=n-bs;++i<auxi;poids+=increment) somme+=histogram[i]*poids;

	/* Partie centrale.
	-------------------*/
	/* i=numberDirections/2-ht; */
	for(auxi=n+ht;i<=auxi;i++) somme+=histogram[i];

	/* Partie laterale droite.
	--------------------------*/
	/* i=numberDirections/2+ht+1; */
	for(auxi=n+bs;i<auxi;i++)
	  {poids-=increment; somme+=histogram[i]*poids;}

	return(somme/(bs+ht));
}




/*==========================================================================
muTriangular | This function may be used by 'degreeOfTruthSmu'. It is from
               [-1,1] (in fact, -1 corresponds to -PI/2 and 1 corresponds
			   to PI/2) to [0,1]. It takes the value 0 at -1 and 1, and
			   the value 1 at 0. It is even, continuous, and decreasing
			   on [0,1] (see [2]).
----------------------------------------------------------------------------
in     | A value 'f' such that: -1.0<=f<=1.0.
----------------------------------------------------------------------------
return | A value 'v' such that: 0.0<=v<=1.0.
==========================================================================*/


double muTriangular (double f) {

	f=fabs(f);
	if(f>=RELAT_PARAM_BS) return(0.0);
	if(f<=RELAT_PARAM_HT) return(1.0);
	return((RELAT_PARAM_BS-f)/(RELAT_PARAM_BS-RELAT_PARAM_HT));
}




/*============================================================================
degreeOfTruthSmu | This function is described at the beginning of the file.
============================================================================*/


double degreeOfTruthSmu (double *histogram,
						 int numberDirections, double alpha,
					     double (*S)(double *,int), double (*mu)(double)) {

	int i, j, auxi, n;
	double dist, moment, forces_pos, forces_neg, auxf;
	double *tab, forces_ccos[6], theta_mp0[3];


	n=numberDirections/4;
	auxi=(int)(0.5+0.5*numberDirections*alpha/PI);
	auxi=(numberDirections/2+auxi)%numberDirections;
	tab=(double *)malloc((numberDirections+1)*sizeof(double));
	for(i=auxi,j=0;i<numberDirections;i++,j++) tab[j]=histogram[i];
	for(i=0;i<=auxi;i++,j++) tab[j]=histogram[i];


	/*-----------------------------------------------------------
	   Partie [-pi;-pi/2] : forces et moment par rapport a -pi/2.
	   forces=0.5tab[0]+tab[1]+...+tab[n]
	   moment=0.5tab[0](n-0.25)+tab[1](n-1)+...+tab[n-1](n-(n-1))
	-------------------------------------------------------------*/
	forces_neg=0.5*tab[0];
	dist=(double)n;
	moment=forces_neg*(dist-0.25);
	for(i=1;i<n;i++) {moment+=tab[i]*--dist; forces_neg+=tab[i];}
	forces_neg+=tab[i];
	forces_ccos[0]=forces_neg;

	/* Le moment precedent engendre une
	   alteration de la partie [-pi/2;pi/2]. 
	-----------------------------------------*/
	auxi=3*n;
	/* i=n; dist=1.0; */
	auxf=tab[++i]*dist;
	while(moment>auxf && i<auxi)
	  {forces_neg+=tab[i]; tab[i]=0.0;
	   moment-=auxf; auxf=tab[++i]*++dist;}
	if(i<auxi) {auxf=moment/dist; forces_neg+=auxf; tab[i]-=auxf;}
	theta_mp0[0]=(-2.0+i/(double)n)*(PI/2.0);
	forces_ccos[2]=forces_neg-forces_ccos[0];

	/*-------------------------------------------------------
	   Partie [pi/2;pi] : forces et moment par rapport a pi/2.
	---------------------------------------------------------*/
	/* auxi=3*n; */
	auxf=0.5*tab[0];
	forces_neg+=auxf;
	dist=(double)n;
	moment=auxf*(dist-0.25);
	for(i=(n<<2)-1;i>auxi;i--)
	  {moment+=tab[i]*--dist; forces_neg+=tab[i];}
	forces_neg+=tab[i];
	forces_ccos[1]=forces_neg-(forces_ccos[0]+forces_ccos[2]);
	
	/* Le moment precedent engendre une
	   alteration de la partie [-pi/2;pi/2].
	-----------------------------------------*/
	auxi=n;
	/* i=3*n; dist=1.0; */
	auxf=tab[--i]*dist;
	while(moment>auxf && i>auxi)
	  {forces_neg+=tab[i]; tab[i]=0.0;
	   moment-=auxf; auxf=tab[--i]*++dist;}
	if(i>auxi) {auxf=moment/dist; forces_neg+=auxf; tab[i]-=auxf;}
	theta_mp0[1]=(-2.0+i/(double)n)*(PI/2.0);
	forces_ccos[3]=forces_neg-
	               (forces_ccos[0]+forces_ccos[1]+forces_ccos[2]);

	/*------------------------------------------
	   Calcul du seuil et des forces favorables.
	   Traitement du cas le plus defavorable.
	--------------------------------------------*/
	auxi=3*n;
	forces_pos=0.0; forces_ccos[4]=forces_ccos[5]=0.0;
	auxf=(*S)(tab,numberDirections);
	for(i=n+1;i<auxi;i++)
	  {
	  forces_pos+=tab[i];
	  if(tab[i]>auxf)
	    {tab[i]-=auxf; forces_ccos[4]+=tab[i]; forces_ccos[5]+=auxf;}
	  else
	    {forces_ccos[5]+=tab[i]; tab[i]=0.0;}
	  }
	if(!forces_pos) {theta_mp0[2]=0.0; return(0.0);}

	/* Calcul du degre dans le cas general.
	---------------------------------------*/
	auxi=n<<1;
	moment=0.0;
	dist=(double)n;
	/* i=3*n; */
	while(--i>auxi) moment+=tab[i]*--dist;
	while(--i>n) moment-=tab[i]*dist++;
	theta_mp0[2]=(PI/2.0)*(moment/(forces_pos*(double)n));
	auxf=(*mu)(moment/(forces_pos*(double)n));
	return(auxf*(forces_pos/(forces_pos+forces_neg)));
}




/*============================================================================
degreeOfTruth | This function is described at the beginning of the file.
============================================================================*/


double degreeOfTruth (double *histogram, int numberDirections, double alpha) {

	double f;

	f=degreeOfTruthSmu(histogram,numberDirections,alpha,
		               &STrapezoid,&muTriangular);
	return(f);
}








/*******************************************************************************
						    EXAMPLE OF main() FUNCTION
*******************************************************************************/






/*==============================================================================
   main | This function is provided as an example and for testing purposes.
--------------------------------------------------------------------------------
   Compile (e.g., "gcc FHistogramDegree.c -lm"), execute ("a.out"),
   and follow the instructions.
==============================================================================*/


void main () {


  	FILE *histogramFile;
	char nameHistogramFile[100];
	int i, alpha, numberOfDirections;
	double *histogram, d;


	/* Choosing the arguments.
	--------------------------*/

	printf("\n\nHistogram.\n---------\n");
	printf("Enter the name of the file it is stored in: ");
	scanf("%s",nameHistogramFile);
	printf("\nDirection.\n---------\n");
	printf("Enter an integer value (angle measure in degrees, between 0 and 360): ");
	scanf("%d",&alpha);
	

	/* Reading the histogram file.
	------------------------------*/

	if((histogramFile=fopen(nameHistogramFile,"rt"))==NULL)
		{printf("\nERROR histogram file\n\n"); exit(1);}
	fscanf(histogramFile,"%d ",&numberOfDirections);
	histogram=(double *)malloc((numberOfDirections+1)*sizeof(double));
	for(i=0;i<=numberOfDirections;i++)
		fscanf(histogramFile,"%lf ",&histogram[i]);
	fclose(histogramFile);


	/* Assessing "A is in direction alpha of B".
	--------------------------------------------*/

	printf("\nDegree of truth of \"A is in direction %d of B\" is:\n",alpha);
	printf("---------------------------------------------------\n");
    d=degreeOfTruthSmu(histogram,numberOfDirections,alpha*PI/180,
		               &SMin,&muTriangular);
	printf("%.2f when choosing SMin\n",d);
    d=degreeOfTruth(histogram,numberOfDirections,alpha*PI/180);
	printf("%.2f when choosing STrapezoid (default choice)\n",d);
    d=degreeOfTruthSmu(histogram,numberOfDirections,alpha*PI/180,
		               &SMax,&muTriangular);
	printf("%.2f when choosing SMax\n\n",d);


	/* Done.
	--------*/

	free(histogram);
}



