



/**************************************************************************************

                          F-HISTOGRAM COMPUTATION (VECTOR DATA)

***************************************************************************************

									   Author
									   ------
		              P. Matsakis (pmatsakis@cecs.missouri.edu).

								      Version
								      -------
								2.1, February 2001.
				     For details (in French), see Section III.

									  Based on
									  --------
     [1] P. Matsakis, Relations spatiales structurelles et interpretation d'images,
		 Ph.D. dissertation, IRIT, Universite Paul Sabatier, Toulouse, France, 1998.
	 [2] P. Matsakis, L. Wendling, "A New Way to Represent the Relative Position of
	     Areal Objects", PAMI, vol. 21, no. 7, pp. 634-643, 1999.

---------------------------------------------------------------------------------------

This module provides the function FRHistogram_CrispVector(). The objects that can be
handled by this function are described in Section I. The function itself is described
in Section II. An example of main() function is presented at the very end of the file.






=======================================================================================
I. POLYGONAL OBJECTS and FILE FORMAT
---------------------------------------------------------------------------------------

This module handles polygonal objects: each object is represented by a set of closed
polylines; each polyline is represented by a set of vertices; each vertex is repres-
ented by its coordinates x and y. The file that defines an object is an ASCII (text)
file. It contains different values, separated by space or line feed characters (' '
or '\n'). 

    1st value:   The total number of vertices (integer).
    2nd value:   The number of vertices of the 1st closed polyline (integer).
    Next values: The coordinates x and y of the 1st vertex of the 1st polyline, 
                 the coordinates x and y of the 2nd vertex of the 1st polyline,
                 etc. These values are floating point numbers. For a square, 8
                 coordinates (x1 y1 x2 y2 x3 y3 x4 y4) should be specified.
    Value i:     The number of vertices of the 2nd polyline (integer).
    Next values: The coordinates x and y of the 1st vertex of the 2nd polyline, 
                 the coordinates x and y of the 2nd vertex of the 2nd polyline,
                 etc.
    And so forth. 

The vertices of a polyline can be listed clockwise or anti-clockwise (it doesn't
matter). Non-convex objects, non-connected objects and objects with holes can be
considered. The only constraints are that: (i) a polyline should not "intersect"
itself (the same vertex cannot appear twice, an edge cannot intersect another
edge), and (ii) two distinct polylines should not intersect either. Note that
here, the y-coordinates go from bottom to top (math convention), and not from
top to bottom.






=======================================================================================
II.

void FRHistogram_CrispVector (double *histogram,
                              int numberDirections, double typeForce,
                              char *objectA, char *objectB);

Computation of the Fr-histogram associated with a pair of
polygonal objects defined by their vertices. The coordinates
of these vertices are stored in text files.
---------------------------------------------------------------------------------------
in  | - typeForce: Can be any real number (0.0 for the computation of
    |   the histogram of constant forces, 2.0 for the computation of
    |   the histogram of gravitational forces, etc.).
    | - numberDirections: It is a positive multiple of 4 (e.g., 16, 32,
	|   64, 120, 180). Forces will be considered in 'numberDirections'
	|   directions.
    | - objectA: The name of the file that defines the argument object.
    | - objectB: The name of the file that defines the referent object.
---------------------------------------------------------------------------------------
out | The histogram of forces. 'numberDirections+1' values will be stored in
    | the array 'histogram'. 'histogram[0]' will be the resultant of forces
	| in direction 0 ("to the right"), 'histogram[numberDirections/2]' the
    | the resultant of forces in direction PI ("to the left"), and the last
    | value, 'histogram[numberDirections]', will be the resultant in direct-
    | ion 2PI (equal to the first value 'histogram[0]').
---------------------------------------------------------------------------------------
   WARNING | - The arguments are not checked.
           | - The memory space pointed by 'histogram'
           |   must be allocated by the calling function.
		   | - Currently, only disjoint objects can be handled.
---------------------------------------------------------------------------------------
      NOTE | In theory:  
		   | - If the interior of A intersects the interior of B, 
		   |   (A,B) is assessable iff 'typeForce<1'.
		   | - If A is adjacent to B (A and B intersect, but the
		   |   interior of A does not intersect the interior of B),
		   |   (A,B) is assessable iff 'typeForce<2'.
		   | - Otherwise (disjoint objects),
		   |   (A,B) is always assessable.






=======================================================================================
III. VERSIONS
---------------------------------------------------------------------------------------

La version 1.0 est "relvect.c" d'avril 1996. La seule difference significative
entre 1.0 et 1.1 concerne la fonction qui s'appelait "Poly" et qui s'appelle
maintenant "FRHistogram_CrispVector".La version 1.2 corrige un bug concernant
le facteur multiplicatif qui rentre en jeu lors du calcul de la resultante des
forces dans la direction theta. La version 2.0 permet le calcul d'histogrammes
quelconques  (seuls les histogrammes F0 et F2 pouvaient etre calcules dans les
versions precedentes) . De plus,  la memoire est proprement geree  (ce n'etait
pas le cas avant :  pour des raisons obscures, les appels a "free" avaient ete
supprimes). Differents couples d'objets peuvent donc etre traites les uns a la
suite des autres sans que le programme ne plante (du moins, je l'espere). Les
fonctions "HistoStretch" et "main" (en fin de fichier) avaient ete rajoutees
dans la version 1.2 afin de tester le module : elles ont ete conservees. Pas
de differences significatives entre 2.0 et 2.1. L'organisation de ce "readme"
integre ainsi que quelques noms ont ete changes ; la premiere valeur de l'his-
togramme n'est plus la resultante des forces dans la direction -PI, mais la
resultante des forces dans la direction 0. Il s'agissait d'uniformiser ce mo-
dule avec le module de calcul des histogrammes associes a des donnees de type
raster (cf. "FHistogramRaster.c"). Les deux modules constituent maintenant
avec "FHistogramDegree.c" un package "FHistogram" plus presentable.

**************************************************************************************/








/* include
************/


#include <stdlib.h>
#include <math.h>
#include <stdio.h>




/* constantes symboliques
***************************/


/*===========================================================================
SEUIL_ZERO est utilise pour tester la nullite de certains parametres (nul-
lite qui conduit a l'instanciation de formules particulieres). CTE_F1 et
CTE_F2 servent au calcul des forces en 1/d et 1/(d*d) (cf 'force_trapezes').
Pour debugger, definir la constante DEBUG.
===========================================================================*/

/* #define         DEBUG           1 */
#define 	PI  		3.1415926535897932384626433832795
#define 	SEUIL_ZERO 	0.001

#define 	CTE_F1 		0.69314718055994530941723212145818 /* ln(2.0) */
#define 	CTE_F2 		0.30685281944005469058276787854182 /* 1+ln(0.5) */




/* types symboliques de donnees
*********************************/


/*===========================================================================
cell_objet | Cellule objet.
-----------------------------------------------------------------------------
Un objet est represente par un ensemble de suites de points (les sommets).
Chaque suite definit une polyligne simple et fermee. L'union des polylignes
constitue la frontiere de l'objet. A tout sommet correspond une cellule ob-
jet et reciproquement. A tout objet correspond un tableau objet : tableau
de cellules objets. Deux objets sont consideres simultanement : l'objet A
et l'objet B. Une relation d'ordre est definie sur l'ensemble des sommets
de ces objets.
-----------------------------------------------------------------------------
x    | Abscisse du sommet.
y    | Ordonnee du sommet.
pred | Indice (dans le tableau objet) du predecesseur (dans la polyligne).
succ | Indice (dans le tableau objet) du successeur (dans la polyligne).
u    | Abscisse du sommet dans le repere utilise pour l'ordination des som-
     | mets : cette abscisse constitue le premier critere de tri.
v    | Ordonnee du sommet dans le repere utilise pour l'ordination des som-
     | mets : cette ordonnee constitue le second critere de tri.
pdt  | Pointeur sur la cellule objet rattachee au sommet precedant (dans la
     | liste triee des sommets).
svt  | Pointeur sur la cellule objet rattachee au sommet suivant (dans la
     | liste triee des sommets).
===========================================================================*/


struct cell_objet
  {
  double x, y;
  int pred, succ;

  double u, v;
  struct cell_objet *pdt, *svt;
  };




/*===========================================================================
cell_trap | Cellule trapeze.
-----------------------------------------------------------------------------
Des droites paralleles distinctes, ordonnees selon u, vont sectionner les
objets en trapezes. Sur une droite de section donnee, les sommets des tra-
pezes d'un objet donne seront ordonnes selon v : chacun de ces sommets se-
ra represente par une cellule trapeze, la suite de ces sommets sera repre-
sentee par une liste trapeze (liste simplement chainee de cellules trapezes).
-----------------------------------------------------------------------------
svt        | Pointeur sur la cellule trapeze suivante
           | (dans la liste trapeze).
v          | Ordonnee du sommet.
est_sommet | 1 si le sommet du trapeze correspond a un sommet de l'objet,
           | 0 sinon.
point      | Si 'est_sommet' vaut 1, indice dans le tableau objet du sommet
           | en question. Sinon, le sommet du trapeze est un point non sommet
           | de l'objet. La polyligne dont il fait partie le relie par contre
           | a deux sommets de l'objet : l'un a deja ete traite (sections an-
           | terieures) l'autre non. 'point' est l'indice dans le tableau ob-
           | jet du sommet non encore traite.
mult       | 0 : faux sommet (a ignorer),
           | 1 : sommet simple (trapeze non degenere),
           | 2 : sommet double (trapeze degenere : triangle).
a          | Si 'est_sommet' vaut 0, a et b permettent de calculer v en
b          | fonction de u : v=au+b (suivi des donnees d'une section a
           | l'autre dans le but de reduire les calculs).
===========================================================================*/


struct cell_trap
  {
  struct cell_trap *svt;
  double v;
  int est_sommet;
  int point;
  int mult;
  double a, b;
  };




/*===========================================================================
cell_yi | Cellule yi (cf la fonction 'init_yi').
-----------------------------------------------------------------------------
A un instant donne deux sections consecutives seront considerees : l'une ap-
pelee 1, l'autre appelee 2. Ces sections definissent un ensemble de trapezes
de A, un ensemble de trapezes de B. Sur chaque section on retrouve un ensem-
ble de 'nb_trapA' bases de trapezes de A, un ensemble de 'nb_trapB' bases de
trapezes de B. Les longueurs des differentes bases seront stockees dans dif-
ferents tableaux de flottants ('xiA1','xiA2','xiB1','xiB2'), tandis que les
positions des bases des trapezes de A relativement aux bases des trapezes de
B le seront dans differents tableaux de cellules yi ('yiA1B1', 'yiA2B2').
===========================================================================*/


struct cell_yi
  {
  double a1, a2;
  double b1, b2;
  double c1, c2;
  double d1, d2;
  };




/* variables a portee locale
******************************/


/*===========================================================================
Les forces sont en 1/d^type_force. Par exemple, elles sont constantes lors-
que 'type_force' vaut 0 et gravitationnelles lorsque 'type_force' vaut 2.0. 
'objetA' et 'objetB' sont les tableaux objets : 'objetA_max-objetA' et
'objetB_max-objetB' sont les nombres d'elements de ces tableaux. 'pc_tri'
va balayer la liste triee des sommets des objets. 'listeA1', 'listeA2',
'listeB1' et 'listeB2' sont les listes trapezes.
===========================================================================*/


static double type_force=0.0;

static struct cell_objet *objetA=NULL, *objetA_max;
static struct cell_objet *objetB=NULL, *objetB_max;
static struct cell_objet *pc_tri;

static struct cell_trap *listeA1=NULL, *listeB1=NULL;
static struct cell_trap *listeA2=NULL, *listeB2=NULL;

static int nb_trapA=0, nb_trapB=0;

static double *xiA1=NULL, *xiB1=NULL;
static double *xiA2=NULL, *xiB2=NULL;

static struct cell_yi *yiA1B1=NULL;
static struct cell_yi *yiA2B2=NULL;




/* declaration des fonctions a portee locale
**********************************************/


static double trier_sommets(double theta, double *umin, double *umax);

static void init_trapezes(double ucrt);
static double maj_trapezes();
static void inserer_sommets(double ucrt, int *p_multA, int *p_multB);

static void detruire(struct cell_trap *liste);
static struct cell_trap *copie(struct cell_trap *liste1);

static double force_trapezes();

static void init_xi(char nom_objet, int num_liste);
static void init_yi(int num_liste);

static void debug_xi(char nom_objet, int num_liste);
static void debug_yi(int num_liste);
static void debug_cell_trap (struct cell_trap *liste);
static void debug_pc_tri();








/****************************************************************************
								LOCAL FUNCTIONS
****************************************************************************/






/*===========================================================================
trier_sommets | Constitution de la liste de tri : mise a jour des
		champs 'u', 'v', 'svt', 'pdt' de chaque cellule objet.
-----------------------------------------------------------------------------
in    | Un element de ]-PI;PI]. Le tri est effectue d'abord selon une
	direction perpendiculaire a 'theta', ensuite selon 'theta'.
inout | Les deux objets ne coexistent dans la direction 'theta' que pour
	des valeurs de u appartenant a l'intervalle [*umin;*umax] (pas de
	coexistence si : *umin>*umax).
out   | La constante 'e_theta', egale a + ou - cos ou sin de 'theta'.
-----------------------------------------------------------------------------
En sortie, le pointeur 'pc_tri' de la liste de
tri pointe sur la premiere cellule de la liste.
===========================================================================*/


static double trier_sommets(double theta, double *umin, double *umax)
  {
  int i, param_tri2;
  double f, g, param_tri1, e_theta, uAmin, uAmax, uBmin, uBmax;
  struct cell_objet *sommet_crt, *sommet_fin, *pc_tri_min;

  /* Initialisations, notamment des parametres de tri.
     Un sommet (x,y) sera range dans la liste de tri d'une part en fonction
     de l'une des valeurs y-x.tan(theta) et y.cotan(theta)-x, d'autre part
     en fonction de la position de son projete sur l'un des 4 axes orientes.
  -------------------------------------------------------------------------*/

  pc_tri=NULL;
  if((f=PI*0.75)<=theta || theta<-f)
    {
    e_theta=-cos(theta); param_tri1=tan(theta); param_tri2=2;
    uAmin=uAmax=objetA->y-param_tri1*objetA->x;
    uBmin=uBmax=objetB->y-param_tri1*objetB->x;
    }
  else if(theta<-(f=PI*0.25))
    {
    e_theta=-sin(theta); param_tri1=tan(PI/2.0-theta); param_tri2=-1;
    uAmin=uAmax=param_tri1*objetA->y-objetA->x;
    uBmin=uBmax=param_tri1*objetB->y-objetB->x;
    }
  else if(theta<f)
    {
    e_theta=cos(theta); param_tri1=tan(theta); param_tri2=0;
    uAmin=uAmax=objetA->y-param_tri1*objetA->x;
    uBmin=uBmax=objetB->y-param_tri1*objetB->x;
    }
  else
    {
    e_theta=sin(theta); param_tri1=tan(PI/2.0-theta); param_tri2=1;
    uAmin=uAmax=param_tri1*objetA->y-objetA->x;
    uBmin=uBmax=param_tri1*objetB->y-objetB->x;
    }

  /* Il va falloir considerer chaque objet.
  -----------------------------------------*/
  for(sommet_crt=objetA,sommet_fin=objetA_max,i=0;
      i<2;
      sommet_crt=objetB,sommet_fin=objetB_max,i++)

    /* Il va falloir considerer chacun des sommets de l'objet.
    ----------------------------------------------------------*/
    for(;sommet_crt<sommet_fin;sommet_crt++)
      {
      /* Mise a jour des champs 'u' et 'v'.
      -------------------------------------*/
      switch(param_tri2)
	{
	case -1 :
	  sommet_crt->u=param_tri1*sommet_crt->y-sommet_crt->x;
	  sommet_crt->v=-sommet_crt->y;
	  break;
	case  0 :
	  sommet_crt->u=sommet_crt->y-param_tri1*sommet_crt->x;
	  sommet_crt->v=sommet_crt->x;
	  break;
	case  1 :
	  sommet_crt->u=param_tri1*sommet_crt->y-sommet_crt->x;
	  sommet_crt->v=sommet_crt->y;
	  break;
	case  2 :
	  sommet_crt->u=sommet_crt->y-param_tri1*sommet_crt->x;
	  sommet_crt->v=-sommet_crt->x;
	  break;
	}

      /* Mise a jour des valeurs 'uAmin', 'uAmax', 'uBmin', 'uBmax'.
      --------------------------------------------------------------*/
      f=sommet_crt->u;
      g=sommet_crt->v;
      if(i) {if(f<uBmin) uBmin=f; else if(f>uBmax) uBmax=f;}
      else {if(f<uAmin) uAmin=f; else if(f>uAmax) uAmax=f;}

      /* Mise a jour des champs 'pdt' et 'svt' (tri).
      -----------------------------------------------*/
      if(!pc_tri)
	{
	/* C'est la premiere mise a jour effectuee... */
	sommet_crt->svt=sommet_crt->pdt=NULL;
	pc_tri=pc_tri_min=sommet_crt;
	}
      else
	{
	/* Ce n'est pas la premiere mise a jour effectuee... */
	while((f>pc_tri->u || (f==pc_tri->u && g>pc_tri->v)) && pc_tri->svt)
	  pc_tri=pc_tri->svt;
	while((f<pc_tri->u || (f==pc_tri->u && g<pc_tri->v)) && pc_tri->pdt)
	  pc_tri=pc_tri->pdt;
	if(f>pc_tri->u || (f==pc_tri->u && g>pc_tri->v))
	  {
	  sommet_crt->pdt=pc_tri;
	  sommet_crt->svt=pc_tri->svt;
	  pc_tri->svt=sommet_crt;
	  if(sommet_crt->svt) (sommet_crt->svt)->pdt=sommet_crt;
	  }
	else
	  {
	  sommet_crt->pdt=pc_tri->pdt;
	  sommet_crt->svt=pc_tri;
	  pc_tri->pdt=sommet_crt;
	  if(sommet_crt->pdt) (sommet_crt->pdt)->svt=sommet_crt;
	  else pc_tri_min=sommet_crt;
	  }
	}
      }

  *umin=uAmin<uBmin?uBmin:uAmin;
  *umax=uAmax<uBmax?uAmax:uBmax;
  pc_tri=pc_tri_min;
  return(e_theta);
  }




/*===========================================================================
init_trapezes | Constitue les listes 2 associees a un rang donne, met a jour
		'nb_trapA' et 'nb_trapB'.
-----------------------------------------------------------------------------
in | Le rang de l'un des points de l'un des objets.
-----------------------------------------------------------------------------
Lorsque 'init_trapeze' est appelee, 'pc_tri' est suppose pointer sur la pre-
miere cellule de la liste de tri (qui doit donc etre deja constituee), les
anciennes listes 2 sont detruites. Apres l'appel, 'pc_tri' pointe sur la pre-
miere cellule de tri dont la valeur du champ 'u' est strictement superieure
a 'ucrt' (ou est le pointeur NULL s'il n'existe pas de telle cellule).
===========================================================================*/


static void init_trapezes(double ucrt)
  {
  int i, l, *p_nb;
  struct cell_objet *objet;
  struct cell_trap *c_liste, *nc_liste, **p_liste;

  detruire(listeA2);
  detruire(listeB2);
  listeA2=listeB2=NULL;
  nb_trapA=nb_trapB=0;
  for(;pc_tri->u<ucrt;pc_tri=pc_tri->svt)
    {
    /* A quel objet est rattache le sommet traite ?
    -----------------------------------------------*/
    if(pc_tri>=objetA && pc_tri<objetA_max)
      {objet=objetA; p_liste=&listeA2; p_nb=&nb_trapA;}
    else
      {objet=objetB; p_liste=&listeB2; p_nb=&nb_trapB;}

    /* Ou sont situes les predecesseur et successeur de ce sommet ?
    ---------------------------------------------------------------*/
    for(i=pc_tri->pred,l=0;l<2;l++,i=pc_tri->succ)
      if(objet[i].u>ucrt)
	{
	/* Creation puis initialisation d'une cellule trapeze.
	------------------------------------------------------*/
	nc_liste=(struct cell_trap *)malloc(sizeof(struct cell_trap));
	nc_liste->point=i;
	nc_liste->est_sommet=0;
	nc_liste->mult=1;
	nc_liste->a=(objet[i].v-pc_tri->v)/(objet[i].u-pc_tri->u);
	nc_liste->b=pc_tri->v-nc_liste->a*pc_tri->u;
	nc_liste->v=nc_liste->a*ucrt+nc_liste->b;
	(*p_nb)++;

	/* Insertion de la cellule dans la liste 2 adequate.
	----------------------------------------------------*/
	if(!*p_liste || (*p_liste)->v>nc_liste->v)
	  {
	  nc_liste->svt=*p_liste;
	  *p_liste=nc_liste;
	  }
	else
	  {
	  for(c_liste=*p_liste;
	      c_liste->svt && c_liste->svt->v<nc_liste->v;
	      c_liste=c_liste->svt);
	  nc_liste->svt=c_liste->svt;
	  c_liste->svt=nc_liste;
	  }
	}
    }

  /* Il ne reste plus qu'a inserer les sommets de rang 'ucrt'.
  ------------------------------------------------------------*/
  inserer_sommets(ucrt,&nb_trapA,&nb_trapB);
  nb_trapA>>=1; nb_trapB>>=1;
  if(xiA2) free(xiA2);
  init_xi('A',2);
  if(xiB2) free(xiB2);
  init_xi('B',2);
  if(yiA2B2) free(yiA2B2);
  init_yi(2);
  }




/*===========================================================================
maj_trapezes | Constitue toutes les listes 2 associees au rang courant
	       ('listeA2', 'listeB2', 'xiA2', 'xiB2' et 'yiA2B2'), toutes
	       les listes 1 associees au rang qui precede le rang courant
	       ('listeA1', 'listeB1', 'xiA1', 'xiB1' et 'yiA1B1') et met a
	       jour 'nb_trapA' et 'nb_trapB'.
-----------------------------------------------------------------------------
out | Le rang u qui etait le rang courant au moment de l'appel.
-----------------------------------------------------------------------------
Lorsque 'maj_trapeze' est appelee, 'listeA2' et 'listeB2' sont supposees etre
les listes associees au rang qui precede le rang courant u, 'pc_tri' est sup-
pose pointer sur la premiere cellule de la liste de tri associee a un sommet
de rang u. Apres l'appel, 'pc_tri' pointe sur la premiere cellule de tri as-
sociee a un sommet de rang strictement superieur a u (est le pointeur NULL
si une telle cellule n'existe pas).
===========================================================================*/


static double maj_trapezes()
  {
  double u;
  int i, j, k, l, s, t;
  struct cell_objet *objet;
  struct cell_trap *c1_liste, *c2_liste, *c_liste, **p_liste;

  u=pc_tri->u;

  /*-------------------------------------
     Constitution des nouvelles listes 1.
  ---------------------------------------*/

  /* Apres modification des champs 'mult', les anciennes 'listeA2'
     et 'listeB2' vont devenir les nouvelles 'listeA1' et 'listeB1'.
  ------------------------------------------------------------------*/
  detruire(listeA1);
  detruire(listeB1);

  listeA1=listeA2;
  listeB1=listeB2;
  t=0;
  for(objet=objetA,c1_liste=listeA1,l=1;
      l<3;
      l++,c1_liste=listeB1,objet=objetB)
    {
    for(s=0;c1_liste;c1_liste=c1_liste->svt)
      if(c1_liste->est_sommet)
	{
	i=objet[c1_liste->point].succ;
	j=objet[c1_liste->point].pred;
	/* attention : u n'est pas le rang correspondant. */
	if(objet[i].u>=u && objet[j].u>=u) k=2;
	else if(objet[i].u<u && objet[j].u<u) k=0;
	else k=1;
	if(k!=c1_liste->mult) {t|=l; c1_liste->mult=k;}
	s+=k;
	}
      else s++;
    if(l==1) nb_trapA=s>>1;
    else nb_trapB=s>>1;
    }

  /* S'il y a effectivement eu des modifications, il
     faut les repercuter sur 'xiA1', 'xiB1' et 'yiA1B1'.
     Sinon, il suffit de reprendre les anciennes valeurs.
  -------------------------------------------------------*/

  if(xiA1) free(xiA1);
  if(xiB1) free(xiB1);
  if(yiA1B1) free(yiA1B1);
  if(t&1) {init_xi('A',1); if(xiA2) free(xiA2);} else xiA1=xiA2;
  if(t&2) {init_xi('B',1); if(xiB2) free(xiB2);} else xiB1=xiB2;
  if(t) {init_yi(1); if(yiA2B2) free(yiA2B2);} else yiA1B1=yiA2B2;

  /*----------------------------------
     Les nouvelles listes 2 vont elles
     aussi se deduire des anciennes...
  ------------------------------------*/

  listeA2=NULL;
  listeB2=NULL;

  /* On va parcourir chaque ancienne liste...
  -------------------------------------------*/
  for(l=0,c1_liste=listeA1,p_liste=&listeA2,objet=objetA;
      l<2;
      l++,c1_liste=listeB1,p_liste=&listeB2,objet=objetB)

    /* ... chaque cellule de chaque ancienne liste.
    -----------------------------------------------*/
    for(;c1_liste;c1_liste=c1_liste->svt)
      {
      /* Chacune peut donner naissance a au plus
	 deux cellules de la nouvelle liste 2...
      ------------------------------------------*/
      s=c1_liste->point;
      if(c1_liste->est_sommet) {i=objet[s].succ; j=objet[s].pred; t=2;}
      else {i=s; t=1;}
      for(k=i;t;t--,k=j)
	if(objet[k].u>u)
	  {
	  /* L'inegalite precedente est stricte, car les
	     sommets seront traites plus tard (en vrac).
	     Creation puis initialisation d'une cellule trapeze.
	  ------------------------------------------------------*/
	  c2_liste=(struct cell_trap *)malloc(sizeof(struct cell_trap));
	  c2_liste->point=k;
	  c2_liste->est_sommet=0;
	  c2_liste->mult=1;
	  if(c1_liste->est_sommet)
	    {
	    c2_liste->a=(objet[k].v-objet[s].v)/(objet[k].u-objet[s].u);
	    c2_liste->b=objet[s].v-c2_liste->a*objet[s].u;
	    }
	  else
	    {
	    c2_liste->a=c1_liste->a;
	    c2_liste->b=c1_liste->b;
	    }
	  c2_liste->v=c2_liste->a*u+c2_liste->b;

	  /* Insertion de la cellule dans la liste 2 adequate.
	  ----------------------------------------------------*/
	  if(!*p_liste || (*p_liste)->v>c2_liste->v)
	    {c2_liste->svt=*p_liste; *p_liste=c2_liste;}
	  else
	    {
	    for(c_liste=*p_liste;
		c_liste->svt && c_liste->svt->v<c2_liste->v;
		c_liste=c_liste->svt);
	    c2_liste->svt=c_liste->svt;
	    c_liste->svt=c2_liste;
	    }
	  }
      }

  /* Il ne reste plus qu'a inserer les sommets de rang
     courant puis mettre a jour les autres listes 2.
  ----------------------------------------------------*/
  inserer_sommets(u,&i,&j);
  init_xi('A',2);
  init_xi('B',2);
  init_yi(2);
  return(u);
  }




/*===========================================================================
inserer_sommets | Insere dans les listes 'listeA2' et 'listeB2' associees a
		  un rang donne les sommets de ce rang (s'il y en a).
-----------------------------------------------------------------------------
in    | Le rang de l'un des points de l'un des objets.
inout | La multiplicite de chaque point insere est sommee dans la zone poin-
	tee par 'p_multA' ou 'p_multB' (selon l'origine dudit point). Le con-
	tenu de la zone avant l'appel constitue le premier terme de la somme.
-----------------------------------------------------------------------------
Lorsque 'inserer_sommets' est appelee, les listes 'listeA2' et 'listeB2' sont
supposees etre en cours de constitution, 'pc_tri' est suppose pointer sur la
premiere cellule de la liste de tri associee a un sommet dont le rang est su-
perieur ou egal a 'ucrt'. Apres l'appel, 'pc_tri' pointe sur la premiere cel-
lule de tri associee a un sommet de rang strictement superieur a 'ucrt'.
===========================================================================*/


static void inserer_sommets(double ucrt, int *p_multA, int *p_multB)
  {
  int i, j, *p_mult;
  struct cell_objet *objet;
  struct cell_trap *nc_liste, *c_liste;
  struct cell_trap *c_listeA, *c_listeB, **p_liste, **pc_liste;

  c_listeA=listeA2;
  c_listeB=listeB2;
  for(;pc_tri && pc_tri->u==ucrt;pc_tri=pc_tri->svt)
    {
    /* A quel objet est rattache le sommet traite ?
    -----------------------------------------------*/
    if(pc_tri>=objetA && pc_tri<objetA_max)
      {objet=objetA; p_liste=&listeA2; pc_liste=&c_listeA; p_mult=p_multA;}
    else
      {objet=objetB; p_liste=&listeB2; pc_liste=&c_listeB; p_mult=p_multB;}

    /* Creation puis initialisation de la cellule trapeze.
    ------------------------------------------------------*/
    nc_liste=(struct cell_trap *)malloc(sizeof(struct cell_trap));
    nc_liste->point=pc_tri-objet;
    nc_liste->est_sommet=1;
    nc_liste->v=pc_tri->v;
    i=pc_tri->pred;
    j=pc_tri->succ;
    if(objet[i].u<ucrt && objet[j].u<ucrt) nc_liste->mult=2;
    else if(objet[i].u>=ucrt && objet[j].u>=ucrt) nc_liste->mult=0;
    else nc_liste->mult=1;
    *p_mult+=nc_liste->mult;

    /* Insertion dans la liste trapeze.
    -----------------------------------*/
    c_liste=*pc_liste;
    if(!c_liste || c_liste->v>nc_liste->v)
      {
      nc_liste->svt=c_liste;
      *p_liste=*pc_liste=nc_liste;
      }
    else
      {
      for(;c_liste->svt && c_liste->svt->v<nc_liste->v;c_liste=c_liste->svt);
      nc_liste->svt=c_liste->svt;
      *pc_liste=c_liste->svt=nc_liste;
      }
    }
  }




/*===========================================================================
detruire | Destruction d'une liste de cellules trapezes.
-----------------------------------------------------------------------------
in | Pointeur sur la premiere cellule de la liste.
===========================================================================*/


static void detruire(struct cell_trap *liste)
  {
  struct cell_trap *c_liste;
  while(liste) {c_liste=liste; liste=liste->svt; free(c_liste);}
  }




/*===========================================================================
copie | Copie d'une liste de cellules trapezes.
-----------------------------------------------------------------------------
in  | Pointeur sur la premiere cellule de la liste a copier.
out | Pointeur sur la premiere cellule de la copie.
===========================================================================*/


static struct cell_trap *copie(struct cell_trap *liste1)
  {
  struct cell_trap *liste2, *c_liste1, *c_liste2;

  if(!liste1) return(NULL);
  liste2=(struct cell_trap *)malloc(sizeof(struct cell_trap));
  *liste2=*liste1;
  c_liste1=liste1->svt;
  c_liste2=liste2;
  while(c_liste1)
    {
    c_liste2->svt=(struct cell_trap *)malloc(sizeof(struct cell_trap));
    *c_liste2->svt=*c_liste1;
    c_liste1=c_liste1->svt;
    c_liste2=c_liste2->svt;
    }
  c_liste2->svt=NULL;
  return(liste2);
  }




/*===========================================================================
force_trapezes | Calcul de la resultante des forces d'interactions
		 entre les deux ensembles de trapezes definis par
		 'xiA1', 'xiA2', 'xiB1', 'xiB2', 'yiA1B1', 'yiA2B2'.
-----------------------------------------------------------------------------
out | La force resultante (non corrigee par une valeur multiplica-
      tive du type : puissance de e_theta * epsilon * constante).
-----------------------------------------------------------------------------
'force_trapezes' consulte les listes precedentes
ainsi que les variables 'nb_trapA' et 'nb_trapB'.
===========================================================================*/


static double force_trapezes()
  {
  int iA, iB;
  double f, g, h;
  struct cell_yi *yi1, *yi2;

  f=0.0;

  for(yi1=yiA1B1,yi2=yiA2B2,iA=0;iA<nb_trapA;iA++)
    for(iB=0;iB<nb_trapB;iB++,yi1++,yi2++)
      if(type_force==0.0) {
        /*----- Forces constantes, objets supposes disjoints. -----*/
	    if(yi1->b1>0.0) /* y1>0 (et donc par hypothese y2>0 aussi) */
		{
	    f+=(xiA1[iA]+xiA2[iA])*(xiB1[iB]+xiB2[iB]);
	    f+=yi1->a1;
	    f+=yi2->a1;
	    }
	  } else if(type_force==1.0) {
	    /*----- Forces en 1/d, objets supposes disjoints. -----*/
	    if(yi1->b1>0.0) /* y1>0 (et donc par hypothese y2>0 aussi) */
		{
	    g=yi2->a1-yi1->a1;
		if(fabs(g)<SEUIL_ZERO) {h=yi2->a1+yi1->a1; f-=h*(log(h)-CTE_F1);}
	    else f-=(yi2->a2-yi1->a2)/g;
	    /*--------------*/
	    g=yi2->b1-yi1->b1;
	    if(fabs(g)<SEUIL_ZERO) {h=yi2->b1+yi1->b1; f+=h*(log(h)-CTE_F1);}
	    else f+=(yi2->b2-yi1->b2)/g;
	    /*--------------*/
	    g=yi2->c1-yi1->c1;
	    if(fabs(g)<SEUIL_ZERO) {h=yi2->c1+yi1->c1; f-=h*(log(h)-CTE_F1);}
	    else f-=(yi2->c2-yi1->c2)/g;
	    /*--------------*/
	    g=yi2->d1-yi1->d1;
	    if(fabs(g)<SEUIL_ZERO) {h=yi2->d1+yi1->d1; f+=h*(log(h)-CTE_F1);}
	    else f+=(yi2->d2-yi1->d2)/g;
	    }
	  } else if(type_force==2.0) {
	    /*----- Forces gravitationnelles, objets supposes disjoints. -----*/
	    if(yi1->b1>0.0) /* y1>0 (et donc par hypothese y2>0 aussi) */
		{
	    g=yi2->a1-yi1->a1;
	    if(fabs(g)<SEUIL_ZERO) f+=log(yi2->a1+yi1->a1)+CTE_F2;
	    else f+=(yi2->a2-yi1->a2)/g;
	    /*--------------*/
	    g=yi2->b1-yi1->b1;
	    if(fabs(g)<SEUIL_ZERO) f-=log(yi2->b1+yi1->b1)+CTE_F2;
	    else f-=(yi2->b2-yi1->b2)/g;
	    /*--------------*/
	    g=yi2->c1-yi1->c1;
	    if(fabs(g)<SEUIL_ZERO) f+=log(yi2->c1+yi1->c1)+CTE_F2;
	    else f+=(yi2->c2-yi1->c2)/g;
	    /*--------------*/
	    g=yi2->d1-yi1->d1;
	    if(fabs(g)<SEUIL_ZERO) f-=log(yi2->d1+yi1->d1)+CTE_F2;
	    else f-=(yi2->d2-yi1->d2)/g;
	    }
	  } else if(type_force==3.0) {
	  	/*----- Dernier cas particulier, objets supposes disjoints. -----*/
	    if(yi1->b1>0.0) /* y1>0 (et donc par hypothese y2>0 aussi) */
		{
	    g=yi2->a1-yi1->a1;
	    if(fabs(g)<SEUIL_ZERO) f-=2.0/(yi2->a1+yi1->a1);
	    else f-=(yi2->a2-yi1->a2)/g;
	    /*--------------*/
	    g=yi2->b1-yi1->b1;
	    if(fabs(g)<SEUIL_ZERO) f+=2.0/(yi2->b1+yi1->b1);
	    else f+=(yi2->b2-yi1->b2)/g;
	    /*--------------*/
	    g=yi2->c1-yi1->c1;
	    if(fabs(g)<SEUIL_ZERO) f-=2.0/(yi2->c1+yi1->c1);
	    else f-=(yi2->c2-yi1->c2)/g;
	    /*--------------*/
	    g=yi2->d1-yi1->d1;
	    if(fabs(g)<SEUIL_ZERO) f+=2.0/(yi2->d1+yi1->d1);
	    else f+=(yi2->d2-yi1->d2)/g;
	    }
	  } else {
		/*---- Cas general (mais objets toujours supposes disjoints). ----*/
	    if(yi1->b1>0.0) /* y1>0 (et donc par hypothese y2>0 aussi) */
		{
	    g=yi2->a1-yi1->a1;
		if(fabs(g)<SEUIL_ZERO)
 		  f-=(3.0-type_force)*pow((yi2->a1+yi1->a1)/2.0,2.0-type_force);
	    else f-=(yi2->a2-yi1->a2)/g;
	    /*--------------*/
	    g=yi2->b1-yi1->b1;
	    if(fabs(g)<SEUIL_ZERO)
		  f+=(3.0-type_force)*pow((yi2->b1+yi1->b1)/2.0,2.0-type_force);
	    else f+=(yi2->b2-yi1->b2)/g;
	    /*--------------*/
	    g=yi2->c1-yi1->c1;
	    if(fabs(g)<SEUIL_ZERO)
		  f-=(3.0-type_force)*pow((yi2->c1+yi1->c1)/2.0,2.0-type_force);
	    else f-=(yi2->c2-yi1->c2)/g;
	    /*--------------*/
	    g=yi2->d1-yi1->d1;
	    if(fabs(g)<SEUIL_ZERO)
		  f+=(3.0-type_force)*pow((yi2->d1+yi1->d1)/2.0,2.0-type_force);
	    else f+=(yi2->d2-yi1->d2)/g;
	    }
	  }

  return(f);
  }




/*===========================================================================
init_xi | Calcul des xi associes a une liste donnee.
-----------------------------------------------------------------------------
in | Le nom de l'objet ('A' ou 'B') auquel se rapporte
     la liste, le numero de celle-ci (1 ou 2).
-----------------------------------------------------------------------------
'init_xi' consultera la liste specifiee ('listeA1', 'listeA2', 'listeB1' ou
'listeB2') et l'une des variables 'nb_trapA' et 'nb_trapB' (le nombre d'xi
doit donc etre calcule avant l'appel) afin de reconstituer 'xiA1', 'xiA2',
'xiB1' ou 'xiB2'. ATTENTION : 'init_xi' alloue de l'espace memoire pour cet-
te liste d'xi mais ne s'occupe pas de liberer celui qui aurait pu etre reser-
ve anterieurement pour cette meme liste.
===========================================================================*/


static void init_xi(char nom_objet, int num_liste)
  {
  int i, n, k;
  double *xi, v;
  struct cell_trap *liste;

  /* De quel objet, de quelle liste s'agit-il ?
  ---------------------------------------------*/
  if(nom_objet=='A')
    {
    n=nb_trapA;
    xi=(double *)calloc(n,sizeof(double));
    if(num_liste==1) {xiA1=xi; liste=listeA1;}
    else {xiA2=xi; liste=listeA2;}
    }
  else
    {
    n=nb_trapB;
    xi=(double *)calloc(n,sizeof(double));
    if(num_liste==1) {xiB1=xi; liste=listeB1;}
    else {xiB2=xi; liste=listeB2;}
    }

  /* Calcul des xi.
  -----------------*/
  k=liste->mult;
  for(i=0;i<n;i++)
    {
    if(!k) {do liste=liste->svt; while(!liste->mult); k=liste->mult;}
    v=liste->v;
    k--;
    if(!k) {do liste=liste->svt; while(!liste->mult); k=liste->mult;}
    xi[i]=liste->v-v;
    k--;
    }
  }




/*===========================================================================
init_yi | Calcul des yi et autres parametres associes aux listes 1 ou 2.
-----------------------------------------------------------------------------
in | Le numero des listes auxquelles se referer (1 ou 2).
-----------------------------------------------------------------------------
'init_yi' consultera les listes specifiees ('listeA1' et 'listeB1', ou bien
'listeA2' et 'listeB2') et les variables 'nb_trapA' et 'nb_trapB' afin de
reconstituer 'yiA1B1' ou 'yiA2B2'. ATTENTION : 'init_yi' alloue de l'espace
memoire pour cette liste d'yi mais ne s'occupe pas de liberer celui qui au-
rait pu etre reserve anterieurement pour cette meme liste.
===========================================================================*/


static void init_yi(int num_liste)
  {
  int iA, iB, kA, kB;
  double aux=3.0-type_force; /* pour cas general */
  double v1A, v2A, v1B, v2B;
  struct cell_yi *yi;
  struct cell_trap *c_listeA, *c_listeB, *listeB;

  /*------------------------------
     De quelles listes s'agit-il ?
  --------------------------------*/

  yi=(struct cell_yi *)calloc(nb_trapA*nb_trapB,sizeof(struct cell_yi));
  if(num_liste==1) {yiA1B1=yi; c_listeA=listeA1; listeB=listeB1;}
  else {yiA2B2=yi; c_listeA=listeA2; listeB=listeB2;}

  /*------------------------------------
     Calcul des yi et autres parametres.
  --------------------------------------*/

  kA=c_listeA->mult;
  for(iA=nb_trapA;iA;iA--)
    {
    /* Parcours des trapezes de A.
    ------------------------------*/
    if(!kA)
      {
      do c_listeA=c_listeA->svt; while(!c_listeA->mult);
      kA=c_listeA->mult;
      }
    v1A=c_listeA->v;
    kA--;
    if(!kA)
      {
      do c_listeA=c_listeA->svt; while(!c_listeA->mult);
      kA=c_listeA->mult;
      }
    v2A=c_listeA->v;
    kA--;

    c_listeB=listeB;
    kB=c_listeB->mult;
    for(iB=nb_trapB;iB;iB--,yi++)
      {
      /* Parcours des trapezes de B.
      ------------------------------*/
      if(!kB)
	{
	do c_listeB=c_listeB->svt; while(!c_listeB->mult);
	kB=c_listeB->mult;
	}
      v1B=c_listeB->v;
      kB--;
      if(!kB)
	{
	do c_listeB=c_listeB->svt; while(!c_listeB->mult);
	kB=c_listeB->mult;
	}
      v2B=c_listeB->v;
      kB--;

      /* Le calcul des parametres.
      ----------------------------*/
      if(type_force==0.0) {
	    /*----- Forces constantes. -----*/
	    if((yi->b1=v1A-v2B)>0.0) /* si y>0 ... */
	      yi->a1=(v2A-v1A)*(v2B-v1B); /* ...calcul de xz */
	  } else if(type_force==1.0) {
	    /*----- Forces en 1/d. -----*/
		if((yi->b1=v1A-v2B)>0.0) { /* si y>0 ... */
	      yi->a1=v2A-v2B; /* x+y */
	      yi->c1=v1A-v1B; /* y+z */
	      yi->d1=v2A-v1B; /* x+y+z */
	      yi->a2=yi->a1*yi->a1*(log(yi->a1)-0.5);
	      yi->b2=yi->b1*yi->b1*(log(yi->b1)-0.5);
	      yi->c2=yi->c1*yi->c1*(log(yi->c1)-0.5);
	      yi->d2=yi->d1*yi->d1*(log(yi->d1)-0.5);
	    }
	  } else if(type_force==2.0) {
	    /*----- Forces gravitationnelles. -----*/
		if((yi->b1=v1A-v2B)>0.0) { /* si y>0 ... */
  	      yi->a1=v2A-v2B; /* x+y */
	      yi->c1=v1A-v1B; /* y+z */
	      yi->d1=v2A-v1B; /* x+y+z */
	      yi->a2=yi->a1*log(yi->a1);
	      yi->b2=yi->b1*log(yi->b1);
	      yi->c2=yi->c1*log(yi->c1);
	      yi->d2=yi->d1*log(yi->d1);
	    }
	  } else if(type_force==3.0) {
	    /*----- Dernier cas particulier. -----*/
		if((yi->b1=v1A-v2B)>0.0) { /* si y>0 ... */
  	      yi->a1=v2A-v2B; /* x+y */
	      yi->c1=v1A-v1B; /* y+z */
	      yi->d1=v2A-v1B; /* x+y+z */
	      yi->a2=log(yi->a1);
	      yi->b2=log(yi->b1);
	      yi->c2=log(yi->c1);
	      yi->d2=log(yi->d1);
	    }
	  } else {
	    /*----- Cas general. -----*/
		if((yi->b1=v1A-v2B)>0.0) { /* si y>0 ... */
  	      yi->a1=v2A-v2B; /* x+y */
	      yi->c1=v1A-v1B; /* y+z */
	      yi->d1=v2A-v1B; /* x+y+z */
	      yi->a2=pow(yi->a1,aux);
	      yi->b2=pow(yi->b1,aux);
	      yi->c2=pow(yi->c1,aux);
	      yi->d2=pow(yi->d1,aux);
	    }
	  }

      }
    }
  }








/****************************************************************************
							    DEBUG FUNCTIONS
****************************************************************************/






/*===========================================================================
debug_xi | Affichage de l'une des listes 'xiA1', 'xiA2', 'xiB1', 'xiB2'.
===========================================================================*/


#ifdef DEBUG
static void debug_xi(char nom_objet, int num_liste)
  {
  char *s;
  int i, n;
  double *x;

  if(nom_objet=='A')
    {
    n=nb_trapA;
    if(num_liste==1) {x=xiA1; s="xiA1";} else {x=xiA2; s="xiA2";}
    }
  else
    {
    n=nb_trapB;
    if(num_liste==1) {x=xiB1; s="xiB1";} else {x=xiB2; s="xiB2";}
    }

  printf("%s   : ",s);
  for(i=0;i<n;i++) {printf("%3.2lf ",x[i]); if(x[i]<0.0) printf("x<0!!! ");}
  printf("\n");
  }
#endif




/*===========================================================================
debug_yi | Affichage des yi de l'une des listes 'yiA1B1', 'yiA2B2'.
===========================================================================*/


#ifdef DEBUG
static void debug_yi(int num_liste)
  {
  int i;
  char *s;
  struct cell_yi *x;

  if(num_liste==1) {x=yiA1B1; s="yiA1B1";}
  else {x=yiA2B2; s="yiA2B2";}

  printf("%s : ",s);
  for(i=0;i<nb_trapA*nb_trapB;i++) printf("%3.2lf ",x[i].b1);
  printf("\n");
  }
#endif




/*===========================================================================
debug_cell_trap | Affichage des elements d'une liste 'cell_trap'.
===========================================================================*/


#ifdef DEBUG
static void debug_cell_trap (struct cell_trap *liste)
  {
  printf("DEBUT liste\n");
  while(liste) {
	  printf("v=%f est_sommet=%d point=%d mult=%d a=%f b=%f\n",
		  liste->v,liste->est_sommet,liste->point,
		  liste->mult,liste->a,liste->b);
	  liste=liste->svt;
  }
  printf("FIN liste\n");
}
#endif




/*===========================================================================
debug_pc_tri | Affichage des elements de la liste 'pc_tri'.
===========================================================================*/


#ifdef DEBUG
static void debug_pc_tri() {
	struct cell_objet* p=pc_tri;
	if(p) {
	     printf("u=%f v=%f x=%f y=%f pred=%d succ=%d\n",
			p->u, p->v, p->x, p->y, p->pred, p->succ);
	     p=p->svt;
	}
	while(p && p!=pc_tri) {
	     printf("u=%f v=%f x=%f y=%f pred=%d succ=%d\n",
			p->u, p->v, p->x, p->y, p->pred, p->succ);
	     p=p->svt;
	}
}
#endif








/****************************************************************************
							   GLOBAL FUNCTIONS
****************************************************************************/






/*===========================================================================
setTypeOfForce | Choix du type de force a considerer
		  (initialisation de la variable 'type_force').
-----------------------------------------------------------------------------
in | La valeur de r (0.0 pour forces constantes,
	 2.0 pour forces gravitationnelles, etc).
===========================================================================*/


void setTypeOfForce (double r)
  {
  double aux=SEUIL_ZERO*0.5;

  if(r>0.0-aux && r<0.0+aux) type_force=0.0;
  else if(r>1.0-aux && r<1.0+aux) type_force=1.0;
  else if(r>2.0-aux && r<2.0+aux) type_force=2.0;
  else if(r>3.0-aux && r<3.0+aux) type_force=3.0;
  else type_force=r;
  }




/*===========================================================================
readObject | Initialisation de 'objetA' et 'objetA_max', ou 'objetB'
             et 'objetB_max', a partir de donnees lues dans un fichier
			 (cf. Section I.)
-----------------------------------------------------------------------------
in | Le nom de l'objet ('A' ou 'B'), le nom du fichier qui le represente.
-----------------------------------------------------------------------------
L'ancien objet de meme nom est detruit. Seuls les champs 'x', 'y', 'pred'
et 'succ' sont initialises (les autres le seront par 'trier_sommets').
===========================================================================*/


void readObject (char nom_objet, char *nom_fichier)
  {
  FILE *f;
  int i, n, s, t;
  struct cell_objet *objet;

  f=fopen(nom_fichier,"rt");
  fscanf(f,"%d",&n);
  if(nom_objet=='A')
    {
	if(objetA) free(objetA);
    objet=objetA=(struct cell_objet *)calloc(n,sizeof(struct cell_objet));
    objetA_max=objetA+n;
    }
  else
    {
	if(objetB) free(objetB);
    objet=objetB=(struct cell_objet *)calloc(n,sizeof(struct cell_objet));
    objetB_max=objetB+n;
    }

  for(s=0;s<n;s=t)
    {
    fscanf(f,"%d",&i);
    t=s+i;
    for(i=s;i<t;i++)
      {
      fscanf(f,"%lf %lf",&objet[i].x,&objet[i].y);
      objet[i].pred=i-1;
      objet[i].succ=i+1;
      }
    objet[s].pred=t-1;
    objet[t-1].succ=s;
    }

  fclose(f);
  }




/*===========================================================================
forcesInDirection | Calcul de la resultante des forces exercees par les
               points de l'objet A sur ceux de l'objet B et tendant cha-
			   cune a deplacer B dans une direction predefinie.
-----------------------------------------------------------------------------
in  | La direction (en radians) : ce doit etre un element de ]-PI;PI].
out | La force resultante.
===========================================================================*/


double forcesInDirection (double theta)
  {
  double f, g, umax, u1, u2, e_theta;

  /* On se ramene a des ensembles d'objets trapezoidaux.
  ------------------------------------------------------*/
  f=0.0;
  e_theta=trier_sommets(theta,&u1,&umax);

  if(u1>=umax)
    /* Les deux objets ne coexistent pas dans la direction theta. */
    return(0.0);
  init_trapezes(u1);
  do
    {
    u2=maj_trapezes();
#ifdef DEBUG
g=(u2-u1)*force_trapezes();
f+=g;
printf("\nu1=%3.2lf u2=%3.2lf ",u1,u2);
if(u2<u1) printf("u2<u1 !!");
if(u2==u1) printf("u2==u1 !!");
printf("\n-----------------------\n");
debug_xi('A',1);
debug_xi('B',1);
debug_yi(1);
debug_xi('A',2);
debug_xi('B',2);
debug_yi(2);
#else
    f+=(u2-u1)*force_trapezes();
#endif
    u1=u2;
    }
  while(u2<umax);

  /* La force calculee doit etre corrigee par une constante
     multiplicative. Ne pas oublier de TOUJOURS multiplier
     par e_theta, puisque l'axe des u utilise pour trier
     les sommets n'est pas normalise.
  -----------------------------------*/
    if(type_force==0.0) {
      /* Forces constantes, objets supposes disjoints. Le facteur
         multiplicatif est : e_theta*[1/(6*e_theta*e_theta)] */
      return(f/(6.0*e_theta)); 
    }
    if(type_force==1.0) {
      /* Forces en 1/d, objets supposes disjoints.
         Le facteur multiplicatif est : e_theta*[1/(2*e_theta)] */
      return(f*0.5);
    }
    if(type_force==2.0) {
      /* Forces gravitationnelles, objets supposes disjoints.
         Le facteur multiplicatif est : e_theta*1 */
      return(f*e_theta);
    }
    if(type_force==3.0) {
      /* Dernier cas particulier, objets supposes disjoints.
         Le facteur multiplicatif est : e_theta*[e_theta/2] */
      return(f*e_theta*e_theta*0.5);
    }
    /* Cas general (mais objets toujours supposes disjoints). */
	g=pow(e_theta,type_force-1.0);
	g/=(1.0-type_force)*(2.0-type_force)*(3.0-type_force);
    return(f*g);
  }




/*===========================================================================
FRHistogram_CrispVector | Lire Sections I. et II. en debut de ce fichier.
===========================================================================*/


void FRHistogram_CrispVector (double *histogram,
                            int numberOfDirections, double typeOfForce,
                            char *argument, char *referent) {
  double step, angle;
  int index1, index2, i;
  setTypeOfForce(typeOfForce);
  readObject('A',argument);
  readObject('B',referent);
  histogram[numberOfDirections/2]=forcesInDirection(PI);
  histogram[0]=histogram[numberOfDirections]=forcesInDirection(0.0);
  step=2*PI/numberOfDirections; 
  for (i=1,angle=step-PI,index1=numberOfDirections/2+1,
	                  index2=numberOfDirections/2-1;
       index2;
       i++,angle=step*i-PI,index1++,index2--)
    {
      histogram[index1]=forcesInDirection(angle);
      histogram[index2]=forcesInDirection(-angle);
    }
} 








/*******************************************************************************
						    EXAMPLE OF main() FUNCTION
*******************************************************************************/






/*==============================================================================
   main | This function is provided as an example and for testing purposes.
--------------------------------------------------------------------------------
   Compile (e.g., "gcc FHistogramVector.c -lm"), execute ("a.out"),
   and follow the instructions.
==============================================================================*/


/*void main () {


  	FILE *histogramFile;
	char nameObjectA[100], nameObjectB[100], nameHistogramFile[100];
	double *histogram, typeForce;
	int i, numberOfDirections;


	/* Choosing the arguments.
	--------------------------*/

/*	printf("\nEnter arguments.");
	printf("\nYou are supposed to know the types and domains.");

	printf("\n\nHistogram will be stored in file: ");
	scanf("%s",nameHistogramFile);
	printf("Number of directions to be considered is: ");
	scanf("%d",&numberOfDirections);
	printf("Type of force is: ");
	scanf("%lf",&typeForce);
	printf("Argument object is defined by the text file: ");
	scanf("%s",nameObjectA);
	printf("Referent object is defined by the text file: ");
	scanf("%s",nameObjectB);
	

	/* Allocating memory for the histogram.
       Computing the histogram.
	---------------------------*/

/*	histogram=(double *)malloc((numberOfDirections+1)*sizeof(double));
    FRHistogram_CrispVector(histogram,numberOfDirections,typeForce,
		                    nameObjectA,nameObjectB);


	/* Writing the histogram file.
	------------------------------*/

/*	if((histogramFile=fopen(nameHistogramFile,"wt"))==NULL)
		{printf("\nERROR histogram file\n\n"); exit(1);}
	fprintf(histogramFile,"%d\n",numberOfDirections);
	for(i=0;i<=numberOfDirections;i++)
		fprintf(histogramFile,"%f\n",histogram[i]);
	fclose(histogramFile);


	/* Done.
	--------*/

/*	free(histogram);
}

*/
