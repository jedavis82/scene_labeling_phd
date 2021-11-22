/*
 * HoF_Vector.cpp
 *
 *  Created on: Dec 13, 2010
 *      Author: Ozy_2
 */

#include "HoF_Vector.hpp"

/*===============================================================================
 int computeHistogram (double *histogram, int scheme, int numberOfDirections,
                       double r, struct file_object *argA, struct file_object *refB)
 ---------------
 Computation of resultant forces in a number of evenly distributed directions.
 ---------------
 in     | - scheme: 1 for the single sum scheme, 2 for the double sum scheme
        | - numberOfDirections: Any positive integer multiple of 4
        |   (to cover the horizontal and vertical directions)
        | - r: The type of force, i.e., any real number
        |   (e.g., 0.0 for the computation of a histogram of constant forces,
        |   2.0 for the computation of a histogram of gravitational forces)
        | - argA: A pointer to the argument object (as returned by 'getObject')
        | - refB: A pointer to the referent object (as returned by 'getObject')
 ---------------
 out    | - histogram: The histogram of forces. 'numberOfDirections+1' values are stored
        |   in the array 'histogram'. 'histogram[0]' is the resultant force in direction
        |   0 (to the right) and 'histogram[numberOfDirections]' is the resultant force
        |   in direction 2PI (equal to 'histogram[0]').
        |   NOTE: The memory space pointed by 'histogram'
        |         must be allocated by the calling function.
 ---------------
 return | 0 if only finite forces are encountered, -1 if infinite
        | forces are encountered (in which case the values stored
        | in the array 'histogram' are arbitrary).
 =================================================================================*/


int hof::HoF_Vector::computeHistogram (	double *histogram, int scheme, int numberOfDirections, double typeOfForce,
										struct file_object *argument, struct file_object *referent)
{
	double step, angle;
	int index1, index2, i;

	setTypeOfForce(typeOfForce);
	setArgAndRef(argument,referent);

	if(scheme==1)
	{
		/*
		 * Single Scheme
		 */

		if((histogram[numberOfDirections/2]=computeForce_SingleScheme(PI))<0) return(-1);
		if((histogram[0]=histogram[numberOfDirections]=computeForce_SingleScheme(0.0))<0) return(-1);
		step=2*PI/numberOfDirections;
		for (i=1,angle=step-PI,index1=numberOfDirections/2+1,
			 index2=numberOfDirections/2-1;
			 index2;
			 i++,angle=step*i-PI,index1++,index2--)
	    {
			if((histogram[index1]=computeForce_SingleScheme(angle))<0) return(-1);
			if((histogram[index2]=computeForce_SingleScheme(-angle))<0) return(-1);
	    }
	}
	else
	{
		/*
		 * Double Scheme
		 */

		if((histogram[numberOfDirections/2]=computeForce_DoubleScheme(PI))<0) return(-1);
		if((histogram[0]=histogram[numberOfDirections]=computeForce_DoubleScheme(0.0))<0) return(-1);
		step=2*PI/numberOfDirections;
		for (i=1,angle=step-PI,index1=numberOfDirections/2+1,
			 index2=numberOfDirections/2-1;
			 index2;
			 i++,angle=step*i-PI,index1++,index2--)
	    {
			if((histogram[index1]=computeForce_DoubleScheme(angle))<0) return(-1);
			if((histogram[index2]=computeForce_DoubleScheme(-angle))<0) return(-1);
	    }
	}

	return(0);
}

/*
 *
 */

int hof::HoF_Vector::computeHistogram( 	std::vector<double>& histogram, int scheme, int numberOfDirection, double typeOfForce,
										HoF_Vector_Object& argument, HoF_Vector_Object& referent )
{
	struct file_object *arg = argument.createFileObjectStruct();
	struct file_object *ref = referent.createFileObjectStruct();

	histogram.clear();
	histogram.assign( numberOfDirection+1, 0.0 );

	int err = computeHistogram ( &histogram[0], scheme, numberOfDirection, typeOfForce, arg, ref);

	freeObject( arg );
	freeObject( ref );
}

/*
 *
 ---------------------------------------------------------------------------------------

	 double computeForce_SingleScheme (double theta)
	 ---------------
	 Calculation of the sum of all the elementary forces of type r exerted by the points
	 of A on those of B and that tend to move B in direction theta. The objects A and B
	 must first be set by 'setArgAndRef', and r must first be set by 'setTypeOfForce'.
	 The single sum scheme is used.
	 ---------------
	 in     | The direction theta, in radians. It must belong to (-PI;PI].
	 ---------------
	 return | The resultant force, or -1.0 if the resultant force is infinite.

 ---------------------------------------------------------------------------------------
 */
double hof::HoF_Vector::computeForce_SingleScheme (double theta)
{
	double f, g, umax, u1, u2, e_theta=0, aA=0, aB=0;
	int iA, iB, i, invalid;

	/* A and B are divided into trapezoids.
	 ---------------------------------------*/
	f=0.0;
	for(i=0;i<nalphas;i++){
		for(iA=0;iA<numA;iA++){
			if(obAlist[obAindex[iA]].fuzz<alphas[i]) continue; /* only consider A^i */
			else{aA=obAlist[obAindex[iA]].fuzz; break;}
		}
		for(iB=0;iB<numB;iB++){
			if(obBlist[obBindex[iB]].fuzz<alphas[i]) continue; /* only consider B^i */
			else{aB=obBlist[obBindex[iB]].fuzz; break;}
		}
		/* for each component in image A */
		for(iA=0;iA<numA;iA++){
			if(fabs(aA-obAlist[obAindex[iA]].fuzz)>SEUIL_ZERO) continue;
			copyOb('A',iA);
			/* for each component in image B */
			for(iB=0;iB<numB;iB++){
				if(fabs(aB-obBlist[obBindex[iB]].fuzz)>SEUIL_ZERO) continue;
				copyOb('B',iB);

				e_theta=trier_sommets(theta,&u1,&umax);

				if(u1<umax){
					/* The two objects do not coexist in direction theta. */
					init_trapezes(u1);
					do
					{
						u2=maj_trapezes();
						invalid=0;
#ifdef DEBUG
						g=(u2-u1)*force_trapezes(u1,u2,&invalid);
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
						if(fabs(u2-u1)>SEUIL_ZERO) f+=(alphas[i]-alphas[i+1])*(u2-u1)*force_trapezes(u1,u2,&invalid);
#endif
						detruire(listeA1); listeA1=NULL;
						detruire(listeB1); listeB1=NULL;
						if(invalid) return(-1.0);
						u1=u2;
					}while(u2<umax);

					detruire(listeA1); listeA1=NULL;
					detruire(listeB1); listeB1=NULL;
					detruire(listeA2); listeA2=NULL;
					detruire(listeB2); listeB2=NULL;
				}
				free(objetB);
			}
			free(objetA);
		}
	}
	/* The calculated force must be corrected by a multiplicative factor.
     Do not forget to always multiply by e_theta, since the u-axis
     used to sort the vertices is not normalized.
	 ---------------------------------------------------------------------*/
    if(type_force==0.0) {
		/* Constant forces. It is assumed that the objets are disjoint or equal.
         The multiplicative factor is: e_theta*[1/(6*e_theta*e_theta)] */
		return(f/(6.0*e_theta));
    }
    if(type_force==1.0) {
		/* Forces in 1/d. It is assumed that the objets are disjoint.
         The multiplicative factor is: e_theta*[1/(2*e_theta)] */
		return(f*0.5);
    }
    if(type_force==2.0) {
		/* Gravitational forces. It is assumed that the objets are disjoint.
         The multiplicative factor is: e_theta*1 */
		return(f*e_theta);
    }
    if(type_force==3.0) {
		/* Last particular case. It is assumed that the objets are disjoint.
         The multiplicative factor is: e_theta*[e_theta/2] */
		return(f*e_theta*e_theta*0.5);
    }
    /* General case (but it is still assumed that the objets are disjoint or equal). */
	g=pow(e_theta,type_force-1.0);
	g/=(1.0-type_force)*(2.0-type_force)*(3.0-type_force);
	return(f*g);
}

/*
 * ---------------------------------------------------------------------------------------

	 double computeForce_DoubleScheme (double theta)
	 ---------------
	 Calculation of the sum of all the elementary forces of type r exerted by the points
	 of A on those of B and that tend to move B in direction theta. The objects A and B
	 must first be set by 'setArgAndRef', and r must first be set by 'setTypeOfForce'.
	 The double sum scheme is used.
	 ---------------
	 in     | The direction theta, in radians. It must belong to (-PI;PI].
	 ---------------
	 return | The resultant force, or -1.0 if the resultant force is infinite.

 * ---------------------------------------------------------------------------------------
 */
double hof::HoF_Vector::computeForce_DoubleScheme (double theta)
{
	double f, g, umax, u1, u2, e_theta=0, aA=0, aB=0;
	int iA, iB, i, j, invalid;

	/* A and B are divided into trapezoids.
	 ---------------------------------------*/
	f=0.0;
	for(i=0;i<nalphas;i++){
		for(iA=0;iA<numA;iA++){
			if(obAlist[obAindex[iA]].fuzz<alphas[i]) continue; /* only consider A^i */
			else{aA=obAlist[obAindex[iA]].fuzz; break;}
		}
		for(j=0;j<nalphas;j++){
			for(iB=0;iB<numB;iB++){
				if(obBlist[obBindex[iB]].fuzz<alphas[j]) continue; /* only consider B^j */
				else{aB=obBlist[obBindex[iB]].fuzz; break;}
			}
			/* for each component in image A */
			for(iA=0;iA<numA;iA++){
				if(fabs(aA-obAlist[obAindex[iA]].fuzz)>SEUIL_ZERO) continue;
				copyOb('A',iA);
				/* for each component in image B */
				for(iB=0;iB<numB;iB++){
					if(fabs(aB-obBlist[obBindex[iB]].fuzz)>SEUIL_ZERO) continue;
					copyOb('B',iB);

					e_theta=trier_sommets(theta,&u1,&umax);

					if(u1<umax){
						/* The two objects do not coexist in direction theta. */
						init_trapezes(u1);
						do
						{
							u2=maj_trapezes();
							invalid=0;
#ifdef DEBUG
							g=(alphas[i]-alphas[i+1])*(alphas[j]-alphas[j+1])*(u2-u1)*force_trapezes(u1,u2,&invalid);
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
							if(fabs(u2-u1)>SEUIL_ZERO)f+=(alphas[i]-alphas[i+1])*(alphas[j]-alphas[j+1])*(u2-u1)*force_trapezes(u1,u2,&invalid);
#endif
							detruire(listeA1); listeA1=NULL;
							detruire(listeB1); listeB1=NULL;
							if(invalid) return(-1.0);
							u1=u2;
						}while(u2<umax);

						detruire(listeA1); listeA1=NULL;
						detruire(listeB1); listeB1=NULL;
						detruire(listeA2); listeA2=NULL;
						detruire(listeB2); listeB2=NULL;
					}
					free(objetB);
				}
				free(objetA);
			}
		}
	}
	/* The calculated force must be corrected by a multiplicative factor.
     Do not forget to always multiply by e_theta, since the u-axis
     used to sort the vertices is not normalized.
	 ---------------------------------------------------------------------*/
    if(type_force==0.0) {
		/* Constant forces. It is assumed that the objets are disjoint or equal.
         The multiplicative factor is: e_theta*[1/(6*e_theta*e_theta)] */
		return(f/(6.0*e_theta));
    }
    if(type_force==1.0) {
		/* Forces in 1/d. It is assumed that the objets are disjoint.
         The multiplicative factor is: e_theta*[1/(2*e_theta)] */
		return(f*0.5);
    }
    if(type_force==2.0) {
		/* Gravitational forces. It is assumed that the objets are disjoint.
         The multiplicative factor is: e_theta*1 */
		return(f*e_theta);
    }
    if(type_force==3.0) {
		/* Last particular case. It is assumed that the objets are disjoint.
         The multiplicative factor is: e_theta*[e_theta/2] */
		return(f*e_theta*e_theta*0.5);
    }
    /* General case (but it is still assumed that the objets are disjoint or equal). */
	g=pow(e_theta,type_force-1.0);
	g/=(1.0-type_force)*(2.0-type_force)*(3.0-type_force);
	return(f*g);
}


/*
 * ---------------------------------------------------------------------------------------

	void setTypeOfForce (double r)
	---------------
	Sets the type of force r.
	---------------
	in     | Type of force (e.g., 0.0 for constant forces, 2.0
		   | for gravitational forces). It can be any real number.

 *---------------------------------------------------------------------------------------
 */
void hof::HoF_Vector::setTypeOfForce (double r)
{
	double aux=SEUIL_ZERO*0.5;

	if(r>0.0-aux && r<0.0+aux) type_force=0.0;
	else if(r>1.0-aux && r<1.0+aux) type_force=1.0;
	else if(r>2.0-aux && r<2.0+aux) type_force=2.0;
	else if(r>3.0-aux && r<3.0+aux) type_force=3.0;
	else type_force=r;
}

/*
 *---------------------------------------------------------------------------------------

 void setArgAndRef (struct file_object *argA, struct file_object *refB)
 ---------------
 Sets the argument object A and the referent object B.
 ---------------
 in     | Pointers (as returned by 'getObject') to the two objects in memory.

 *---------------------------------------------------------------------------------------
 */
void hof::HoF_Vector::setArgAndRef(struct file_object *arg, struct file_object *ref)
{
	int i, j;
	double *a;

	obAlist=arg->ob;
	obAindex=arg->index;
	numA=arg->numob;

	obBlist=ref->ob;
	obBindex=ref->index;
	numB=ref->numob;

	nalphas=arg->nalphas+ref->nalphas;

	free(alphas);
	alphas=(double *)malloc(sizeof(double)*nalphas);
	a=arg->alphas;
	for(i=0,j=0,a=arg->alphas;i<nalphas;i++,j++){
		if(i==arg->nalphas){j=0;a=ref->alphas;}
		alphas[i]=a[j];
	}

	sortAlphas();
}


 /*
  * ---------------------------------------------------------------------------------------

	int writeHistogram (char *histogramFileName, double *histogram, int numberOfDirections,
						char *argObjectName, char *refObjectName, char *op)
	---------------
	Outputs a force histogram to a text file.
	---------------
	in     | Name of the output file, force histogram (as calculated by
		   | 'computeHistogram'), number of evenly distributed directions,
		   | names of the two objects, and the option for fopen ("w" or "a").
	---------------
	return | 0 if the histogram was written correctly, 1 otherwise.
	---------------
	  NOTE | The first line of the output file contains the name of the argument object,
		   | the name of the referent object, the number of directions, the type of force;
		   | the next lines contain the histogram values in directions 0 to 2PI.

 *	---------------------------------------------------------------------------------------
 */

int hof::HoF_Vector::writeHistogram(char *nameHistogramFile, double *histogram, int numberOfDirections, char *nameObjectA, char *nameObjectB, char *op)
{
	int i;
	FILE *histogramFile;

	if((histogramFile=fopen(nameHistogramFile,op))==NULL)
		return(1);

	fprintf(histogramFile,"%s, %s, %d, %.2lf\n", nameObjectA, nameObjectB, numberOfDirections, type_force);
	for(i=0;i<=numberOfDirections;i++)
		fprintf(histogramFile,"%f\n",histogram[i]);
	fprintf(histogramFile,"\n");
	fclose(histogramFile);
	return(0);
}

/*
 * ---------------------------------------------------------------------------------------

	 int readObjectFiles (struct file_object **objectPointers,
						  char** fileNames, int numberOfFiles)
	 ---------------
	 Reads in objects as described in files, and stores them in memory.
	 ---------------
	 in     | An array of filenames, and the number of files to be read.
	 ---------------
	 out    | An array of pointers to the objects in memory.
			| NOTE: Memory for this array must be allocated
			|       by the calling function.
	 ---------------
	 return | The number of successfully read files.
			| NOTE: The function returns either after
			|       all files have been read or after
			|       it encounters a file it cannot open.

 * ---------------------------------------------------------------------------------------
 */
int hof::HoF_Vector::readObjectFiles(struct file_object **file, char** filenames, int nfiles)
{
	FILE *f;
	int i, j, k, n=0, s, t, *index=NULL, *num, numa;
	struct cell_objet *objet=NULL;

	/* for each file */
	for(k=0;k<nfiles;k++){
		n=0;
		f=fopen(filenames[k],"rt");
		if(f==NULL) return(k);

		file[k]=(struct file_object *)malloc(sizeof(struct file_object));
		(file[k])->name=(char *)malloc(sizeof(char)*(strlen(filenames[k])+1));
		strcpy((file[k])->name, filenames[k]);

		int fscanf_ret;

		/* Read # of alpha cut */
		fscanf_ret = fscanf(f, "%d", &(file[k])->nalphas);

		numa=(file[k])->nalphas;
		(file[k])->alphas=(double *)malloc(sizeof(double)*numa);

		/* scan in alphas */
		for(i=0;i<numa;i++){
			fscanf_ret = fscanf(f, "%lf", &((file[k])->alphas)[i]);
		}

		s=0;
		num=&((file[k])->numob);
		*num=0;

		/**************************************************************************************************************
		 * VAR GLOSSARY:
		 *
		 * Info: an object may contain layers of alpha cuts.
		 *       Each alpha cut may be described using multiple polygons.
		 *       Each polygon in an alpha cut is described as an ordered list of vertices given as an (x,y) coordinate.
		 *
		 * nfiles = number of files to read (each file contains the definition of one object)
		 * k      = index to file array 'file', [0:k:nfiles-1]
		 * num    = numob
		 * numa   = nalphas ( number of alpha values )
		 * n      = total number of vertices in the object ( sum of number of vertices from all alpha values )
		 * i      = total number of vertices in a polygon for a given alpha.
		 * s      = accumulator for keeping the total number of vertices in the object as we read
		 */

		/* for each alpha */
		for(j=0;j<numa;j++)
		{
			/*
			 * Read # of vertices in this alpha cut,
			 * it is the sum of vertices in all polygons that make the object
			 */
			fscanf_ret = fscanf(f,"%d",&i);
			n+=i;
			/* Reallocate the array of all vertices in the object */
			objet=(file[k])->ob=(struct cell_objet *)realloc(objet,sizeof(struct cell_objet)*n);
			index=(file[k])->index=(int *)realloc(index,sizeof(int)*(n+1));

			index[0]=0;

			/* for each polygon in the alpha cut */
			for(;s<n;s=t)
			{
				/* Read # of vertices in this polygon */
				fscanf_ret = fscanf(f,"%d",&i);

				/*
				 * i = vertex count for this polygon
				 * s is set to equal t after each iteration
				 * s keeps track of the vertex count from all polygons in this alpha cut
				 * eventually t == s == n
				 */
				t=s+i;

				/*
				 * 's' still holds the vertex count from previously read polygons,
				 *   i.e., it is a starting index for inserting new vertices for the current polygon.
				 * 't' holds the vertex count from all polygons (both current and previous).
				 * 'j' is the index to alpha list
				 * 'objet' is a pointer to the beginning of an array of cell_objet structure used to store vertex info for an alpha cut.
				 * 'index' is an array of start and end index for each polygon in a given alhpa cuts.
				 * index[0] is always zero (meaning the first vertex in the first polygon is always stored first in objet array).
				 * index[1] stores the number of vertices of the 1st polygon.
				 * So we can access the vertices of the 1st polygon using objet[i] where index[0] <= i < index[1].
				 * index[2] is obtained by adding the number of vertices of the 2nd polygon to index[1].
				 * So we can access the vertices of the 2nd polygon using objet[i] where index[1] <= i < index[2].
				 * So the vertices of the k-th polygon are accessible using objet[i] where index[k-1] <= i < index[k]
				 * This is the reason why index must be allocated for n+1 where n is the total number of vertices for the alpha cut.
				 */
				for(i=s;i<t;i++)
				{
					/* Read (x,y) of this vertex */
					fscanf_ret = fscanf(f,"%lf %lf",&objet[i].x,&objet[i].y);

					/*
					 * check if next component is contained in last component
					 * If current vertex is contained, do what???
					 * and decrement the total count for polygon in alpha cut.
					 * In other words, do we reject such polygon? But why do we still
					 * add its vertices to the list? Moreover, since we do not
					 * terminate this loop if isContained returns 1, then we are
					 * running a serious risk of reducing
					 */
					if(i!=0 && i==s &&
						((file[k])->alphas)[j]==objet[i-1].fuzz &&
						isContained(objet[i].x,objet[i].y, *num, index, objet))
						(*num)--;
					objet[i].fuzz=((file[k])->alphas)[j];
					objet[i].pred=i-1;
					objet[i].succ=i+1;
				}
				objet[s].pred=t-1;
				objet[t-1].succ=s;
				index[++(*num)]=t; /* update index array */
			}

			index=(int *)realloc(index,sizeof(int)*(*num+1));
		}

		objet=NULL;
		index=NULL;
		fclose(f);
	}
	return(nfiles);
}

/*
 * ---------------------------------------------------------------------------------------

	 struct file_object *getObject (struct file_object **objectPointers,
									int numberOfObjects,
									char *fileName)
	 ---------------
	 Returns a pointer to the object that is described in a given file.
	 ---------------
	 in     | An array of pointers to objects (as initialized by 'readObjectFiles'),
			| the number of objects, and the name of the file
			| that describes the desired object.
	 ---------------
	 return | A pointer to the object if the object is found, NULL otherwise.

 * ---------------------------------------------------------------------------------------
 */

struct file_object* hof::HoF_Vector::getObject(struct file_object **files, int nfiles, char *filename)
{
	int i;
	for(i=0;i<nfiles;i++)
		if(files[i] && strcmp((files[i])->name,filename)==0)return(files[i]);
	return(NULL);
}

/*
 * ---------------------------------------------------------------------------------------

	void freeObjects (struct file_object **objectPointers, int numberOfObjects)
	---------------
	Frees the memory allocated to objects.
	---------------
	in     | An array of pointers to objects (as initialized by 'readObjectFiles'),
		   | and the number of objects.
		   | NOTE: The memory allocated to each object is freed, and the
		   |       memory allocated to the array of pointers is freed as well.
 * ---------------------------------------------------------------------------------------
 */
void hof::HoF_Vector::freeObjects(struct file_object **files, int n)
{
	int i;

	if(files==NULL) return;

	for(i=0;i<n;i++){
		if(files[i]){
			free(files[i]->name);
			free(files[i]->ob);
			free(files[i]->index);
			free(files[i]->alphas);
			free(files[i]);
		}
	}
	free(files);
}

/*
 *
 */

void hof::HoF_Vector::freeObject(struct file_object *files)
{
	if(files==NULL) return;

	free(files->name);
	free(files->ob);
	free(files->index);
	free(files->alphas);
	free(files);
}

/*
 * PRIVATE methods
 */

/*
 *
 */
double hof::HoF_Vector::trier_sommets(double theta, double *umin, double *umax)
{
	int i, param_tri2;
	double f, g, param_tri1, e_theta, uAmin, uAmax, uBmin, uBmax;
	struct cell_objet *sommet_crt, *sommet_fin, *pc_tri_min=NULL;

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

/*
 *
 */
void hof::HoF_Vector::init_trapezes(double ucrt)
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
						c_liste=c_liste->svt)
					{
							;
					}
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

/*
 *
 */
double hof::HoF_Vector::maj_trapezes()
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
							c_liste=c_liste->svt)
						{
							;
						}
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

/*
 *
 */
void hof::HoF_Vector::inserer_sommets(double ucrt, int *p_multA, int *p_multB)
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
			for(;c_liste->svt && c_liste->svt->v<nc_liste->v;c_liste=c_liste->svt)
			{
				;
			}
			nc_liste->svt=c_liste->svt;
			*pc_liste=c_liste->svt=nc_liste;
		}
    }
}


/*
 *
 */
void hof::HoF_Vector::detruire(struct cell_trap *liste)
{
	struct cell_trap *c_liste;
	while(liste) {c_liste=liste; liste=liste->svt; free(c_liste);}
}

/*
 *
 */
struct cell_trap* hof::HoF_Vector::copie(struct cell_trap *liste1)
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


/*=================================================================================
 overlapping_force_trapezes | Calculation of the the forces between two overlapping
 trapezoids passed in as yi1 and yi2.
 -----------------------------------------------------------------------------------
 in  | Current rank of the trapezoids, the two 'cell_yi' structures associated with
 the trapezoids.
 out | The resultant force (not corrected by a multiplicative factor of the form:
 power of e_theta * epsilon * constant).
 ----------------------------------------------------------------------------------
 'overlapping_force_trapezes' splits the trapezoids into sections based on the
 intersection of their non parallel sides. These sections are then reduced further
 into smaller trapezoids. These smaller trapezoids are either disjoint, touching,
 or equal.
 ==================================================================================*/

double hof::HoF_Vector::overlapping_force_trapezes(double u1, double u2, struct cell_yi *yi1, struct cell_yi *yi2)
{
	double f=0, f_temp, g, h;
	double ints[5]; /* array to hold intersects (4 possible intersects + 1 space for u2) */
	int i=0, j, k, l;
	double u=u2-u1;
	double nu1, nu2; /* upper and lower bounds for sections */
	double lista[6][2], listb[6][2]; /* lists to hold trapezoid sections */
	double slista[6][3][2], slistb[6][3][2]; /* lists to hold trapezoid sections after splitting */
	double x1, x2, z1, z2, y1, y2, xyz1, xyz2;
	double aux=3.0-type_force; /* for the general case */

	/* find intersects */
	if((g=yi1->x1-yi1->z1)*(h=yi2->x1-yi2->z1)<0){ /* intersect */
		ints[i++]=fabs(g)*u/(fabs(g)+fabs(h))+u1;
	}
	if((g=yi1->x1-yi1->z2)*(h=yi2->x1-yi2->z2)<0){ /* intersect */
		ints[i++]=fabs(g)*u/(fabs(g)+fabs(h))+u1;
	}
	if((g=yi1->x2-yi1->z1)*(h=yi2->x2-yi2->z1)<0){ /* intersect */
		ints[i++]=fabs(g)*u/(fabs(g)+fabs(h))+u1;
	}
	if((g=yi1->x2-yi1->z2)*(h=yi2->x2-yi2->z2)<0){ /* intersect */
		ints[i++]=fabs(g)*u/(fabs(g)+fabs(h))+u1;
	}
	ints[i]=u2;

	/* sort intersects by ascending y values */
	for(j=0;j<i;j++){
		g=ints[j];
		for(k=j+1;k<i;k++){
			if(ints[k]<g){
				g=ints[k];
				ints[k]=ints[j];
				ints[j]=g;
			}
		}
	}

	/* split trapezoids sections based on intersects */
	lista[0][0]=yi1->x1;
	lista[0][1]=yi1->x2;
	listb[0][0]=yi1->z1;
	listb[0][1]=yi1->z2;
	for(j=1;j<=i;j++){
		lista[j][0]=(yi2->x1-yi1->x1)/u*(ints[j-1]-u1)+yi1->x1;
		lista[j][1]=(yi2->x2-yi1->x2)/u*(ints[j-1]-u1)+yi1->x2;
		listb[j][0]=(yi2->z1-yi1->z1)/u*(ints[j-1]-u1)+yi1->z1;
		listb[j][1]=(yi2->z2-yi1->z2)/u*(ints[j-1]-u1)+yi1->z2;
	}
	lista[j][0]=yi2->x1;
	lista[j][1]=yi2->x2;
	listb[j][0]=yi2->z1;
	listb[j][1]=yi2->z2;

	/* split overlapping sections into smaller trapezoids */
	for(j=0;j<=i+1;j++){
		slista[j][0][0]=lista[j][0];
		slistb[j][0][0]=listb[j][0];
		if(listb[j][1]<lista[j][0]){ /* lines are disjoint, A is ahead of B in direction theta */
			slista[j][0][1]=lista[j][0];
			slista[j][1][0]=slista[j][1][1]=lista[j][0];
			slista[j][2][0]=lista[j][0];
			slista[j][2][1]=lista[j][1];

			slistb[j][0][1]=listb[j][1];
			slistb[j][1][0]=slistb[j][1][1]=listb[j][1];
			slistb[j][2][0]=slistb[j][2][1]=listb[j][1];
		}
		else if(listb[j][0]>lista[j][1]){ /* lines are disjoint, B is ahead of A in direction theta */
			slista[j][0][1]=lista[j][1];
			slista[j][1][0]=slista[j][1][1]=lista[j][1];
			slista[j][2][0]=slista[j][2][1]=lista[j][1];

			slistb[j][0][1]=listb[j][0];
			slistb[j][1][0]=slistb[j][1][1]=listb[j][0];
			slistb[j][2][0]=listb[j][0];
			slistb[j][2][1]=listb[j][1];
		}
		else{ /* lines overlap */
			slista[j][0][1]=((listb[j][0]+SEUIL_ZERO)>=lista[j][0] && listb[j][0]<=(lista[j][1]+SEUIL_ZERO))?listb[j][0]:lista[j][0];
			slista[j][1][0]=slista[j][0][1];
			slista[j][1][1]=((listb[j][1]+SEUIL_ZERO)>=slista[j][0][1] && listb[j][1]<=(lista[j][1]+SEUIL_ZERO))?listb[j][1]:lista[j][1];
			slista[j][2][0]=slista[j][1][1];
			slista[j][2][1]=lista[j][1];

			slistb[j][0][1]=((lista[j][0]+SEUIL_ZERO)>=listb[j][0] && lista[j][0]<=(listb[j][1]+SEUIL_ZERO))?lista[j][0]:listb[j][0];
			slistb[j][1][0]=slistb[j][0][1];
			slistb[j][1][1]=((lista[j][1]+SEUIL_ZERO)>=slistb[j][0][1] && lista[j][1]<=(listb[j][1]+SEUIL_ZERO))?lista[j][1]:listb[j][1];
			slistb[j][2][0]=slistb[j][1][1];
			slistb[j][2][1]=listb[j][1];
		}
	}


	/* calculate forces between new sections */
	nu1=u1;
	for(j=0;j<=i;j++){
		nu2=ints[j];
		f_temp=0.0;
		for(k=0;k<3;k++)
			for(l=0;l<3;l++){
				if(type_force==0){
					if((slista[j][k][0]-slistb[j][l][1]+SEUIL_ZERO)>=0.0 &&
					   (slista[j+1][k][0]-slistb[j+1][l][1]+SEUIL_ZERO)>=0.0){ /* disjoint or touching */
						x1=slista[j][k][1]-slista[j][k][0];
						x2=slista[j+1][k][1]-slista[j+1][k][0];
						z1=slistb[j][l][1]-slistb[j][l][0];
						z2=slistb[j+1][l][1]-slistb[j+1][l][0];

						f_temp+=(x1+x2)*(z1+z2);
						f_temp+=x1*z1;
						f_temp+=x2*z2;
					}
					else if(fabs(slista[j][k][0]-slistb[j][l][0])<SEUIL_ZERO &&
							fabs(slista[j][k][1]-slistb[j][l][1])<SEUIL_ZERO &&
							fabs(slista[j+1][k][0]-slistb[j+1][l][0])<SEUIL_ZERO &&
							fabs(slista[j+1][k][1]-slistb[j+1][l][1])<SEUIL_ZERO){ /* equal */
						x1=slista[j][k][1]-slista[j][k][0];
						x2=slista[j+1][k][1]-slista[j+1][k][0];

						f_temp+=x1*x2;
						f_temp+=x1*x1;
						f_temp+=x2*x2;
					}
				}
				else{
					/*---- General case (but it is still assumed that the objects are disjoint or equal). ----*/
					if((y1=(slista[j][k][0]-slistb[j][l][1]))>=-SEUIL_ZERO &&
					   (y2=(slista[j+1][k][0]-slistb[j+1][l][1]))>=-SEUIL_ZERO) /* disjoint or touching */
					{
						x1=slista[j][k][1]-slistb[j][l][1]; /* x+y */
						z1=slista[j][k][0]-slistb[j][l][0]; /* y+z */
						xyz1=slista[j][k][1]-slistb[j][l][0]; /* x+y+z */
						x2=slista[j+1][k][1]-slistb[j+1][l][1]; /* x+y */
						z2=slista[j+1][k][0]-slistb[j+1][l][0]; /* y+z */
						xyz2=slista[j+1][k][1]-slistb[j+1][l][0]; /* x+y+z */

						/* very small values are corected to 0 to avoid problems later */
						x1=(fabs(x1)<SEUIL_ZERO)?0:x1;
						z1=(fabs(z1)<SEUIL_ZERO)?0:z1;
						y1=(fabs(y1)<SEUIL_ZERO)?0:y1;
						xyz1=(fabs(xyz1)<SEUIL_ZERO)?0:xyz1;
						x2=(fabs(x2)<SEUIL_ZERO)?0:x2;
						z2=(fabs(z2)<SEUIL_ZERO)?0:z2;
						y2=(fabs(y2)<SEUIL_ZERO)?0:y2;
						xyz2=(fabs(xyz2)<SEUIL_ZERO)?0:xyz2;

						if((x1>SEUIL_ZERO || x2>SEUIL_ZERO) && (z1>SEUIL_ZERO || z2>SEUIL_ZERO)){ /* both trapezoids must have a non-zero area */
							g=x2-x1;
							if(fabs(g)<SEUIL_ZERO)
								f_temp-=(aux)*pow((x2+x1)/2.0,2.0-type_force);
							else f_temp-=(pow(x2,aux)-pow(x1,aux))/g;
							/*--------------*/
							g=y2-y1;
							if(fabs(g)<SEUIL_ZERO)
								f_temp+=(aux)*pow((y2+y1)/2.0,2.0-type_force);
							else f_temp+=(pow(y2,aux)-pow(y1,aux))/g;
							/*--------------*/
							g=z2-z1;
							if(fabs(g)<SEUIL_ZERO)
								f_temp-=(aux)*pow((z2+z1)/2.0,2.0-type_force);
							else f_temp-=(pow(z2,aux)-pow(z1,aux))/g;
							/*--------------*/
							g=xyz2-xyz1;
							if(fabs(g)<SEUIL_ZERO)
								f_temp+=(aux)*pow((xyz2+xyz1)/2.0,2.0-type_force);
							else f_temp+=(pow(xyz2,aux)-pow(xyz1,aux))/g;
						}
					}
					else if((type_force<1) &&
							fabs(slista[j][k][0]-slistb[j][l][0])<SEUIL_ZERO &&
							fabs(slista[j][k][1]-slistb[j][l][1])<SEUIL_ZERO &&
							fabs(slista[j+1][k][0]-slistb[j+1][l][0])<SEUIL_ZERO &&
							fabs(slista[j+1][k][1]-slistb[j+1][l][1])<SEUIL_ZERO) /* equal */
					{
						x1=slista[j][k][1]-slista[j][k][0];
						x2=slista[j+1][k][1]-slista[j+1][k][0];

						/* very small values are corected to 0 to avoid problems later */
						x1=(fabs(x1)<SEUIL_ZERO)?0:x1;
						x2=(fabs(x2)<SEUIL_ZERO)?0:x2;

						g=pow(x1,aux);
						h=pow(x2,aux);
						if(fabs(x2-x1)<SEUIL_ZERO)
							f_temp+=(aux)*pow(x1, 2-type_force);
						else
							f_temp+=(g-h)/(x1-x2);
					}
				}
			}
		f+=f_temp*(nu2-nu1);
		nu1=nu2;
	}
	return(f/u); /* divided by u to account for multiplication in computeForce */
}

/*=================================================================================
 force_trapezes | Calculation of the sum of all the forces between the two sets of
 trapezoids defined by 'xiA1', 'xiA2', 'xiB1', 'xiB2', 'yiA1B1', 'yiA2B2'.
 -----------------------------------------------------------------------------------
 in    | The current rank of the trapezoids.
 inout | The integer invalid is set to 1 if an invalid type of force was selected
 for the current objects, otherwise it is set to 0.
 out   | The resultant force (not corrected by a multiplicative factor of the form:
 power of e_theta * epsilon * constant).
 ----------------------------------------------------------------------------------
 'force_trapezes' takes a look at the lists described earlier and at the variables
 'nb_trapA' and 'nb_trapB'.
 ==================================================================================*/

double hof::HoF_Vector::force_trapezes(double u1, double u2, int *invalid)
{
	int iA, iB;
	double f, g, h;
	struct cell_yi *yi1, *yi2;

	*invalid=0;
	f=0.0;

	for(yi1=yiA1B1,yi2=yiA2B2,iA=0;iA<nb_trapA;iA++)
		for(iB=0;iB<nb_trapB;iB++,yi1++,yi2++)
			if(type_force==0.0) {
				/*----- Constant forces -----*/
				if(yi1->b1>=0.0 && yi2->b1>=0.0) /* y1>=0 and y2>=0 */
				{
					f+=(xiA1[iA]+xiA2[iA])*(xiB1[iB]+xiB2[iB]);
					f+=yi1->a1;
					f+=yi2->a1;
				}
				else if(fabs(yi1->b1+xiA1[iA])<SEUIL_ZERO && fabs(yi2->b1+xiA2[iA])<SEUIL_ZERO &&
						fabs(xiA1[iA]-xiB1[iB])<SEUIL_ZERO && fabs(xiA2[iA]-xiB2[iB])<SEUIL_ZERO) /* Objects are equal considering assumption. */
				{
					f+=(xiA1[iA])*(xiA2[iA]);
					f+=yi1->a1;
					f+=yi2->a1;
				}
				else if((yi1->b1+xiA1[iA]+xiB1[iB])>0 || (yi2->b1+xiA2[iA]+xiB2[iB])>0) /* trapezoids overlap */
				{
					f+=overlapping_force_trapezes(u1,u2,yi1,yi2);
				}
			} else if(type_force==1.0) {
				/*----- Forces in 1/d -----*/
				if(yi1->b1>=0.0 && yi2->b1>=0.0) /* y1>=0 and y2>=0 */
				{
					g=yi2->a1-yi1->a1;
					if(fabs(g)<SEUIL_ZERO) {h=yi2->a1+yi1->a1; f-=(h==0)?0:h*(log(h)-CTE_F1);}
					else f-=(yi2->a2-yi1->a2)/g;
					/*--------------*/
					g=yi2->b1-yi1->b1;
					if(fabs(g)<SEUIL_ZERO) {h=yi2->b1+yi1->b1; f+=(h==0)?0:h*(log(h)-CTE_F1);}
					else f+=(yi2->b2-yi1->b2)/g;
					/*--------------*/
					g=yi2->c1-yi1->c1;
					if(fabs(g)<SEUIL_ZERO) {h=yi2->c1+yi1->c1; f-=(h==0)?0:h*(log(h)-CTE_F1);}
					else f-=(yi2->c2-yi1->c2)/g;
					/*--------------*/
					g=yi2->d1-yi1->d1;
					if(fabs(g)<SEUIL_ZERO) {h=yi2->d1+yi1->d1; f+=(h==0)?0:h*(log(h)-CTE_F1);}
					else f+=(yi2->d2-yi1->d2)/g;
				}
				else if((yi1->b1<0 && (yi1->b1+xiA1[iA]+xiB1[iB])>SEUIL_ZERO) ||
						(yi2->b1<0 && (yi2->b1+xiA2[iA]+xiB2[iB])>SEUIL_ZERO)){ *invalid=1; return(0);}
			} else if(type_force==2.0) {
				/*----- Gravitational forces, and it is assumed that the objects are disjoint. -----*/
				if(yi1->b1>0.0 && yi2->b1>0) /* y1>0 (and therefore y2>0, considering the assumption) */
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
				else if((yi1->b1<=0 && (yi1->b1+xiA1[iA]+xiB1[iB])>SEUIL_ZERO) ||
						(yi2->b1<=0 && (yi2->b1+xiA2[iA]+xiB2[iB])>SEUIL_ZERO)){ *invalid=1; return(0);}

			} else if(type_force==3.0) {
				/*----- Last particular case, and it is assumed that the objects are disjoint. -----*/
				if(yi1->b1>0.0) /* y1>0 (and therefore y2>0, considering the assumption) */
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
				else if((yi1->b1<=0 && (yi1->b1+xiA1[iA]+xiB1[iB])>SEUIL_ZERO) ||
						(yi2->b1<=0 && (yi2->b1+xiA2[iA]+xiB2[iB])>SEUIL_ZERO)){ *invalid=1; return(0);}
			} else {
				/*---- General case ----*/
				if((yi1->b1>0.0 && yi2->b1>0.0) ||
				   (type_force<2 && yi1->b1==0.0 && yi2->b1==0.0))/* y1>=0 and y2>=0 */
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
				else if((type_force<1) && fabs(yi1->b1+xiA1[iA])<SEUIL_ZERO && fabs(yi2->b1+xiA2[iA])<SEUIL_ZERO &&
						fabs(xiA1[iA]-xiB1[iB])<SEUIL_ZERO && fabs(xiA2[iA]-xiB2[iB])<SEUIL_ZERO) /* Objects are equal. */
				{
					if(fabs(xiA1[iA]-xiA2[iA])<SEUIL_ZERO)
						f+=(3-type_force)*pow(xiA1[iA], 2-type_force);
					else
						f+=(yi1->a1-yi2->a1)/(xiA1[iA]-xiA2[iA]);
				}
				else if((type_force<1) && ((yi1->b1+xiA1[iA]+xiB1[iB])>0 || (yi2->b1+xiA2[iA]+xiB2[iB])>0)) /* trapezoids overlap */
				{
					f+=overlapping_force_trapezes(u1,u2,yi1,yi2);
				}
				else if((yi1->b1<=0 && (yi1->b1+xiA1[iA]+xiB1[iB])>SEUIL_ZERO) ||
						(yi2->b1<=0 && (yi2->b1+xiA2[iA]+xiB2[iB])>SEUIL_ZERO)){ *invalid=1; return(0);}
			}

	return(f);
}

/*
 *
 */
void hof::HoF_Vector::init_xi(char nom_objet, int num_liste)
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


/*=========================================================================================
 init_yi | Calculation of the yi cells associated with one of the two longitudinal sections.
 -------------------------------------------------------------------------------------------
 in | 1 or 2, depending on the section we are focusing on.
 -------------------------------------------------------------------------------------------
 'init_yi' takes a look at the relevant lists ('listeA1' and 'listeB1', or 'listeA2' and
 'listeB2') and at the variables 'nb_trapA' and 'nb_trapB' in order to calculate 'yiA1B1'
 or 'yiA2B2'. WARNING --- 'init_yi' allocates memory space for this list of yi cells but
 doesn't free the memory space that might have been previously allocated for this same list.
 ==========================================================================================*/

void hof::HoF_Vector::init_yi(int num_liste)
{
	int iA, iB, kA, kB;
	double aux=3.0-type_force; /* for the general case */
	double v1A, v2A, v1B, v2B;
	struct cell_yi *yi;
	struct cell_trap *c_listeA, *c_listeB, *listeB;

	/*----------------------------------
     Which lists are we talking about?
	 ------------------------------------*/

	yi=(struct cell_yi *)calloc(nb_trapA*nb_trapB,sizeof(struct cell_yi));
	if(num_liste==1) {yiA1B1=yi; c_listeA=listeA1; listeB=listeB1;}
	else {yiA2B2=yi; c_listeA=listeA2; listeB=listeB2;}

	/*-----------------------------
     Calculation of the yi cells
	 -------------------------------*/

	kA=c_listeA->mult;
	for(iA=nb_trapA;iA;iA--)
    {
		/* We scan the set of A-trapezoids
		 ----------------------------------*/
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
			/* We scan the set of B-trapezoids
			 ----------------------------------*/
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

			/* Calculation of the yi cell
			 ---------------------------------------------------------------
			 NOTE --- v1A and v2A are the coordinates on the longitudinal
			 section of the vertices of the A-trapezoid (and v2A is after
			 v1A in direction theta). v1B and v2B are the coordinates of
			 the vertices of the B-trapezoid.
			 --------------------------------------------------------------*/

			yi->x1=v1A;
			yi->x2=v2A;
			yi->z1=v1B;
			yi->z2=v2B;

			yi->b1=v1A-v2B;
			yi->b1=(fabs(yi->b1)<SEUIL_ZERO)?0:yi->b1;

			if(type_force==0.0) {
				/*----- Constant forces -----*/
				if((yi->b1)>=0.0){ /* if y>=0 ... */
					yi->a1=(v2A-v1A)*(v2B-v1B); /* ...calculation of xz */
				}
				else if(fabs(v1A-v1B)<SEUIL_ZERO && fabs(v2A-v2B)<SEUIL_ZERO){
					yi->a1=(v2A-v1A)*(v2A-v1A); /* ...calculation of x^2 */
				}
			} else if(type_force==1.0) {
				/*----- Forces in 1/d -----*/
				if((yi->b1)>=0.0) { /* if y>=0 ... */
					yi->a1=fabs(v2A-v2B)<SEUIL_ZERO?0:v2A-v2B; /* x+y */
					yi->c1=fabs(v1A-v1B)<SEUIL_ZERO?0:v1A-v1B; /* y+z */
					yi->d1=fabs(v2A-v1B)<SEUIL_ZERO?0:v2A-v1B; /* x+y+z */
					yi->a2=(yi->a1==0)?0:yi->a1*yi->a1*(log(yi->a1)-0.5);
					yi->b2=(yi->b1==0)?0:yi->b1*yi->b1*(log(yi->b1)-0.5);
					yi->c2=(yi->c1==0)?0:yi->c1*yi->c1*(log(yi->c1)-0.5);
					yi->d2=(yi->d1==0)?0:yi->d1*yi->d1*(log(yi->d1)-0.5);
				}
			} else if(type_force==2.0) {
				/*----- Gravitational forces -----*/
				if((yi->b1)>0.0) { /* if y>0 ... */
					yi->a1=fabs(v2A-v2B)<SEUIL_ZERO?0:v2A-v2B; /* x+y */
					yi->c1=fabs(v1A-v1B)<SEUIL_ZERO?0:v1A-v1B; /* y+z */
					yi->d1=fabs(v2A-v1B)<SEUIL_ZERO?0:v2A-v1B; /* x+y+z */
					yi->a2=yi->a1*log(yi->a1);
					yi->b2=yi->b1*log(yi->b1);
					yi->c2=yi->c1*log(yi->c1);
					yi->d2=yi->d1*log(yi->d1);
				}
			} else if(type_force==3.0) {
				/*----- Last particular case -----*/
				if((yi->b1)>0.0) { /* if y>0 ... */
					yi->a1=fabs(v2A-v2B)<SEUIL_ZERO?0:v2A-v2B; /* x+y */
					yi->c1=fabs(v1A-v1B)<SEUIL_ZERO?0:v1A-v1B; /* y+z */
					yi->d1=fabs(v2A-v1B)<SEUIL_ZERO?0:v2A-v1B; /* x+y+z */
					yi->a2=log(yi->a1);
					yi->b2=log(yi->b1);
					yi->c2=log(yi->c1);
					yi->d2=log(yi->d1);
				}
			} else {
				/*----- General case -----*/
				if((yi->b1)>=0.0) { /* if y>=0 ... */
					yi->a1=fabs(v2A-v2B)<SEUIL_ZERO?0:v2A-v2B; /* x+y */
					yi->c1=fabs(v1A-v1B)<SEUIL_ZERO?0:v1A-v1B; /* y+z */
					yi->d1=fabs(v2A-v1B)<SEUIL_ZERO?0:v2A-v1B; /* x+y+z */
					yi->a2=pow(yi->a1,aux);
					yi->b2=pow(yi->b1,aux);
					yi->c2=pow(yi->c1,aux);
					yi->d2=pow(yi->d1,aux);
				}
				else if(type_force<1 && fabs(v1A-v1B)<SEUIL_ZERO && fabs(v2A-v2B)<SEUIL_ZERO){
					yi->a1=fabs(v2A-v1A)<SEUIL_ZERO?0:pow((v2A-v1A), aux); /* x^(3-r) */
				}

			}

		}
    }
}

/*===========================================================================
 isContained | Determines if a coordinate is contained within a component by
 checking for vertical and horizontal intersects.
 -----------------------------------------------------------------------------
 in  | The coordinate values, the number of the component, the component index
 array, and the object list.
 out | 1 if the coordinate is contained in the component, 0 otherwise.
 ===========================================================================*/

int hof::HoF_Vector::isContained (double x, double y, int num, int *index, struct cell_objet *list)
{
	int i, up=0, down=0, right=0, left=0;
	double segv;

	/* for each edge in component */
	for(i=index[num-1];i<(index[num]);i++)
	{
		if((!up || !down) &&
		   (((x+SEUIL_ZERO)>=list[i].x && x<=(list[list[i].succ].x+SEUIL_ZERO)) ||
			(x<=(list[i].x+SEUIL_ZERO) && (x+SEUIL_ZERO)>=list[list[i].succ].x)))
		{
			/* coordinate must be under or over edge */
			if(fabs(list[list[i].succ].y-list[i].y)<SEUIL_ZERO)
			{
				/* edge is horizontal */
			    if(list[i].y>y) up=1;
			    else down=1;
			}
			else if(fabs(list[list[i].succ].x-list[i].x)>SEUIL_ZERO)
			{
				/* check for intersect */
				segv=(list[list[i].succ].y-list[i].y)/(list[list[i].succ].x-list[i].x)*(x-list[i].x)+list[i].y;
				if(segv>y) up=1;
				else down=1;
			 }
		}

		if((!right || !left) &&
		   (((y+SEUIL_ZERO)>=list[i].y && y<=(list[list[i].succ].y+SEUIL_ZERO)) ||
			(y<=(list[i].y+SEUIL_ZERO) && (y+SEUIL_ZERO)>=list[list[i].succ].y))){ /* coordinate mus be to the left or right of edge */
			   if(fabs(list[list[i].succ].x-list[i].x)<SEUIL_ZERO){ /* edge is vertical */
				   if(list[i].x>x) right=1;
				   else left=1;
			   }
			   else if(fabs(list[list[i].succ].y-list[i].y)>SEUIL_ZERO){ /* check for intersect */
				   segv=(y-list[i].y)/((list[list[i].succ].y-list[i].y)/(list[list[i].succ].x-list[i].x))+list[i].x;
				   if(segv>x) right=1;
				   else left=1;
			   }
		   }

	}
	if(up && down && right && left)
		return(1); /* must have four intersects */
	else
		return(0);
}

/*===========================================================================
 sortAlphas | Sorts the list of alphas and removes duplicates.
 ===========================================================================*/

void hof::HoF_Vector::sortAlphas()
{
	int i, j, n=1;
	double g;

	/* sort alphas in descending order */
	for(i=0;i<nalphas;i++){
		g=alphas[i];
		for(j=i+1;j<nalphas;j++){
			if(alphas[j]>g){
				g=alphas[j];
				alphas[j]=alphas[i];
				alphas[i]=g;
			}
		}
	}

	/* remove duplicates */
	for(i=1;i<nalphas;i++){
		if(alphas[i]!=alphas[i-1])
			alphas[n++]=alphas[i];
	}

	alphas=(double *)realloc(alphas, sizeof(double)*(n+1));
	alphas[n]=0; /* zero is always an alpha */
	nalphas=n;
}

/*===========================================================================
 copyOb | Copies the next component from 'obAlist' or 'obBlist' to 'objetA'
 or 'objetB'
 -----------------------------------------------------------------------------
 in  | 'A' for objetA or 'B' for objetB, the number of the component in the
 list
 ===========================================================================*/

void hof::HoF_Vector::copyOb(char object_name, int nobject)
{
	struct cell_objet *obList, *object;
	int size, i, j, *index;

	if(object_name=='A'){
		obList=obAlist;
		listeA1=NULL;
		listeA2=NULL;
		size=obAindex[nobject+1]-obAindex[nobject];
		j=obAindex[nobject];
		index=obAindex;
		object=objetA=(struct cell_objet *)calloc((size),sizeof(struct cell_objet));
		objetA_max=objetA+size;
	}
	else{
		obList=obBlist;
		listeB1=NULL;
		listeB2=NULL;
		size=obBindex[nobject+1]-obBindex[nobject];
		j=obBindex[nobject];
		index=obBindex;
		object=objetB=(struct cell_objet *)calloc((size),sizeof(struct cell_objet));
		objetB_max=objetB+size;
	}

	for(i=0;i<size;i++){
		object[i].x=obList[j+i].x;
		object[i].y=obList[j+i].y;
		object[i].pred=obList[j+i].pred-index[nobject];
		object[i].succ=obList[j+i].succ-index[nobject];
		object[i].fuzz=obList[j+i].fuzz;
	}

}


/*
 * DEBUG methods
 */
void hof::HoF_Vector::debug_xi(char nom_objet, int num_liste)
{
	char *s;
	int i, n;
	double *x;

	if(nom_objet=='A')
    {
		n=nb_trapA;
		if(num_liste==1) {x=xiA1; s=(char*)("xiA1");} else {x=xiA2; s=(char*)("xiA2");}
    }
	else
    {
		n=nb_trapB;
		if(num_liste==1) {x=xiB1; s=(char*)("xiB1");} else {x=xiB2; s=(char*)("xiB2");}
    }

	printf("%s   : ",s);
	for(i=0;i<n;i++) {printf("%3.2lf ",x[i]); if(x[i]<0.0) printf("x<0!!! ");}
	printf("\n");
}

void hof::HoF_Vector::debug_yi(int num_liste)
{
	int i;
	char *s;
	struct cell_yi *x;

	if(num_liste==1) {x=yiA1B1; s=(char*)("yiA1B1");}
	else {x=yiA2B2; s=(char*)("yiA2B2");}

	printf("%s : ",s);
	for(i=0;i<nb_trapA*nb_trapB;i++) printf("%3.2lf ",x[i].b1);
	printf("\n");
}

void hof::HoF_Vector::debug_cell_trap (struct cell_trap *liste)
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

void hof::HoF_Vector::debug_pc_tri()
{
	struct cell_objet* p=pc_tri;
	if(p)
	{
		printf("u=%f v=%f x=%f y=%f pred=%d succ=%d\n",
			   p->u, p->v, p->x, p->y, p->pred, p->succ);
		p=p->svt;
	}
	while(p && p!=pc_tri)
	{
		printf("u=%f v=%f x=%f y=%f pred=%d succ=%d\n",
			   p->u, p->v, p->x, p->y, p->pred, p->succ);
		p=p->svt;
	}
}
}
