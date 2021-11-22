/*
 *  HoF_Vector.hpp
 *
 *  Created on: Dec 13, 2010
 *      Author: Ozy_2
 *
 *  C++ Wrapper for Pascal's HoF code for vector data ver3.0 June 2010.
 */

#ifndef HOF_VECTOR_HPP_
#define HOF_VECTOR_HPP_

#include "HoF_Vector_Structs.hpp"
#include "HoF_Vector_Object.hpp"
#include "../constants.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cstring>

namespace hof{

	/*
	=======================================================================================
	 ABOUT THE OBJECTS
	 ---------------------------------------------------------------------------------------

	 This module handles 2D objects in vector form. An object may have holes and multiple
	 connected components. It may also be fuzzy. It is assumed that the number of member-
	 ship degrees is finite. It is also assumed that any alpha-cut with 0<alpha<=1 can
	 be expressed---using the union and difference set operations---in terms of a finite
	 number of distinct simple polygons; these polygons must be such that an edge of a
	 polygon does not intersect an edge of another polygon. Each object is described in
	 a text file, as shown below. The file contains integer and floating point values
	 separated by space or line feed characters (' ' or '\n'). The object is described
	 as a set of alpha-cuts (sorted by increasing alpha), each alpha-cut is described
	 as a set of polygons (in any order), each polygon as a set of vertices (listed
	 either clockwise or counterclockwise), and each vertex as a pair of coordinates
	 x (from left to right) and y (from bottom to top).

	 ------------
	 First line: number (int) of distinct nonzero membership degrees,
				 followed by these membership degrees (double) in ascending order
	 Next line:  total number (int) of the vertices that define the 1st alpha-cut
	 Next line:  number of vertices (int) of the 1st simple polygon
				 that defines the 1st alpha-cut
	 Next line:  coordinates x and y (double) of the 1st vertex of the 1st simple polygon
				 that defines the 1st alpha-cut, coordinates of the 2nd vertex, etc.
	 Next line:  number of vertices (int) of the 2nd simple polygon
				 that defines the 1st alpha-cut
	 Next line:  coordinates x and y (double) of the 1st vertex of the 2nd simple polygon
				 that defines the 1st alpha-cut, coordinates of the 2nd vertex, etc.
	 ......
	 Next line:  total number (int) of the vertices that define the 2nd alpha-cut
	 Next line:  number of vertices (int) of the 1st simple polygon
				 that defines the 2nd alpha-cut
	 Next line:  coordinates x and y (double) of the 1st vertex of the 1st simple polygon
				 that defines the 2nd alpha-cut, coordinates of the 2nd vertex, etc.
	 Next line:  number of vertices (int) of the 2nd simple polygon
				 that defines the 2nd alpha-cut
	 Next line:  coordinates x and y (double) of the 1st vertex of the 2nd simple polygon
				 that defines the 2nd alpha-cut, coordinates of the 2nd vertex, etc.
	 ......
	 ------------

	 Note that any pair (A,B) of objects and any type r of force can be considered.
	 Be aware, however, that infinite forces are obtained when r>=1 and the interior
	 of A intersects the interior of B, or when r>=2 and A is adjacent to B (A and B
	 intersect, but the interior of A does not intersect the interior of B).

	 */

/* Trying to collect all these HoF functions in a class */

class HoF_Vector{
public:
	HoF_Vector()
	{
		/* initialize class members */

		numA = 0;
		numB = 0;
		obAlist = NULL;
		obAindex = NULL;
		obBindex = NULL;
		alphas = NULL;
		nalphas=0;
		type_force=0.0;
		objetA=NULL;
		objetB=NULL;
		listeA1=NULL;
		listeB1=NULL;
		listeA2=NULL;
		listeB2=NULL;
		nb_trapA=0;
		nb_trapB=0;
		xiA1=NULL;
		xiB1=NULL;
		xiA2=NULL;
		xiB2=NULL;
		yiA1B1=NULL;
		yiA2B2=NULL;
	}

	~HoF_Vector()
	{
		;
	}

	/*
	 * ----------------------------------------------------------------------------------------
	 *	int computeHistogram (double *histogram, int scheme, int numberOfDirections,
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
	 * ------------------------------------------------------------------------------------------
	 */
	int computeHistogram (	double *histogram, int scheme,
					    	int numberOfDirections, double typeOfForce,
					    	struct file_object *argument, struct file_object *referent );

	int computeHistogram( 	std::vector<double>& histogram, int scheme, int numberOfDirection, double typeOfForce,
							HoF_Vector_Object& argument, HoF_Vector_Object& referent );

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
	int writeHistogram(char *nameHistogramFile, double *histogram, int numberOfDirections, char *nameObjectA, char *nameObjectB, char *op);

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
	int readObjectFiles(struct file_object **file, char** filenames, int nfiles);

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
	struct file_object *getObject(struct file_object **files, int nfiles, char *filename);

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
	void freeObjects(struct file_object **files, int n);

	void freeObject(struct file_object *);

private:
	/*
	 * ------------------------------------------------------------------
	 * 	void setTypeOfForce (double r)
		 ---------------
		 Sets the type of force r.
		 ---------------
		 in     | Type of force (e.g., 0.0 for constant forces, 2.0
				| for gravitational forces). It can be any real number.
	 * ------------------------------------------------------------------
	 */
	void setArgAndRef(struct file_object *arg, struct file_object *ref);

	/*
	 * 	 void setTypeOfForce (double r)
		 ---------------
		 Sets the type of force r.
		 ---------------
		 in     | Type of force (e.g., 0.0 for constant forces, 2.0
				| for gravitational forces). It can be any real number.
	 *
	 */
	void setTypeOfForce (double r);

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
	double computeForce_SingleScheme (double theta);

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
	double computeForce_DoubleScheme (double theta);

	double trier_sommets(double theta, double *umin, double *umax);
	void init_trapezes(double ucrt);
	double maj_trapezes();
	void inserer_sommets(double ucrt, int *p_multA, int *p_multB);
	void detruire(struct cell_trap *liste);
	cell_trap* copie(struct cell_trap *liste1);

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
	double force_trapezes(double u1, double u2, int *invalid);

	void init_xi(char nom_objet, int num_liste);

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
	void init_yi(int num_liste);

	void copyOb(char object_name, int object);

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
	double overlapping_force_trapezes(double u1, double u2, struct cell_yi *yi1, struct cell_yi *yi2);

	/*===========================================================================
	 isContained | Determines if a coordinate is contained within a component by
	 checking for vertical and horizontal intersects.
	 -----------------------------------------------------------------------------
	 in  | The coordinate values, the number of the component, the component index
	 array, and the object list.
	 out | 1 if the coordinate is contained in the component, 0 otherwise.
	 ===========================================================================*/
	int isContained (double x, double y, int num, int *index, struct cell_objet *list);

	/*===========================================================================
	 sortAlphas | Sorts the list of alphas and removes duplicates.
	 ===========================================================================*/
	void sortAlphas();

	void debug_xi(char nom_objet, int num_liste);
	void debug_yi(int num_liste);
	void debug_cell_trap (struct cell_trap *liste);
	void debug_pc_tri();

	void freeFileObject( file_object* fob );

private:
	int numA;
	int numB;
	struct cell_objet *obAlist;
	struct cell_objet *obBlist;
	int *obAindex;
	int *obBindex;

	double *alphas;
	int nalphas;

	double type_force;

	struct cell_objet *objetA;
	struct cell_objet *objetA_max;
	struct cell_objet *objetB;
	struct cell_objet *objetB_max;
	struct cell_objet *pc_tri;

	struct cell_trap *listeA1;
	struct cell_trap *listeB1;
	struct cell_trap *listeA2;
	struct cell_trap *listeB2;

	int nb_trapA;
	int nb_trapB;

	double *xiA1;
	double *xiB1;
	double *xiA2;
	double *xiB2;

	struct cell_yi *yiA1B1;
	struct cell_yi *yiA2B2;


};
#endif /* HOF_VECTOR_HPP_ */
