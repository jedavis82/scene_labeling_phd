#ifndef HOF_VECTOR_STRUCTS_HPP_
#define HOF_VECTOR_STRUCTS_HPP_


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

/*
 *
 */

namespace hof{


#ifndef cell_objet
struct cell_objet
{
	double x, y;
	int pred, succ;
	double fuzz;
	double u, v;
	struct cell_objet *pdt, *svt;
};
#endif

/*
 *
 */
#ifndef cell_trap
struct cell_trap
{
	struct cell_trap *svt;
	double v;
	int est_sommet;
	int point;
	int mult;
	double a, b;
};
#endif
/*
 *
 */
#ifndef cell_yi
struct cell_yi
{
	double a1, a2;
	double b1, b2;
	double c1, c2;
	double d1, d2;
	double x1, x2;
	double z1, z2;
};
#endif
/*
 *
 */
#ifndef file_object
struct file_object
{
	char *name;
	struct cell_objet *ob;
	int numob;
	int *index;
	int nalphas;
	double *alphas;
};
#endif

}

#endif // HOF_VECTOR_STRUCTS_HPP_
