#ifndef OBJFUNC_H
#define OBJFUNC_H

#include "imageio.h"

#define Y 0
#define X 1

/*This is the data structure used by angle histogram*/
/*
typedef struct a_obj
{
	int n;
	double * value;
	int * x;
	int * y;
} A_OBJ;
*/


/*
A raster object loaded from a pgm file
*/
typedef struct ras_obj
{
	int rows;                /*number of rows of the host image*/
	int cols;                /*number of columns of the host image*/
	int org_y;               /*the vertical position of the origin*/
	int org_x;               /*the horizontal position of the origin*/
	double ** re;            /*2D array storing the pixels' membership degrees (belong to the range [0,1])*/
	double ** im;            /*Aid 2D array for fourrier transformation*/
} RAS_OBJ;

/*
A 1D vector.
*/
typedef struct comp_vec
{
	int length;              /*The length of the vector*/
	double * re;             /*1D array recording the real part of each complex value*/
	double * im;             /*1D array recording the imaginary part of the compelx value*/
} COMP_VEC;

/*
A minimum bounding rectangle of a raster object.
*/
typedef struct bound_box
{
	int v1[2];               /* up-left corner of the rectangle */
	int v2[2];               /* low-left corner of the rectangle */
	int v3[2];               /* low-right corner of the rectangle */
	int v4[2];               /* up-right corner of the rectangle */
} BOUND_BOX;

/*
A link list recording all possible alpha-cut values extracted from a fuzzy raster object.
*/

typedef struct alpha
{
	double f;                /*One possible membership degrees existing in the raster object*/
	struct alpha * next;     /*Link list*/
}ALPHA;

/*
Link list of alpha cut values
*/

typedef struct a_list
{
	int n;                   /*Number of different alpha cut values in the list*/
	ALPHA * head;            /*The header of the list*/
}A_LIST;

/*
A_OBJ * loadAngleObject(RAS_OBJ * p_obj);
A_OBJ * image2AngleObject(IMAGE * image);
void destroyAngleObject(A_OBJ * p_obj);
*/

RAS_OBJ * createRasObj(int p_rows, int p_cols, int p_org_y, int p_org_x);
RAS_OBJ * createRasObjComplex(int p_rows, int p_cols, int p_org_y, int p_org_x);
double ** createMatrix(int p_rows, int p_cols);

RAS_OBJ * image2RasObj(IMAGE * p_image);
RAS_OBJ * image2RasObjExtract(IMAGE * p_image, int p_mark, int p_mark_joint);

RAS_OBJ * copyRasObj(RAS_OBJ * p_object, int p_only_re);
int copyMatrix(double ** p_org, double ** p_tar, int p_rows, int p_cols);
int extendRasObj(RAS_OBJ * p_object);
RAS_OBJ * resizeRasObj(RAS_OBJ * p_object, int p_t_rows, int p_t_cols, int p_org_y, int p_org_x, int region_y, int region_x, int bo_rows, int bo_cols);
RAS_OBJ * resizeRasObj_boundBox(RAS_OBJ * p_object, BOUND_BOX * p_box, int p_t_rows, int p_t_cols, int p_org_y, int p_org_x);
RAS_OBJ * resizeRasObj_boundBox_region(RAS_OBJ * p_object, BOUND_BOX * p_box, int p_t_rows, int p_t_cols, int p_org_y, int p_org_x, 
									   int bound_rows, int bound_cols);

double ** extendMatrix(double ** p_org, int p_rows, int p_cols, int p_t_rows, int p_t_cols, int p_org_y, int p_org_x);

int reflectRasObj(RAS_OBJ * p_object);
int reflectRasObj_boundBox(RAS_OBJ * p_object, BOUND_BOX * p_box);

int reflectMatrix(double ** p_matrix, int p_rows, int p_cols, int p_org_y, int p_org_x);
int complexRasObj(RAS_OBJ * p_object);
int magnitudeRasObj(RAS_OBJ * p_object);
int normalizRasObj(RAS_OBJ * p_object);
IMAGE * rasObj2Image(RAS_OBJ * p_object, int p_maxval);
void destroyRasObj(RAS_OBJ * p_object);
void destroyMatrix(double ** p_matrix, int p_rows);

COMP_VEC * createCompVec(int p_length);
COMP_VEC * copyCompeVec(COMP_VEC * p_vector);
int copyVec(double * p_org, double * p_tar, int p_length);
int vec2Row(double * p_vector, int p_row, double ** p_matrix, int start, int end);
int vec2Col(double * p_vector, int p_col, double ** p_matrix, int start, int end);
int row2Vec(double ** p_matrix, int p_row, double * p_vector, int start, int end);
int col2Vec(double ** p_matrix, int p_col, double * p_vector, int start, int end);
void destroyCompVec(COMP_VEC * p_vector);

BOUND_BOX * loadBoundBox(RAS_OBJ * p_object);
int boundBox(BOUND_BOX * p_box, RAS_OBJ * p_object);

int randomFuzzify(RAS_OBJ * p_object, int seed, int p_number);

A_LIST * createAList();
int insertAList(A_LIST * p_list, double p_f);
void printfAList(A_LIST * p_list);
void destroyAList(A_LIST * p_list);
int detectAlphas(RAS_OBJ * p_object, A_LIST * p_list);
int getAlphaCut(RAS_OBJ * p_object, RAS_OBJ * p_cut, double p_f);

int histogramWrite(FILE * p_fp, double * p_histogram, int p_d);

void swap(int * a, int *b);
void bresenham(unsigned char ** pixels, VERTEX * start, VERTEX * end, unsigned char color);
VERTEX * init_koch_snowflake(int n);
void iteration_koch_snowflake(int iter, VERTEX * list);
void draw_koch_snowflake(int rows, int cols, unsigned char ** pixels, VERTEX * list, unsigned b_color, unsigned char f_color);
int compare_image(int rows, int cols, unsigned char ** p1, unsigned char ** p2);

#endif
