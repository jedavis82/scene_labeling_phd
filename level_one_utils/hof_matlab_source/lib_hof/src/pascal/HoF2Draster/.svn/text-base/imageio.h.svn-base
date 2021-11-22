#ifndef IMAGEIO_H
#define IMAGEIO_H

/*
Data structure for storing a pgm (raw) image.
*/

typedef struct image
{
	int width;                   /*image width*/
	int height;                  /*image height*/
	int max_grey;                /*image mas grey level, usually the value is 255*/
	unsigned char ** pixels;     /*2D array for all the pixel grey values*/
} IMAGE;

/*
Link list for storing all the vertexes of a vector region.
*/
typedef struct vertex
{
	int y;                        /*vertical position of the vertex*/
	int x;                        /*horizontal position of the vertex*/
	struct vertex * next;         /*link list*/
} VERTEX;


IMAGE * createImage(int p_width, int p_height, int p_max_grey);
IMAGE * loadImage(char * p_filename);
IMAGE * loadImageBound(char * p_filename, int p_background);
IMAGE * copyImage(IMAGE * p_image);
int writeImage(char * p_filename, IMAGE * p_image);
void destroyImage(IMAGE * p_image);
IMAGE * ball_shell(int w, int h, int r1, int r2, int r3, unsigned char b, unsigned char s);
IMAGE * shell_shell(int w, int h, int r1, int r2, int r3, int r4, unsigned char b, unsigned char s);


#endif
