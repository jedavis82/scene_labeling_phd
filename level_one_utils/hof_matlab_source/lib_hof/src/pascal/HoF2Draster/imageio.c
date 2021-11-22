/*
Implemented a number of functions for image operations.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "imageio.h"

/*
Create a empty Image (defined in imageio.h).
*/
IMAGE * createImage(int p_width, int p_height, int p_max_grey)
{
	int j;
	unsigned char ** pixels = NULL;
	IMAGE * image = (IMAGE *) calloc (1, sizeof(IMAGE));

	image->width = p_width;
	image->height = p_height;
	image->max_grey = p_max_grey;
	pixels = (unsigned char **) calloc (p_height, sizeof(unsigned char *));
	
	for(j = 0; j < p_height; j++)
		pixels[j] = (unsigned char *) calloc (p_width, sizeof(unsigned char));

	image->pixels = pixels;
	return image;
}

/*
Load an entire pgm image and save all the pixels into the created Image data structure.
*/
IMAGE * loadImage(char * p_filename)
{
	FILE* fp = NULL;
	IMAGE * image = NULL;
	unsigned char ** pixels = NULL;
    char line[256];
    int maxval, w, h;
    int binary;
    int i,j, int_tmp;

    if ((fp = fopen(p_filename, "rb")) == NULL)
		return NULL;

    fgets(line, 256, fp);
    if (strncmp(line,"P5", 2)) 
	{
	    if (strncmp(line,"P2", 2)) 
		{
		    fclose(fp);
			return NULL;
	    } else binary = 0;
    } else binary = 1;

    fgets(line, 256, fp);
    while (line[0] == '#') fgets(line, 256, fp);

    sscanf(line,"%d %d", &w, &h);
    fgets(line, 256, fp);
    sscanf(line, "%d", &maxval);

    image = createImage(w,h,maxval);

    pixels = image->pixels;

	if (binary)
	{
		for(j = 0; j < h; j++)
			fread((void*)(pixels[j]), sizeof(unsigned char), w, fp);
	}
    else
	      for (j = 0; j < h; j++) 
		  {
			  for(i = 0; i < w; i++)
			  {
					fscanf(fp,"%d", &int_tmp);
					pixels[j][i] = (unsigned char) int_tmp;
			  }
	      }


    fclose(fp);
    return(image);
}

/*
Only load the minimum bounding rectangle which contains the object (pixels grey values != p_background)
into a newly created Image data stucture.
*/

IMAGE * loadImageBound(char * p_filename, int p_background)
{
	int rows, cols, n_rows, n_cols, min_i, max_i, min_j, max_j, i, j;
	unsigned char ** org_pixels;
	unsigned char ** new_pixels;
	IMAGE * image = loadImage(p_filename);
	if(image == NULL)
		return NULL;
	rows = image->height;
	cols = image->width;
	org_pixels = image->pixels;
	min_i = rows - 1;
	max_i = 0;
	min_j = cols - 1;
	max_j = 0;
	for( i = 0; i < rows; i++)
	{
		for(j = 0; j < cols; j++)
		{
			if(org_pixels[i][j] != p_background)
			{
				if(i < min_i)
					min_i = i;
				if(i > max_i)
					max_i = i;
				if(j < min_j)
					min_j = j;
				if(j > max_j)
					max_j = j;
			}
		}
	}

	n_rows = max_i - min_i + 1;
	n_cols = max_j - min_j + 1;
	new_pixels = (unsigned char **)calloc(n_rows,sizeof(unsigned char*));
	for( i = 0; i < n_rows; i++)
	{
		new_pixels[i] = (unsigned char *) calloc (n_cols,sizeof(unsigned char));
		for(j = 0; j < n_cols; j++)
			new_pixels[i][j] = org_pixels[i + min_i][j + min_j];
	}
	for( i = 0; i < rows; i++)
		free(org_pixels[i]);
	free(org_pixels);
	image->height = n_rows;
	image->width = n_cols;
	image->pixels = new_pixels;
	return image;
}

/*
Copy an Image from p_image.
*/

IMAGE * copyImage(IMAGE * p_image)
{
	int w,h;
	int i,j;
	IMAGE * tar_image = NULL;
	unsigned char ** org_pixels = NULL;
	unsigned char ** tar_pixels = NULL;

	w = p_image->width;
	h = p_image->height;
	org_pixels = p_image->pixels;

	tar_image = createImage(w,h,p_image->max_grey);
	
	tar_pixels = tar_image->pixels;

	for (j = 0; j < h; j++) 
	{
		for (i = 0; i < w; i++)
		{
			tar_pixels[j][i] = org_pixels[j][i];
		}
	}
	return tar_image;
}

/*
Write a Image into a pgm file for presenting.
*/

int writeImage(char * p_filename, IMAGE * p_image)
{
	FILE* fp = NULL;
	unsigned char ** pixels = NULL;
	int i,j;
	int w,h;
	
	if ((fp = fopen(p_filename, "w")) == NULL) 
		return 0;

	w = p_image->width;
	h = p_image->height;
	pixels = p_image->pixels;

	fprintf( fp, "P2\n" );
	fprintf( fp, "%d %d\n", w, h );
	fprintf( fp, "%d\n", p_image->max_grey);

	for (j = 0; j < h; j++) 
	{
		for (i = 0; i < w; i++)
		{
			fprintf(fp,"%d ",(int)(pixels[j][i]));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	return 1;
}

/*
Release the memory of p_image.
*/
void destroyImage(IMAGE * p_image)
{
	int i, h;
	unsigned char ** pixels = NULL;
	h = p_image->height;
	pixels = p_image->pixels;
	for(i = 0; i < h; i++)
		free(pixels[i]);
	free(pixels);
	free(p_image);
}

/*
Create an Image containing two concentric objects, one ball (with radius r1)
surrounded by a shell (with inner radius r2 and outer radius r3).
*/

IMAGE * ball_shell(int w, int h, int r1, int r2, int r3, unsigned char b, unsigned char s)
{
	int i,j;
	int o_i,o_j;
	double r;
	double r1_2, r2_2, r3_2;
	IMAGE * image = createImage(w, h, 255);
	r1_2 = pow(r1,2);
	r2_2 = pow(r2,2);
	r3_2 = pow(r3,2);
	
	o_i = (h+1)/2;
	o_j = (w+1)/2;
	for(i=0;i<h;i++)
	{
		for(j=0;j<w;j++)
		{
			r = pow(i-o_i,2) + pow(j-o_j,2);
			if(r <= r1_2)
				image->pixels[i][j] = b;
			else if(r >= r2_2 && r <= r3_2)
				image->pixels[i][j] = s;
			else
				image->pixels[i][j] = 255;
		}
	}
	return image;
}

/*Create an Image containing two concentric objects, one shell (with inner radius r1 and outer radius r2)
surrounded by another shell (with inner radius r3 and outer radius r4).*/

IMAGE * shell_shell(int w, int h, int r1, int r2, int r3, int r4, unsigned char b, unsigned char s)
{
	int i,j;
	int o_i,o_j;
	double r;
	double r1_2, r2_2, r3_2, r4_2;
	IMAGE * image = createImage(w, h, 255);
	r1_2 = pow(r1,2);
	r2_2 = pow(r2,2);
	r3_2 = pow(r3,2);
	r4_2 = pow(r4,2);
	
	o_i = (h+1)/2;
	o_j = (w+1)/2;
	for(i=0;i<h;i++)
	{
		for(j=0;j<w;j++)
		{
			r = pow(i-o_i,2) + pow(j-o_j,2);
			if(r >= r1_2 && r <= r2_2)
				image->pixels[i][j] = b;
			else if(r>=r3_2 && r<=r4_2)
				image->pixels[i][j] = s;
			else
				image->pixels[i][j] = 255;
		}
	}
	return image;
}
