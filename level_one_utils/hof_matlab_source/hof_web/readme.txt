This code implements the histograms of forces model described in:

P. Matsakis and L. Wendling, “A new way to represent the relative position between areal objects,” 
IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 27, pp. 634–643, 1999.

and used for linguistic description and pose estimation in:

Matsakis, P., Keller, J., Wendling, L., Marjamaa, J. and Sjahputera, O., "Linguistic Description of Relative Positions of Objects in Images", 
IEEE Transactions on Systems, Man, and Cybernetics, Vol. 31, No. 4, 2001, pp. 573-588.

Matsakis, P., Keller, J., Sjahputera, O., and Marjamaa, J. “The Use of Force Histograms for Affine-Invariant Relative Position Description”, 
IEEE Transactions on Pattern Analysis and Machine Intelligence, Vol. 26, No. 1, 2004, pp.1-18.

If you use this code and write a paper, please reference the above three manuscripts (the pdf versions available in IEEE Xplore.

Additionally, if you are using the code in an NGA application area, you might find the following papers of interest:

Sledge, I., Keller, J., Skubic, M., “Mapping Natural Language to Imagery: Placing Objects Intelligently”, 
Proceedings, IEEE International Conference on Fuzzy Systems, Jeju Island, Korea, August, 2009, pp. 518-524.

A. Buck, J. Keller, and M. Skubic, “A Modified Genetic Algorithm for Matching Building Sets with the Histograms of Forces”, 
to appear, Proceedings IEEE Congress on Evolutionary Computation, World Congress on Computational Intelligence, Barcelona, Spain,  July, 2010.


We have provided two ways in which you can use this code.  The first is through MATLAB MEX-files.  The second is through static C/C++ libraries.



**************************************************
MATLAB MEX-Files
**************************************************


**************************************************
Overview:
If you are unfamiliar with MATLAB MEX files, we recommend that you first visit the following:

http://www.mathworks.com/support/tech-notes/1600/1605.html

In short, MATLAB MEX-files are dynamically linked subroutines produced from C, C++ or Fortran source code that, when compiled, 
can be run from within MATLAB in the same way as MATLAB M-files or built-in functions.  
We have provided 3 sets of MEX-files for the histograms of forces, each of which have been compiled for a particular operating system and architecture:

32-bit Windows OS:
hof_vector.mexw32
hof_raster.mexw32

64-bit Windows OS:
hof_vector.mexw64
hof_raster.mexw64

Intel Mac OS X:
hof_vector.mexmaci
hof_raster.mexmaci

We have tested the 32- and 64-bit versions of these MEX-files on versions of Windows Vista and Windows 7 using MATLAB 2007a and MATLAB 2009b.  
They should, however, work under 32- and 64-bit versions of Windows, e.g., Windows XP, and with other versions of MATLAB.  
We have tested the Intel Mac OS X versions of the MEX-files on OS X 10.5 and 10.6 using MATLAB 2009b.  
Again, they should work under other versions of OS X and MATLAB.  

If you would like MEX-files for a different architecture or OS, please contact Professor Jim Keller at kellerj@missouri.edu and 
we will see if we can accomodate your request.


**************************************************
Installation and Use:
To use these MEX-files, place them in your MATLAB active directory or in a folder that is included to the active directory search path.  
For more information on how to do this, please consult: 

http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_env/br7ppws-1.html

There are two different functions you can call to compute the histograms of forces between two objects.  
The first is hof_vector(), which has the following function prototype:

h1 = hof_vector(objA,objB);

Here, objA and objB are, respectively, 2*n and 2*m arrays of object vertices, where 2 is the dimensionality of the objects, 
n is the number of vertices of object A and m is the number of vertices of object B.  
Each column in objA represents a new vertex in an object A, where the first row of objA are the x-coordinates and the second row are the y-coordinates.  
For example: objA = [3 3; 4 3; 4 4; 3 4]';.  The same is true for objB, i.e., each column in objB represents a new vertex in an object B, 
where the first row of objB are the x-coordinates and the second row are the y-coordinates.  
The vertices of objA must be sorted in either counter-clockwise or clockwise order and the vertices of objB must be sorted in the same manner.  
As such, you should not sort the vertices of objA in one order, e.g., counter-clockwise order, and the vertices of objB in the other one, e.g., clockwise order).  
Note that hof_vector() can handle convex and non-convex objects.  
However, hof_vector() does not currently support vector objects with holes in them, fuzzy vector objects, or disjoint objects.  
For this functionality, we recommend using hof_raster().

The second way to use hof_raster(), which has the following function prototype:

h1 = hof_raster(imgA,imgB);

Here, imgA and imgB are m*n single channel grayscale images with integer values ranging between 0 and 255.  
Although imgA and imgB only contain integer values, their data type should be double (and not uint8, short, etc.).  
Values of 0 in imgA denote that a particular pixel is not part of object A.  Values of 255 in imgA denote that a pixel perfectly belongs to object A.  
Similarly, values of 0 in imgB denote that a particular pixel is not part of object B and values of 255 in imgB denote that a pixel perfectly belongs to object B.  
Note that imgA and imgB may contain integer pixel values between 0 and 255.  These are used to represent fuzzy objects.  
Also note that hof_raster() can consider convex and non-convex objects which have holes in them and objects with disjoint components.

For both hof_vector and hof_raster, h1 will be a 2*181 array of the histogram of constant forces (first row of h1)
 and the histogram of gravitational forces (second row of h1).  
We have explicitly fixed the number of bins (directions) to be 181.  Note that in the event object A and object B intersect, 
only the histogram of constant forces will be defined.

We have provided some example test code and outputs with this file.



**************************************************
Static C++ Libraries
**************************************************

**************************************************
Overview:

If you are unfamiliar with static libraries, we commend that you consult the following:

http://en.wikipedia.org/wiki/Static_library
http://msdn.microsoft.com/en-us/library/ms235627%28VS.80%29.aspx

In short, a static library is a set of routines, external functions and variables which are resolved in a calling function at compile-time 
and copied into a target application by a compiler, linker, or binder, producing an object file and a stand-alone executable.  

Currently, we only provide a 32-bit static library for Windows.  We have tested this static library under 32- and 64-bit versions of Windows Vista and Windows 7.

If you would like a static library for a different architecture/OS, please contact Professor Jim Keller at kellerj@missouri.edu 
and we will see if we can accomodate your request.


**************************************************
Installation and Use:

To use the static libraries with applications built using Microsoft Visual Studio, you should place hof_vector.lib and hof_raster.lib in:
 C:\Program Files\Microsoft Visual Studio *\VC\lib.  
They should then be linked by going under Project -> Properties -> Configuration Properties -> Linker -> Input
 and adding hof_vector.lib and hof_raster.lib to the Additional Dependencies field.

There are two different functions you can call to compute the histograms of forces between two objects.  
The first is hof_vector(), which has the following function prototype:

h1 ** hof_vector(int ** objA, int ** objB, int numVertA, int numVertB);

Here, objA and objB are pointers to, respectively, numVertA*2 and numVertB*2 arrays of object vertices, where 2 is the dimensionality of the objects 
and numVertA is the number of vertices in object A and numVertB are the number of vertices of object B.  
Each row in objA corresponds to a new vertex of object A, where the first column denotes the x-coordinate of the vertex and the second column the y-coordinate.  
Similarly, each row in objB corresponds to a new vertex of object B, where the first column denotes the x-coordinate of the vertex and the second column the y-coordinate.  The vertices of objA must be sorted in either counter-clockwise or clockwise order and the vertices of objB must be sorted in the same manner.  As such, you should not sort the vertices of objA in one order, e.g., counter-clockwise order, and the vertices of objB in the other one, e.g., clockwise order).  Note that hof_vector() can handle convex and non-convex objects.  However, hof_vector() does not currently support vector objects with holes in them, fuzzy vector objects, or disjoint objects.  For this functionality, we recommend using hof_raster().  Also note that hof_vector() does NOT free up the memory used by either objA or objB, so you will need to do this yourself.

The second way to use hof_raster(), which has the following function prototype:

h1 ** hof_raster(int ** imgA, int ** imgB, int m, int n);

Here, imgA and imgB are pointers to m*n arrays with integer values ranging between 0 and 255.  Values of 0 in imgA denote that a particular pixel is not part of object A.  Values of 255 in imgA denote that a pixel perfectly belongs to object A.  Similarly, values of 0 in imgB denote that a particular pixel is not part of object B and values of 255 in imgB denote that a pixel perfectly belongs to object B.  Note that imgA and imgB may contain integer pixel values between 0 and 255.  These are used to represent fuzzy objects.  Also note that hof_raster() can consider convex and non-convex objects which have holes in them and objects with disjoint components.   As well, hof_raster() does NOT free up the memory used by either imgA or imgB, so you will need to do this yourself.

For both hof_vector and hof_raster, h1 will be a pointer to a 181*2 array of the histogram of constant forces (first column of h1) 
and the histogram of gravitational forces (second column of h1).  We have explicitly fixed the number of bins (directions) to be 181.  
Note that in the event object A and object B intersect, only the histogram of constant forces will be defined.

We have provided some example test code and outputs with this file.



**************************************************
Disclaimer
**************************************************
The code contained herein is to be taken "AS IS".  While we have tested the libraries, we do not guarantee that they are bug free.  If you discover a bug, 
please let Jim Keller know about it and we will try to provide a fix.
