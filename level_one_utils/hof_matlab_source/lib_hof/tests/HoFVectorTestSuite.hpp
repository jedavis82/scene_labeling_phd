#ifndef HOFVECTORTESTSUITE_H_
#define HOFVECTORTESTSUITE_H_

#include <cxxtest/TestSuite.h>

#include "../src/vector/vector.hpp"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <exception>
#include <gdal_priv.h>
#include <cpl_string.h>

class SetupFixture : public CxxTest::GlobalFixture
{
	public:
		bool setUpWorld()
		{
			std::cout << "HoFRasterTestSetup" << std::endl;
			tearDownWorld();
			GDALAllRegister();
			return true;
		}

		bool tearDownWorld() { return true; }

		bool setUp() { return true; }

		bool tearDown() { return true; }
};

static SetupFixture loggingFixture;

class HoFRasterTestSuite : public CxxTest::TestSuite
{
	public:

		void setUp()
		{
			init_data();
		}

		void tearDown()
		{
		}

		std::vector<unsigned char> copy_object( std::vector<unsigned char>& canvas, unsigned int canvas_size_x, unsigned int canvas_size_y,
							std::vector<unsigned char>& object, unsigned int object_size_x, unsigned int object_size_y,
							unsigned int x_offset, unsigned int y_offset )
		{
			std::vector<unsigned char> target( canvas );

			for( unsigned int y=0; y<object_size_y; ++y)
				for( unsigned int x=0; x<object_size_x; ++x)
					target[((y_offset + y) * canvas_size_x) + (x_offset + x )] = object[(y*object_size_x) + x];

			return(target);
		}

		void init_data()
		{
			/* Create a canvas 256x256 pixels */

			_canvas_size_x = 256;
			_canvas_size_y = 256;
			_object_size_x = 32;
			_object_size_y = 32;
			_p0 = 0.01;
			_p1 = 3.0;

			_black_canvas.resize( _canvas_size_x*_canvas_size_y, 0 );

			/* Square object 32x32 pixels */

			_square.resize( _object_size_x*_object_size_y, 255 );

			_obj_A = copy_object(_black_canvas, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128, 128);

			_obj_B_right = copy_object(_black_canvas, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128+64, 128);

			_obj_B_right_above = copy_object(_black_canvas, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128+64, 128-96);

			_obj_B_above = copy_object(_black_canvas, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128, 128-96);

			_obj_B_left_above = copy_object(_black_canvas, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128-96, 128-96);

			_obj_B_left = copy_object(_black_canvas, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128-96, 128);

			_obj_B_left_below = copy_object(_black_canvas, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128-96, 128+64);

			_obj_B_below = copy_object(_black_canvas, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128, 128+64);

			_obj_B_right_below = copy_object(_black_canvas, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128+64, 128+64);

			_obj_B_adjacent = copy_object(_black_canvas, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128+32, 128);

			_obj_B_half_intersect = copy_object(_black_canvas, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128+16, 128);

			_obj_B_surround = copy_object(_black_canvas, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128+64, 128);
			_obj_B_surround = copy_object(_obj_B_surround, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128+64, 128-32);
			_obj_B_surround = copy_object(_obj_B_surround, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128+64, 128-64);
			_obj_B_surround = copy_object(_obj_B_surround, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128+32, 128-64);
			_obj_B_surround = copy_object(_obj_B_surround, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128, 128-64);
			_obj_B_surround = copy_object(_obj_B_surround, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128-32, 128-64);
			_obj_B_surround = copy_object(_obj_B_surround, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128-64, 128-64);
			_obj_B_surround = copy_object(_obj_B_surround, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128-64, 128-32);
			_obj_B_surround = copy_object(_obj_B_surround, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128-64, 128);
			_obj_B_surround = copy_object(_obj_B_surround, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128-64, 128+32);
			_obj_B_surround = copy_object(_obj_B_surround, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128-64, 128+64);
			_obj_B_surround = copy_object(_obj_B_surround, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128-32, 128+64);
			_obj_B_surround = copy_object(_obj_B_surround, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128, 128+64);
			_obj_B_surround = copy_object(_obj_B_surround, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128+32, 128+64);
			_obj_B_surround = copy_object(_obj_B_surround, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128+64, 128+64);
			_obj_B_surround = copy_object(_obj_B_surround, _canvas_size_x, _canvas_size_y,_square, _object_size_x, _object_size_y, 128+64, 128+32);

			_numberDirections = 32;

			std::cout << "_numberDirections = " << _numberDirections << std::endl;

			std::cout << "canvas_size = " << _canvas_size_x << "x" << _canvas_size_y << std::endl;

			std::cout << "object_size = " << _object_size_x << "x" << _object_size_y << std::endl;

			std::cout << "p0 = " << _p0 << " p1=" << _p1 << std::endl;
		}

		/*
		 * Print histogram to screen;
		 */
		std::string printHistogram( const std::vector<double>& h )
		{
			std::stringstream ss;

			unsigned int i = 0;
			for(; i<h.size()-1; ++i)
				ss << h[i] << ",";
			ss << h[i];

			return ss.str();
		}


		std::vector<unsigned char> UnionImages( const std::vector<unsigned char> imgA, const std::vector<unsigned char> imgB )
		{
			std::vector<unsigned char> U(imgA);

			for(unsigned int i=0; i<U.size(); ++i)
			{
				if ( U[i] == 0  && imgB[i] > 0 )
					U[i] = imgB[i];
			}

			return U;
		}



		void writeImage( std::vector<unsigned char> imgA, std::vector<unsigned char> imgB, int size_x, int size_y, std::string filename )
		{
			std::string drivername("GTiff");
			GDALDriver *driver = NULL;
			driver = GetGDALDriverManager()->GetDriverByName(drivername.c_str());
			if (driver == NULL)
			{
				std::cout << "Failed to get GTiff driver for creating file '" << filename << "'" << std::endl;
				return;
			}
			char **options = NULL;
			GDALDataset *ds = (GDALDataset*) driver->Create( filename.c_str(), size_x, size_y, 1, GDT_Byte, options);
			if ( ds == NULL )
			{
				std::cout << "Failed to create image '" << filename << "'" << std::endl;
			}
			else
			{
				std::vector<unsigned char> U = UnionImages(imgA,imgB);
				ds->RasterIO(GF_Write, 0,0, size_x, size_y, &U[0], size_x, size_y, GDT_Byte, 1, NULL, 0,0,0);
				ds->FlushCache();
				delete ds;
			}
		}

		///
		///
		///

		void test_Raster_Disjoint_F0()
		{
			hof::HoF_Raster H;

			std::cout << "test_Raster_Disjoint_F0()" << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 0, &_obj_A[0], &_obj_B_right[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F0 Right: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 0, &_obj_A[0], &_obj_B_right_above[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F0 Right-Above: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 0, &_obj_A[0], &_obj_B_above[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F0 Above: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 0, &_obj_A[0], &_obj_B_left_above[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F0 Left-Above: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 0, &_obj_A[0], &_obj_B_left[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F0 Left: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 0, &_obj_A[0], &_obj_B_left_below[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F0 Left-Below: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 0, &_obj_A[0], &_obj_B_below[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F0 Below: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 0, &_obj_A[0], &_obj_B_right_below[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F2 Right-Below: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 0, &_obj_A[0], &_obj_B_adjacent[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F0 Adjacent-Right: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 0, &_obj_A[0], &_obj_B_half_intersect[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F0 Half-Intersect-Right: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 0, &_obj_A[0], &_obj_B_surround[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F0 Surround: " << printHistogram(_histogram) << std::endl;
		}

		void test_Raster_Disjoint_F2()
		{
			hof::HoF_Raster H;

			std::cout << "test_Raster_Disjoint_F2()" << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 2.0, &_obj_A[0], &_obj_B_right[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F2 Right: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 2.0, &_obj_A[0], &_obj_B_right_above[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F2 Right-Above: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 2.0, &_obj_A[0], &_obj_B_above[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F2 Above: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 2.0, &_obj_A[0], &_obj_B_left_above[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F2 Left-Above: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 2.0, &_obj_A[0], &_obj_B_left[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F2 Left: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 2.0, &_obj_A[0], &_obj_B_left_below[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F2 Left-Below: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 2.0, &_obj_A[0], &_obj_B_below[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F2 Below: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 2.0, &_obj_A[0], &_obj_B_right_below[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F2 Right-Below: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 2.0, &_obj_A[0], &_obj_B_adjacent[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F2 Adjacent-Right: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 2.0, &_obj_A[0], &_obj_B_half_intersect[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F2 Half-Intersect-Right: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.FRHistogram_CrispRaster (&_histogram[0], _numberDirections, 2.0, &_obj_A[0], &_obj_B_surround[0], _canvas_size_x, _canvas_size_y);
			std::cout << "F2 Surround: " << printHistogram(_histogram) << std::endl;
		}

		void test_Raster_Disjoint_FH()
		{
			hof::HoF_Raster H;

			std::cout << "test_Raster_Disjoint_FH()" << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.F02Histogram_CrispRaster (&_histogram[0], _numberDirections, &_obj_A[0], &_obj_B_right[0], _canvas_size_x, _canvas_size_y, _p0, _p1);
			std::cout << "F02 Right: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.F02Histogram_CrispRaster (&_histogram[0], _numberDirections, &_obj_A[0], &_obj_B_right_above[0], _canvas_size_x, _canvas_size_y, _p0, _p1);
			std::cout << "F02 Right-Above: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.F02Histogram_CrispRaster (&_histogram[0], _numberDirections, &_obj_A[0], &_obj_B_above[0], _canvas_size_x, _canvas_size_y, _p0, _p1);
			std::cout << "F02 Above: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.F02Histogram_CrispRaster (&_histogram[0], _numberDirections, &_obj_A[0], &_obj_B_left_above[0], _canvas_size_x, _canvas_size_y, _p0, _p1);
			std::cout << "F02 Left-Above: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.F02Histogram_CrispRaster (&_histogram[0], _numberDirections, &_obj_A[0], &_obj_B_left[0], _canvas_size_x, _canvas_size_y, _p0, _p1);
			std::cout << "F02 Left: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.F02Histogram_CrispRaster (&_histogram[0], _numberDirections, &_obj_A[0], &_obj_B_left_below[0], _canvas_size_x, _canvas_size_y, _p0, _p1);
			std::cout << "F02 Left-Below: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.F02Histogram_CrispRaster (&_histogram[0], _numberDirections, &_obj_A[0], &_obj_B_below[0], _canvas_size_x, _canvas_size_y, _p0, _p1);
			std::cout << "F02 Below: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.F02Histogram_CrispRaster (&_histogram[0], _numberDirections, &_obj_A[0], &_obj_B_right_below[0], _canvas_size_x, _canvas_size_y, _p0, _p1);
			std::cout << "F02 Right-Below: " << printHistogram(_histogram) << std::endl;

		}

		void test_Raster_Overlap_FH()
		{
			hof::HoF_Raster H;

			std::cout << "test_Raster_Overlap_FH()" << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.F02Histogram_CrispRaster (&_histogram[0], _numberDirections, &_obj_A[0], &_obj_B_adjacent[0], _canvas_size_x, _canvas_size_y, _p0, _p1);
			std::cout << "F02 Adjacent-Right: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.F02Histogram_CrispRaster (&_histogram[0], _numberDirections, &_obj_A[0], &_obj_B_half_intersect[0], _canvas_size_x, _canvas_size_y, _p0, _p1);
			std::cout << "F02 Half-Intersect-Right: " << printHistogram(_histogram) << std::endl;

			_histogram.resize( static_cast<unsigned int>(_numberDirections+1), 0 );
			H.F02Histogram_CrispRaster (&_histogram[0], _numberDirections, &_obj_A[0], &_obj_B_surround[0], _canvas_size_x, _canvas_size_y, _p0, _p1);
			std::cout << "F02 Surround: " << printHistogram(_histogram) << std::endl;
		}

		void test_write_image()
		{
			writeImage( _obj_A, _obj_B_right,_canvas_size_x, _canvas_size_y, "right.tif" );
			writeImage( _obj_A, _obj_B_right_above, _canvas_size_x, _canvas_size_y, "right_above.tif" );
			writeImage( _obj_A, _obj_B_above, _canvas_size_x, _canvas_size_y, "above.tif" );
			writeImage( _obj_A, _obj_B_left_above, _canvas_size_x, _canvas_size_y, "left_above.tif" );
			writeImage( _obj_A, _obj_B_left, _canvas_size_x, _canvas_size_y, "left.tif" );
			writeImage( _obj_A, _obj_B_left_below, _canvas_size_x, _canvas_size_y, "left_below.tif" );
			writeImage( _obj_A, _obj_B_below, _canvas_size_x, _canvas_size_y, "below.tif" );
			writeImage( _obj_A, _obj_B_right_below, _canvas_size_x, _canvas_size_y, "right_below.tif" );
			writeImage( _obj_A, _obj_B_adjacent, _canvas_size_x, _canvas_size_y, "adjacent.tif" );
			writeImage( _obj_A, _obj_B_half_intersect, _canvas_size_x, _canvas_size_y, "half_intersect.tif" );
			writeImage( _obj_A, _obj_B_surround, _canvas_size_x, _canvas_size_y, "surround.tif" );
		}

	private:
		std::vector<double> _histogram;
		int _numberDirections;
		int _canvas_size_x;
		int _canvas_size_y;
		int _object_size_x;
		int _object_size_y;
		double _p0;
		double _p1;
		std::vector<unsigned char> _square;
		std::vector<unsigned char> _black_canvas;
		std::vector<unsigned char> _obj_A;
		std::vector<unsigned char> _obj_B_right;
		std::vector<unsigned char> _obj_B_right_above;
		std::vector<unsigned char> _obj_B_above;
		std::vector<unsigned char> _obj_B_right_below;
		std::vector<unsigned char> _obj_B_below;
		std::vector<unsigned char> _obj_B_left_below;
		std::vector<unsigned char> _obj_B_left;
		std::vector<unsigned char> _obj_B_left_above;
		std::vector<unsigned char> _obj_B_adjacent;
		std::vector<unsigned char> _obj_B_half_intersect;
		std::vector<unsigned char> _obj_B_surround;
};


#endif /* TESTESTSUITE_H_*/
