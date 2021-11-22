#ifndef HOF_VECTOR_OBJECT_HPP_
#define HOF_VECTOR_OBJECT_HPP_

#include "HoF_Vector_Structs.hpp"
#include "../constants.hpp"

#include <string>
#include <vector>
#include <map>
#include <utility>
#include <fstream>
#include <sstream>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <cstring>
#include <cstdlib>
#include <cmath>

namespace hof{

class HoF_Vector_Object
{
public:
	HoF_Vector_Object()
	{
		;
	}

	~HoF_Vector_Object()
	{

	}

	/*
	 * Performs similar task to HoF_Vector::readObjectFile().
	 * Reading and loading object information from file following the
	 * format specified in HoF_Vector_Structs.hpp which is identical to
	 * Pascal's FHistogramVectorV3.0.c
	 */
	void readObjectFile( const std::string& filename );

	/*
	 * Given a file_objet struct that holds a definition of an object,
	 * copy all information to this class.
	 */
	void readFileObjectStruct( hof::file_object* file_obj );

	/*
	 * Create a copy of object information in the forms of hof::File_Object struct.
	 * This allows us to use Pascal's HoF methods that expect object definitions
	 * given in hof::File_Object struct.
	 */
	hof::file_object* createFileObjectStruct();

	/*
	 * User can add one alpha layer at a time. All polygons and their vertices must be provided.
	 */
	unsigned int addAlphaLayer( const double& alpha, const std::vector< std::vector< std::vector < std::pair<double,double> > > >& polygon_vertices );

	/*
	 *
	 */
	void deleteAlphaLayer( const double& alpha );

	/*
	 *
	 */
	unsigned int addPolygonGroup( const double& alpha, std::vector<std::vector<std::pair<double,double> > > polygon_group );

	/*
	 * User can add one polygon to a specific alpha layer. If an alpha layer is not
	 * available, a new layer will be created. Otherwise, the polygon will be added to
	 * the existing list of polygons in the specified alpha layer.
	 */
	void addPolygonGroupToAlphaLayer ( const double& alpha, const std::vector<std::vector< std::pair<double,double> > >& polygon );

	/*
	 * Get all polygons in an alpha layer.
	 */
	std::vector<std::vector< std::vector< std::pair<double,double> > > > getAlphaLayer( const double* alpha );

	/*
	 * Check if an alpha layer exists. Returns 0 if it does not exist.
	 * Else, return the number of polygon if it does.
	 */

	bool checkAlpha( const double& alpha )
	{
		if ( _object.count( alpha ) > 0 )
			return static_cast<unsigned int>( _object[alpha].size() );
		else
			return 0;
	}

	/*
	 * Get alpha count
	 */

	unsigned int getAlphaCount()
	{
		return static_cast<unsigned int> (_object.size());
	}

	/*
	 * Set object name.
	 */
	void setName ( const std::string& name ) { _name = name; }

	/*
	 * Get object name.
	 */
	std::string getName() { return _name; }

	/*
	 *
	 */
	unsigned int getTotalVertices();

	/*
	 *
	 */
	unsigned int getTotalPolygonGroup();

	/*
	 * Print vector object definition to stdout
	 */
	void print();

	/*
	 * Write vector object definition to text file.
	 */
	void write(std::string outpath);

private:
	bool polygonIsContainedInPolygonGroup ( const std::vector<std::pair<double,double> >& polygon_contained,
											const std::vector<std::vector<std::pair<double,double> > >& polygon_group );
	bool polygonIsContainedInPolygon ( const std::vector<std::pair<double,double> >& polygon_contained, const std::vector<std::pair<double,double> >& polygon_container );
	bool vertexIsContainedInPolygon( const double& x, const double& y, const std::vector<std::pair<double,double> >& polygon );

private:
	std::string _name;

	/*
	 * object[alpha].at(g).at(i).at(j) = (x,y)
	 * is the coordinate of the j-th vertex in the i-th polygon belonging to the
	 * g-th polygon group of alpha layer.
	 *
	 * A polygon group is used to describe a complex polygon that may contain
	 * smaller polygons. This is used to represents an object with holes, that may contain
	 * another polygon inside the holes. See D.vdata for example and read the file
	 * src/pascal/HoF2Dvector/examples.pdf.
	 */

	std::map<double, std::vector< std::vector< std::vector< std::pair<double,double> > > > > _object;

	std::vector<double> _alpha;
};


}
#endif
