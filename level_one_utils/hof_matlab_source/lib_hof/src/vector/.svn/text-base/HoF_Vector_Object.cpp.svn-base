#include "HoF_Vector_Object.hpp"

void hof::HoF_Vector_Object::readObjectFile( const std::string& filename )
{
	std::ifstream ifs;

	try
	{
		ifs.open(filename.c_str(), std::ifstream::in);
	}
	catch(std::exception &e)
	{
		std::stringstream ss;
		ss << "[hof::HoF_Vector_Object::readObjectFile] failed to open filename=" << filename;
		std::cout << "ERROR " << ss.str() << std::endl;
		throw std::runtime_error(ss.str());
	}

	/*
	 * Read how many alpha layers in this object.
	 */
	unsigned int alpha_count;
	ifs >> alpha_count;

	_alpha.clear();

	/*
	 * Scan alpha values.
	 */
	for(unsigned int alpha_idx=0; alpha_idx < alpha_count; alpha_idx++)
	{
		double alpha;
		ifs >> alpha;
		_alpha.push_back(alpha);
	}

	if ( alpha_count != _alpha.size() )
	{
		std::stringstream ss;
		ss << "[hof::HoF_Vector_Object::readObjectFile]  while reading filename=" << filename
				<< " alpha_count=" << alpha_count << " does not match with _alpha.size()=" << _alpha.size();
		std::cout << "ERROR " << ss.str() << std::endl;
		throw std::runtime_error( ss.str() );
	}

	/*
	 * Scan the polygon vertices for each alpha layer.
	 */

	_object.clear();

	for(unsigned int alpha_idx=0; alpha_idx < alpha_count; alpha_idx++)
	{
		/*
		 * Read total number of vertices from all polygons in this alpha layer.
		 */

		unsigned int total_vertex_count;
		ifs >> total_vertex_count;

		std::vector<std::vector<std::vector<std::pair<double,double> > > > alpha_layer_polygons;

		std::vector<std::vector<std::pair<double,double> > > polygon_group_vertices;

		unsigned int vertex_counter = 0;

		while( vertex_counter < total_vertex_count )
		{
			unsigned int polygon_vertex_counter;
			ifs >> polygon_vertex_counter;

			vertex_counter += polygon_vertex_counter;

			std::vector<std::pair<double,double> > polygon_vertices;

			for( unsigned int vertex_idx=0; vertex_idx < polygon_vertex_counter; vertex_idx++)
			{
				double x, y;
				ifs >> x >> y;
				polygon_vertices.push_back( std::pair<double,double>(x,y) );
			}

			/*
			 * Check if the new polygon is contained within the first polygon in the current
			 * polygon group.
			 */

			/*
			 * polygon_group is empty, insert current polygon into group
			 * no need to check for containment.
			 */
			if ( polygon_group_vertices.empty() )
				polygon_group_vertices.push_back( polygon_vertices );
			else
			{
				/*
				 * Current polygon is contained within at least one of the polygons
				 * in the polygon group. Add current polygon to the group.
				 */
				if ( polygonIsContainedInPolygonGroup(polygon_vertices, polygon_group_vertices ) )
					polygon_group_vertices.push_back( polygon_vertices );
				else
				{
					/*
					 * Current polygon is not contained within any of the polygons
					 * in the polygon group. Save polygon group. And start a new
					 * polygon group with the current polygon as its first member.
					 */
					alpha_layer_polygons.push_back( polygon_group_vertices );
					polygon_group_vertices.clear();
					polygon_group_vertices.push_back(polygon_vertices);
				}
			}
		}

		alpha_layer_polygons.push_back( polygon_group_vertices );

		_object[_alpha[alpha_idx]] = alpha_layer_polygons;
	}
	ifs.close();
	setName( filename );
}

/***************************************************************************************
 *
 */

void hof::HoF_Vector_Object::readFileObjectStruct( hof::file_object* file_obj )
{
	_alpha.clear();

	setName( std::string(file_obj->name) );

	/*
	 * Transfer the number of alpha and the alpha values.
	 */

	unsigned int alpha_count = static_cast<unsigned int>( file_obj->nalphas );

	for( unsigned int i=0; i<alpha_count; i++ )
		_alpha.push_back( file_obj->alphas[i] );

	/*
	 * File_Object.numob indicates the number of polygon group in the object across all alpha cuts.
	 * File_Object.index length is numob+1. It stores indices of vertex that belong
	 * to each polygon group. There can be multiple polygons in a polygon group.
	 * But there is no separation for each polygon. One way to detect:
	 * The beginning polygon vertices (predecessor index is larger than current):
	 * File_Object.ob[i].pred > i
	 *
	 * The end of polygon vertices (successor index is smaller than current):
	 * File_Object.ob[i].succ < i
	 *
	 * Total number of elements in array File_Object.ob can be found in:
	 *
	 * File_Object.index[ File_Object.numob ]
	 *
	 */

	_object.clear();

	std::vector<std::vector<std::pair<double,double> > > polygon_group_vertices;

	double curr_alpha = 0;

	/*
	 * Scan vertices in each polygon group.
	 */

	for(int i=0; i < file_obj->numob; i++ )
	{
		int start_idx = file_obj->index[i];
		int end_idx   = file_obj->index[i+1];

		if ( file_obj->ob[start_idx].fuzz == file_obj->ob[end_idx-1].fuzz )
		{
			curr_alpha = file_obj->ob[start_idx].fuzz;
		}
		else
		{
			std::stringstream ss;
			ss << "[hof::HoF_Vector_Object::readFileObjectStruct] alpha value mismatch for ob[" << start_idx
			   << "]=" << file_obj->ob[start_idx].fuzz << " and ob[" << end_idx << "]=" << file_obj->ob[end_idx].fuzz;
			std::cout << "ERROR " << ss.str() << std::endl;
			throw std::runtime_error( ss.str() );
		}

		/*
		 * A polygon group could contain multiple simple polygon.
		 */
		std::vector<std::pair<double,double> > polygon_vertices;

		/*
		 * Scan the vertices, watch for the end of a simple polygon.
		 */
		for( int j=start_idx; j<end_idx; j++)
		{
			polygon_vertices.push_back( std::pair<double,double>(file_obj->ob[j].x, file_obj->ob[j].y) );

			/*
			 * This is the end of a simple polygon. Push it to the polygon group and
			 * reset the polygon vertices.
			 */
			if ( file_obj->ob[j].succ < j )
			{
				polygon_group_vertices.push_back( polygon_vertices );
				polygon_vertices.clear();
			}
		}

		/*
		 * Add this polygon group to an alpha layer.
		 */
		if ( _object.count(curr_alpha) > 0 )
		{
			std::vector< std::vector< std::vector< std::pair<double,double> > > > new_alpha_layer;
			new_alpha_layer.push_back( polygon_group_vertices );
			_object[curr_alpha] = new_alpha_layer;
		}
		else
			_object[curr_alpha].push_back( polygon_group_vertices );

		polygon_group_vertices.clear();
	}
}

/***************************************************************************************
 *
 */

hof::file_object* hof::HoF_Vector_Object::createFileObjectStruct()
{
	hof::file_object *fob = (struct hof::file_object *) malloc(sizeof(struct hof::file_object));

	fob->name = (char *)malloc(sizeof(char)*(_name.size()+1));
	strcpy(fob->name, _name.c_str());

	fob->nalphas = static_cast<int>( _alpha.size() );
	fob->alphas = (double *)malloc(sizeof(double)*_alpha.size());

	for(int i=0; i<fob->nalphas; ++i)
		fob->alphas[i] = _alpha.at(i);

	unsigned int total_vertices = getTotalVertices();
	fob->ob = (struct cell_objet *)malloc(sizeof(struct cell_objet)*total_vertices);

	unsigned int total_polygon_group = getTotalPolygonGroup();
	fob->numob = total_polygon_group;
	fob->index = (int *)malloc(sizeof(int)*(total_polygon_group+1));
	fob->index[0] = 0;
	/*
	 * Polygon group index.
	 */
	unsigned int g = 1;
	/*
	 * Vertices index.
	 */
	unsigned int v = 0;

	/*
	 * For each alpha
	 */
	for(unsigned int i=0; i<_alpha.size(); i++)
	{
		double alpha = _alpha[i];
		/*
		 * For each polygon group in alpha, calculate the number of vertices.
		 */

		for(unsigned int j=0; j<_object[alpha].size(); j++)
		{
			int vertices_counter = 0;
			/*
			 * For each simple polygon in the polygon group.
			 */
			for(unsigned int k=0; k<_object[alpha].at(j).size(); k++)
			{
				/*
				 * Get the number of vertices in this simple polygon.
				 */
				vertices_counter += static_cast<int>(_object[alpha].at(j).at(k).size());
				/*
				 * For each vertex in this simple polygon, transfer (x,y) coordinate.
				 */
				int v_start = v;
				for(unsigned int l=0; l<_object[alpha].at(j).at(k).size(); l++, v++)
				{
					fob->ob[v].x = _object[alpha].at(j).at(k).at(l).first;
					fob->ob[v].y = _object[alpha].at(j).at(k).at(l).second;
					fob->ob[v].fuzz = alpha;
					fob->ob[v].pred = v-1;
					fob->ob[v].succ = v+1;
				}

				fob->ob[v_start].pred = v-1;
				fob->ob[v-1].succ = v_start;
			}
			fob->index[g] = fob->index[g-1] + vertices_counter;
		}
	}

	return fob;
}

/***************************************************************************************
 *
 */

unsigned int hof::HoF_Vector_Object::addAlphaLayer( const double& alpha, const std::vector<std::vector< std::vector < std::pair<double,double> > > >& polygon_group_vertices )
{
	_object.erase( _object.find(alpha) );
	_object[alpha] = polygon_group_vertices;

	/*
	 * Insert new alpha layer to the ordered list
	 * Alpha cuts are ordered in increasing order of alpha.
	 */

	unsigned int inserted_idx = 0;
	unsigned int initial_size = static_cast<unsigned int>(_alpha.size());

	if ( _alpha.front() > alpha )
	{
		/*
		 * Insert new alpha to front.
		 */
		_alpha.insert( _alpha.begin(), alpha );
		return 0;
	}
	else
	if ( _alpha.back() < alpha )
	{
		/*
		 * Push new alpha to back;
		 */
		_alpha.push_back( alpha );
		return ( static_cast<unsigned int> ( _alpha.size() - 1 ) );
	}
	else
	{
		/*
		 * Find the first element in _alpha and insert new alpha before it.
		 */
		for(unsigned int i=0; i<_alpha.size(); ++i)
			if ( _alpha.at(i) > alpha )
			{
				std::vector<double>::iterator it = _alpha.begin();
				_alpha.insert( it+i, alpha );
				inserted_idx = i;
				break;
			}

		if ( inserted_idx == 0 && initial_size == static_cast<unsigned int>(_alpha.size()) )
		{
			std::stringstream ss;
			ss << "<hof::HoF_Vector_Object::addAlphaLayer> alpha does not get inserted to _alpha";
			std::cout << "ERROR " << ss.str() << std::endl;
			throw std::runtime_error(ss.str());
		}
	}
	return inserted_idx;
}

/*
 *
 */
void hof::HoF_Vector_Object::deleteAlphaLayer( const double& alpha )
{
	if ( _object.count(alpha) > 0 )
	{
		_object.erase( _object.find(alpha) );

		for( unsigned int i=0; i<_alpha.size(); ++i)
			if ( _alpha.at(i) == alpha )
				_alpha.erase( _alpha.begin() + i );
	}
}

/*
 *
 */

unsigned int hof::HoF_Vector_Object::addPolygonGroup( const double& alpha, std::vector<std::vector<std::pair<double,double> > > polygon_group )
{
	if ( _object.count(alpha) > 0 )
	{
		_object[alpha].push_back( polygon_group );
		return ( static_cast<unsigned int> ( _object[alpha].size() - 1 ) );
	}
	else
	{
		std::vector<std::vector<std::vector<std::pair<double,double> > > > alpha_layer;
		alpha_layer.push_back( polygon_group );
		addAlphaLayer( alpha, alpha_layer );
		return (0);
	}
}

/*
 * Returns true if a simple polygon is contained within any simple polygon
 * that is a member of a polygon group.
 */

bool hof::HoF_Vector_Object::polygonIsContainedInPolygonGroup ( const std::vector<std::pair<double,double> >& polygon_contained,
														 const std::vector<std::vector<std::pair<double,double> > >& polygon_group )
{
	for( unsigned int i=0; i<polygon_group.size(); i++)
	{
		if ( polygonIsContainedInPolygon( polygon_contained, polygon_group[i]))
			return true;
	}
	return false;
}

/*
 * Returns true if a simple polygon is contained inside another simple polygon.
 *
 */

bool hof::HoF_Vector_Object::polygonIsContainedInPolygon ( const std::vector<std::pair<double,double> >& polygon_contained, const std::vector<std::pair<double,double> >& polygon_container )
{

	for( unsigned int i=0; i<polygon_contained.size(); i++)
	{
		if ( !vertexIsContainedInPolygon( polygon_contained[i].first, polygon_contained[i].second, polygon_container ))
			return false;
	}
	return true;
}

/*
 * Returns true if a vertex (x,y) is contained within a simple polygon.
 */
bool hof::HoF_Vector_Object::vertexIsContainedInPolygon( const double& x, const double& y, const std::vector<std::pair<double,double> >& polygon )
{
	unsigned int i, pred_i, succ_i;
	bool up=false, down=false, right=false, left=false;
	double segv;

	double curr_x, curr_y;
	double succ_x, succ_y;
	double pred_x, pred_y;

	/* for each edge in component */
	for(i=0;i<polygon.size();i++)
	{
		pred_i = i > 0 ? i-1 : polygon.size()-1;
		succ_i = i < (polygon.size()-1) ? i+1 : 0;

		curr_x = polygon[i].first;
		curr_y = polygon[i].second;

		pred_x = polygon[pred_i].first;
		pred_y = polygon[pred_i].second;

		succ_x = polygon[succ_i].first;
		succ_y = polygon[succ_i].second;

		if((!up || !down) &&
		   (((x+hof::SEUIL_ZERO)>=curr_x && x<=(succ_x+hof::SEUIL_ZERO)) ||
			(x<=(curr_x+hof::SEUIL_ZERO) && (x+hof::SEUIL_ZERO)>=succ_x)))
		{
			/* coordinate must be under or over edge */
			if(fabs(succ_y-curr_y)<hof::SEUIL_ZERO)
			{
				/* edge is horizontal */
			    if(curr_y>y)
			    	up=true;
			    else
			    	down=true;
			}
			else if(fabs(succ_x-curr_x)>hof::SEUIL_ZERO)
			{
				/* check for intersect */
				segv=(succ_y-curr_y)/(succ_x-curr_x)*(x-curr_x)+curr_y;
				if(segv>y)
					up=true;
				else
					down=true;
			 }
		}

		if((!right || !left) &&
		   (((y+hof::SEUIL_ZERO)>=curr_y && y<=(succ_y+hof::SEUIL_ZERO)) ||
			(y<=(curr_y+hof::SEUIL_ZERO) && (y+hof::SEUIL_ZERO)>=succ_y))){ /* coordinate mus be to the left or right of edge */
			   if(fabs(succ_x-curr_x)<hof::SEUIL_ZERO){ /* edge is vertical */
				   if(curr_x>x)
					   right=true;
				   else
					   left=true;
			   }
			   else if(fabs(succ_y-curr_y)>hof::SEUIL_ZERO){ /* check for intersect */
				   segv=(y-curr_y)/((succ_y-curr_y)/(succ_x-curr_x))+curr_x;
				   if(segv>x)
					   right=true;
				   else
					   left=true;
			   }
		   }

	}
	if(up && down && right && left)
		return(true); /* must have four intersects */
	else
		return(false);
}

/*
 *
 */

unsigned int hof::HoF_Vector_Object::getTotalVertices()
{
	unsigned int counter = 0;

	std::map<double, std::vector< std::vector< std::vector< std::pair<double,double> > > > >::iterator it;

	/*
	 * For each alpha layer
	 */
	for( it=_object.begin(); it!=_object.end(); ++it )
	{
		/*
		 * For each polygon group in an alpha layer
		 */
		for( unsigned int i=0; i<it->second.size(); ++i )
		{
			/*
			 * For each simple polygon in a polygon group
			 */
			for( unsigned int j=0; j<it->second.at(i).size(); ++j )
			{
				counter += it->second.at(i).at(j).size();
			}
		}
	}
	return ( counter );
}

/*
 *
 */

unsigned int hof::HoF_Vector_Object::getTotalPolygonGroup()
{
	unsigned int counter = 0;

	std::map<double, std::vector< std::vector< std::vector< std::pair<double,double> > > > >::iterator it;

	/*
	 * For each alpha layer
	 */
	for( it=_object.begin(); it!=_object.end(); ++it )
	{
		counter += it->second.size();
	}
	return ( counter );
}

/*
 * Print to stdout.
 */
void  hof::HoF_Vector_Object::print()
{
	/*
	 * How many alpha we have.
	 */
	std::stringstream ss;
	ss << "alpha_count " << _alpha.size();
	std::cout << ss.str() << std::endl;
	ss.str("");
	ss << "alpha";
	for(unsigned int i=0; i<_alpha.size();++i)
		ss << " " << _alpha[i];
	std::cout << ss.str() << std::endl;

	/*
	 * For each alpha layers.
	 */
	for(unsigned int a=0; a<_alpha.size();++a)
	{
		std::cout << "polygon_group_count " << _object[_alpha.at(a)].size() << std::endl;

		/*
		 * For each polygon group in the alpha layer.
		 */

		for(unsigned int g=0; g<_object[_alpha.at(a)].size(); ++g)
		{
			std::cout << g << " polygon_count " << _object[_alpha.at(a)].at(g).size() << std::endl;

			/*
			 * For each simple polygon in the polygon group.
			 */

			for(unsigned int i=0; i<_object[_alpha.at(a)].at(g).size(); ++i)
			{
				std::cout << g << " " << i << " vertex_count " << _object[_alpha.at(a)].at(g).at(i).size() << std::endl;

				ss.str("");

				for(unsigned int j=0; j<_object[_alpha.at(a)].at(g).at(i).size(); ++j)
						ss << _object[_alpha.at(a)].at(g).at(i).at(j).first << " " << _object[_alpha.at(a)].at(g).at(i).at(j).second << " ";

				std::cout << ss.str() << std::endl;
			}
		}
	}
}

/*
 *
 */

void hof::HoF_Vector_Object::write(std::string outpath)
{
	;
}
