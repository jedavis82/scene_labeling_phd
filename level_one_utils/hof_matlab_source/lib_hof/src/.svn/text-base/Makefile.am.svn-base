EXAMPLE_LIBRARY_VERSION=1:0:0

SUBDIRS = raster vector

lib_LTLIBRARIES = lib_hof.la

lib_hof_la_LDFLAGS = -version-info $(EXAMPLE_LIBRARY_VERSION)

lib_hof_la_SOURCES = 

lib_hof_la_LIBADD = 	raster/lib_hof_raster.la \
						vector/lib_hof_vector.la 

tldir = $(includedir)/cgi/hof
tl_HEADERS = hof.hpp \
			 constants.hpp

rasterdir = $(includedir)/cgi/hof/raster
raster_HEADERS = 	raster/raster.hpp \
					raster/HoF_Raster.hpp

vectordir = $(includedir)/cgi/hof/vector
vector_HEADERS = 	vector/vector.hpp \
					vector/HoF_Vector.hpp \
					vector/HoF_Vector_Structs.hpp \
					vector/HoF_Vector_Object.hpp
