AC_PREREQ(2.63)

AC_INIT([lib_hof], [1.0.0], [sjahputerao@missouri.edu])
AC_CONFIG_AUX_DIR(config)

AM_INIT_AUTOMAKE(lib_hof,1.0.0)

AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_PATH_PROG(CXXTEST, [cxxtestgen.pl])
if test -z "$CXXTEST"; then AC_MSG_ERROR([cxxtest not found]); fi

AC_LANG([C++])

AC_CHECK_HEADER([cxxtest/TestSuite.h], [], AC_MSG_ERROR([cxxtest/TestSuite.h not found]))
AC_CHECK_HEADER([gdal_priv.h], [], AC_MSG_ERROR([gdal_priv.h not found]))
AC_CHECK_HEADER([cpl_string.h], [], AC_MSG_ERROR([cpl_string.h not found]))
AC_CHECK_LIB([gdal], [main], [], AC_MSG_ERROR([libgdal not found]))
AM_CXXFLAGS="-Wall -Wextra -rdynamic"
AC_SUBST([AM_CXXFLAGS])

AC_CONFIG_FILES([
                 Makefile
                 doc/Makefile
                 src/Makefile
                 src/raster/Makefile
                 src/vector/Makefile
                 tests/Makefile
                 ])

AC_OUTPUT
