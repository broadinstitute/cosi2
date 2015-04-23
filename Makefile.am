
ACLOCAL_AMFLAGS = -I m4
lib_LTLIBRARIES=
noinst_LTLIBRARIES=
bin_PROGRAMS=
EXTRA_DIST=ext/boost examples doc
TESTS=

BOOST_CPPFLAGS=-isystem $(srcdir)/ext/boost
AM_CPPFLAGS = -I$(srcdir)/.. $(BOOST_CPPFLAGS) -DNDEBUG 
AM_CXXFLAGS = $(CXX11_CXXFLAGS) -O3

include $(srcdir)/ext/boost/libs/make.am
include $(srcdir)/cosi_rand/Makefile.am
include $(srcdir)/cosi/Makefile.am
include $(srcdir)/recosim/Makefile.am
include $(srcdir)/tests/tests.am

#BOOST_PROGRAM_OPTIONS_LIBS = libprogopt.la
#lib_LTLIBRARIES += $(BOOST_PROGRAM_OPTIONS_LIBS) 
#BOOST_DIR=$(srcdir)/ext/boost

#PO=$(BOOST_DIR)/libs/program_options/src
#libprogopt_la_SOURCES = $(PO)/cmdline.cpp $(PO)/config_file.cpp $(PO)/convert.cpp $(PO)/options_description.cpp $(PO)/parsers.cpp $(PO)/positional_options.cpp $(PO)/split.cpp $(PO)/utf8_codecvt_facet.cpp $(PO)/value_semantic.cpp $(PO)/variables_map.cpp $(PO)/winmain.cpp



