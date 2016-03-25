
ACLOCAL_AMFLAGS = -I m4
lib_LTLIBRARIES=
noinst_LTLIBRARIES=
bin_PROGRAMS=
EXTRA_DIST=ext/boost examples doc release-info makedist.sh
TESTS=

AM_CPPFLAGS = -I$(srcdir) -isystem $(srcdir)/ext/boost -DNDEBUG 
AM_CXXFLAGS = $(CXX11_CXXFLAGS) -O3

include $(srcdir)/ext/boost/libs/make.am
include $(srcdir)/cosi/Makefile.am
include $(srcdir)/recosim/Makefile.am
include $(srcdir)/tests/tests.am

dist-hook:
	$(srcdir)/makedist.sh $(srcdir) $(distdir)
