ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS=-I$(srcdir)/.. -O3 -DNDEBUG
AM_CXXFLAGS = -O3 -std=c++0x --param max-inline-insns-single=1000 -Wall -Wno-unused -Werror
AM_LDFLAGS = -L$(srcdir)/../cosi_rand
AM_DEFAULT_SOURCE_EXT = .cc

lib_LTLIBRARIES = libmiscutil.la libmiscutil_pg.la libmiscutil_valg.la libmiscutil_valgnoinl.la
libmiscutil_la_SOURCES = utils.cc

subpkgincludedir = $(includedir)/miscutil
subpkginclude_HEADERS = utils.h typedval.h defs.h


libmiscutil_pg_la_SOURCES = $(libmiscutil_la_SOURCES)
libmiscutil_pg_la_CFLAGS = $(AM_CFLAGS) -pg

libmiscutil_valg_la_SOURCES = $(libmiscutil_la_SOURCES)
libmiscutil_valg_la_CFLAGS = $(AM_CFLAGS) -ggdb3

libmiscutil_valgnoinl_la_SOURCES = $(libmiscutil_la_SOURCES)
libmiscutil_valgnoinl_la_CFLAGS = $(AM_CFLAGS) -ggdb3 -fno-inline





