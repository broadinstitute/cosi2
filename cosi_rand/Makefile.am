ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS=-I$(srcdir)/.. -DNDEBUG

lib_LTLIBRARIES = libcosi_rand.la
libcosi_rand_la_SOURCES = expdev.c gamma.c gammln.c mtwist.c multinom.c poisson.c ranbinom.c rng.c

noinst_LTLIBRARIES = libcosi_rand_noinst.la
libcosi_rand_noinst_la_SOURCES = $(libcosi_rand_la_SOURCES)

subpkgincludedir = $(includedir)/cosi_rand
subpkginclude_HEADERS = expdev.h gamma.h gammln.h mtwist.h multinom.h poisson.h ranbinom.h random.h





