ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(srcdir)/..
AM_LDFLAGS = -L$(srcdir)/../cosi_rand -lm

bin_PROGRAMS = recosimulate get_recomap recomap_hapmap2

LIBRAND = $(abs_top_builddir)/../cosi_rand/libcosi_rand.la

recosimulate_SOURCES = recosim.c
recosimulate_LDADD = $(LIBRAND)

get_recomap_SOURCES = get_recomap.c
get_recomap_LDADD = $(LIBRAND)

recomap_hapmap2_SOURCES = recomap_hapmap2.c
recomap_hapmap2_LDADD = $(LIBRAND)





