

bin_PROGRAMS += recosimulate get_recomap recomap_hapmap2

LIBRAND = $(builddir)/libcosi_rand.la

recosimulate_SOURCES = $(srcdir)/recosim/recosim.c
recosimulate_LDADD = $(LIBRAND)
recosimulate_LDFLAGS = -lm

get_recomap_SOURCES = $(srcdir)/recosim/get_recomap.c
get_recomap_LDADD = $(LIBRAND)
get_recomap_LDFLAGS = -lm

recomap_hapmap2_SOURCES = $(srcdir)/recosim/recomap_hapmap2.c
recomap_hapmap2_LDADD = $(LIBRAND)
recomap_hapmap2_LDFLAGS = -lm






