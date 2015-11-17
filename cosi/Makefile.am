#
# Makefile: cosi makefile
#
# Describes how to build the cosi library, libcosi, and the main cosi program, 'coalescent'.
#
# Run "make doc" to build HTML documentation.
#
# Most of cosi code is compiled into a library, libcosi, with which the main
# executable is then linked.
#

#SUBDIRS=libsigc++-2.3.1

COSI_CPPFLAGS = $(AM_CPPFLAGS) -DCOSI_TRACK_LAST_COAL #-I$(builddir)/libsigc++-2.3.1 -I$(srcdir)/libsigc++-2.3.1 
COSI_CXXFLAGS= $(AM_CXXFLAGS)

#AM_LDFLAGS = -L$(srcdir)/../cosi_rand
AM_DEFAULT_SOURCE_EXT = .cc

################################################
#
# Section: Building cosi
#
################################################

lib_LTLIBRARIES += libcosiutil.la libcosi.la 
bin_PROGRAMS += coalescent

if DEVOPTS
COSI_CXXFLAGS += --param max-inline-insns-single=1000 \
	-Werror -Wall -Wextra -Wno-unused-function -Wuninitialized \
	-Wredundant-decls -Wpointer-arith -Wcast-qual -Wcast-align -Wformat \
	-Woverloaded-virtual -Wstrict-null-sentinel -Wlogical-op \
	-Wmissing-declarations -Wdouble-promotion -Winit-self

endif  # if DEVOPTS

#if COALAPX
# Turn on support for approximating the coalescent by restricting coalescence to nodes whose
# genetic material overlaps or nearly-overlaps.
COSI_CPPFLAGS += -DCOSI_SUPPORT_COALAPX
#endif

#if CONDSNP
# Turn on support for conditioning the generated simulations on those that have a SNP at a
# specified location with specified sample frequencies.
#COSI_CPPFLAGS += -DCOSI_CONDSNP
#endif

#if EHH
# Turn on support for computing EHH and related statistics.
#COSI_CPPFLAGS += -DCOSI_EHH -DCOSI_LEAFSET_BITSET
#endif

SD=$(srcdir)/cosi

COSIUTILHEADERS = $(SD)/utils.h $(SD)/arrival.h $(SD)/generalmath.h $(SD)/cosirand.h $(SD)/gauss_legendre.h \
	$(SD)/valgrind.h $(SD)/evtgrp.h $(SD)/typedval.h \
	$(SD)/order_statistics.hpp $(SD)/random_access_skip_list.h $(SD)/skip_list_detail.h $(SD)/utildefs.h \
	$(SD)/drd.h $(SD)/memcheck.h $(SD)/mempool.h 
libcosiutil_la_SOURCES = $(SD)/utils.cc $(SD)/arrival.cc $(SD)/cosirand.cc $(SD)/gauss_legendre.c $(SD)/mempool.cc \
	$(COSIUTILHEADERS)
libcosiutil_la_CPPFLAGS = $(AM_CPPFLAGS)


COSISOURCES = $(SD)/coalesce.cc $(SD)/coalescent.cc $(SD)/demography.cc $(SD)/file.cc $(SD)/geneconversion.cc \
	$(SD)/historical.cc \
	$(SD)/migrate.cc $(SD)/mutate.cc $(SD)/mutlist.cc $(SD)/node.cc $(SD)/output.cc $(SD)/pop.cc $(SD)/recomb.cc \
	$(SD)/seglist.cc $(SD)/segsumm.cc $(SD)/simulator.cc \
	$(SD)/sweep.cc $(SD)/sweep1.cc $(SD)/sweep2.cc $(SD)/sweep3.cc $(SD)/genmap.cc $(SD)/stats.cc $(SD)/hooks.cc \
	$(SD)/traj.cc $(SD)/cositop.cc $(SD)/mutcontext.cc \
	$(SD)/condsnp.cc $(SD)/hullmgr.cc $(SD)/leafset-tree.cc $(SD)/leafset-counts.cc

libcosi_la_SOURCES = $(SD)/coalesce.cc $(SD)/coalescent.cc $(SD)/demography.cc $(SD)/file.cc $(SD)/geneconversion.cc \
	$(SD)/historical.cc \
	$(SD)/migrate.cc $(SD)/mutate.cc $(SD)/mutlist.cc $(SD)/node.cc $(SD)/output.cc $(SD)/pop.cc $(SD)/recomb.cc $(SD)/seglist.cc $(SD)/segsumm.cc $(SD)/simulator.cc \
	$(SD)/sweep.cc $(SD)/sweep1.cc $(SD)/sweep2.cc $(SD)/sweep3.cc $(SD)/genmap.cc $(SD)/stats.cc $(SD)/hooks.cc $(SD)/traj.cc $(SD)/cositop.cc \
	$(SD)/callgrind.h $(SD)/defs.h $(SD)/coalesce.h $(SD)/demography.h $(SD)/ehh.h $(SD)/ehh.cc $(SD)/file.h $(SD)/geneconversion.h $(SD)/historical.h $(SD)/leafset.h \
	$(SD)/leafset-bitset.h $(SD)/leafset-tree.h $(SD)/leafset-sizeonly.h $(SD)/leafset-sizeonly.cc $(SD)/leafset-counts.h $(SD)/mutcontext.cc \
	$(SD)/migrate.h $(SD)/mutate.h $(SD)/node.h $(SD)/nodefwd.h $(SD)/nodelist.h $(SD)/output.h $(SD)/pop.h $(SD)/popfwd.h $(SD)/recomb.h $(SD)/seglist.h $(SD)/seglistfwd.h \
	$(SD)/seglistfwd2.h $(SD)/segsumm.h $(SD)/simulator.h $(SD)/cositop.h \
	$(SD)/sweep.h $(SD)/sweep1.h $(SD)/sweep2.h $(SD)/sweep3.h $(SD)/recomb.h $(SD)/genmap.h $(SD)/stats.h $(SD)/traj.h $(SD)/decls.h $(SD)/hooks.h \
	$(SD)/coalrate.h $(SD)/condsnp.h $(SD)/condsnp.cc \
	$(SD)/mutcontext.h $(SD)/hullmgr.h $(SD)/hullmgr.cc \
	$(SD)/cosicfg.h $(SD)/coalescent.h $(SD)/mutlist.h $(SD)/output_coal_data_cosi.cc $(SD)/coal_data_cosi.h $(SD)/customstats.h \
	$(SD)/calc_pop_stats_cosi.cc $(SD)/calc_Fst_cosi.cc \
	$(COSIUTILHEADERS)

libcosi_la_CPPFLAGS = $(COSI_CPPFLAGS)
libcosi_la_CXXFLAGS = $(COSI_CXXFLAGS)
libcosi_la_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
libcosi_la_LIBADD = $(builddir)/libcosi_rand.la $(builddir)/libcosiutil.la $(BOOST_PROGRAM_OPTIONS_LIBS) \
	$(BOOST_FILESYSTEM_LIBS) $(BOOST_TIMER_LIBS) #$(builddir)/libsigc++-2.3.1/sigc++/libsigc-2.0.la

coalescent_SOURCES = $(SD)/cosimain.cc
coalescent_CPPFLAGS = $(COSI_CPPFLAGS)
coalescent_CXXFLAGS = $(COSI_CXXFLAGS)
coalescent_LDADD = $(builddir)/libcosi.la $(builddir)/libcosi_rand.la $(BOOST_PROGRAM_OPTIONS_LIBS)	\
	$(BOOST_FILESYSTEM_LIBS) $(BOOST_TIMER_LIBS) $(BOOST_SERIALIZATION_LIBS) #$(builddir)/libsigc++-2.3.1/sigc++/libsigc-2.0.la
coalescent_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_FILESYSTEM_LDFLAGS)

#if EHH
#libcosi_la_SOURCES += leafset-bitset.cc
#else
libcosi_la_SOURCES += $(SD)/leafset-tree.cc $(SD)/leafset-counts.cc
#endif

#COSI_CXX_FLAGS_DEBUG += -O0 -UNDEBUG -D_GLIBCXX_DEBUG -fno-inline-functions -ggdb3

COSISTATICSOURCES=$(COSISOURCES) $(SD)/cosimain.cc $(libprogopt_la_SOURCES) $(libcosifilesystem_la_SOURCES) $(libcosiserialization_la_SOURCES) $(libcositimer_la_SOURCES) $(libcosichrono_la_SOURCES) $(SD)/../cosi_rand/mtwist.c $(SD)/../cosi_rand/poisson.c $(SD)/../cosi_rand/ranbinom.c $(SD)/../cosi_rand/gammln.c $(SD)/../cosi_rand/rng.c $(SD)/../cosi_rand/expdev.c $(SD)/../cosi_rand/gamma.c

coalescent_static: $(COSISTATICSOURCES)
	g++ -v -static -o $@ -I$(SD)/.. -I$(SD)/boost $(COSI_CPPFLAGS) $(CPPFLAGS) $(COSI_CXXFLAGS) $(CXXFLAGS)  -Wno-error -pipe $(COSISTATICSOURCES) -lrt -lpthread

# End section: Building cosi


###########################################################################################################################

#
# Section: Tests
#

#
# Section: Building documentation
#

make-docs:
	$(mkdir_p) $(builddir)/doc/html
	$(srcdir)/doc/NaturalDocs/NaturalDocs -p $(srcdir)/doc/developers/ndproj -i $(SD) -o FramedHTML $(builddir)/doc/html

.PHONY: make-docs

COSIFLAGS = $(builddir)/coalescent.flags.txt

$(COSIFLAGS): $(builddir)/coalescent$(EXEEXT)
	echo "coalescent " $(coalescent_CPPFLAGS) $(CPPFLAGS) $(coalescent_CXXFLAGS) $(CXXFLAGS) $(coalescent_LDFLAGS) $(LDFLAGS) > $@
	echo "libcosi " $(libcosi_la_CPPFLAGS) $(CPPFLAGS) $(libcosi_la_CXXFLAGS) $(CXXFLAGS) $(libcosi_la_LDFLAGS) $(LDFLAGS) >> $@
	$(CXX) --version >> $@
if GITSTATE
	-git status $(srcdir) >> $@
	-git branch --no-abbrev -v
	-git diff $(srcdir)	>> $@
endif


all-local: $(COSIFLAGS)

MOSTLYCLEANFILES = $(COSIFLAGS)


