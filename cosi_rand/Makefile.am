lib_LTLIBRARIES += libcosi_rand.la

RD=$(srcdir)/cosi_rand

libcosi_rand_la_SOURCES = $(RD)/expdev.c $(RD)/gamma.c $(RD)/gammln.c $(RD)/mtwist.c $(RD)/multinom.c $(RD)/poisson.c \
	$(RD)/ranbinom.c $(RD)/rng.c $(RD)/expdev.h $(RD)/gamma.h $(RD)/gammln.h $(RD)/mtwist.h $(RD)/multinom.h \
	$(RD)/poisson.h $(RD)/ranbinom.h $(RD)/random.h
libcosi_rand_la_CPPFLAGS = -I$(srcdir)

noinst_LTLIBRARIES += libcosi_rand_noinst.la
libcosi_rand_noinst_la_SOURCES = $(libcosi_rand_la_SOURCES)
libcosi_rand_noinst_la_CPPFLAGS = $(libcosi_rand_la_CPPFLAGS)







