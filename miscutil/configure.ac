AC_PREREQ([2.68])
AC_INIT([cosi_rand], [2.0], [ilya_shl@alum.mit.edu])
AC_CONFIG_SRCDIR([utils.h])
#AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wall -Werror])
#AC_CONFIG_MACRO_DIR([m4])
AM_PROG_AR
LT_INIT
AM_PROG_CC_C_O
AC_PROG_CXX
AC_PROG_LN_S
AC_CHECK_HEADERS([limits.h stdlib.h sys/time.h sys/timeb.h unistd.h])
AC_CHECK_FUNCS([floor gettimeofday pow sqrt])
AM_MAINTAINER_MODE([enable])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
