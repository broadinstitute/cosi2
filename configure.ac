#
# configure.ac: cosi configuration
#
# Automake/autoconf configuration options for cosi
#

AC_PREREQ([2.68])
AC_INIT([cosi], [2.3.2rc7], [ilya_shl@alum.mit.edu])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([1.9 tar-pax parallel-tests -Wall])

#AM_PATH_PYTHON([2.5],, [:])
#AM_CONDITIONAL([HAVE_PYTHON], [test "$PYTHON" != :])
AX_PYTHON
AM_CONDITIONAL([HAVE_PYTHON], [test x${PYTHON_BIN} != x])

#
# Section: Optional features/modules
#
# Options to enable/disable some cosi features.
#

dnl AC_ARG_ENABLE([ehh],
dnl 		       [  --enable-ehh    Turn on support for EHH computation (experimental)],
dnl      [case "${enableval}" in
dnl        yes) ehh=true ;;
dnl        no)  ehh=false ;;
dnl        *) AC_MSG_ERROR([bad value ${enableval} for --enable-ehh]) ;;
dnl      esac],[ehh=false])
dnl AM_CONDITIONAL([EHH], [test x$ehh = xtrue])

dnl AC_ARG_ENABLE([coalapx],
dnl 		       [  --enable-coalapx    Turn on support for approximating the coalescent],
dnl      [case "${enableval}" in
dnl        yes) coalapx=true ;;
dnl        no)  coalapx=false ;;
dnl        *) AC_MSG_ERROR([bad value ${enableval} for --enable-coalapx]) ;;
dnl      esac],[coalapx=true])
dnl AM_CONDITIONAL([COALAPX], [test x$coalapx = xtrue])

dnl AC_ARG_ENABLE([condsnp],
dnl 		       [  --enable-condsnp    Turn on support for conditioning on the presence of a particular SNP],
dnl      [case "${enableval}" in
dnl        yes) condsnp=true ;;
dnl        no)  condsnp=false ;;
dnl        *) AC_MSG_ERROR([bad value ${enableval} for --enable-condsnp]) ;;
dnl      esac],[condsnp=true])
dnl AM_CONDITIONAL([CONDSNP], [test x$condsnp = xtrue])


AC_ARG_ENABLE([devopts],
		       [  --enable-devopts    (Internal) Enable development-only options],
     [case "${enableval}" in
       yes) devopts=true ;;
       no)  devopts=false ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-devopts]) ;;
     esac],[devopts=false])
AM_CONDITIONAL([DEVOPTS], [test x$devopts = xtrue])

AC_ARG_ENABLE([gitstate],
		       [  --enable-gitstate    Turn on saving of git state used to build cosi (devel-only)],
     [case "${enableval}" in
       yes) gitstate=true ;;
       no)  gitstate=false ;;
       *) AC_MSG_ERROR([bad value ${enableval} for --enable-gitstate]) ;;
     esac],[gitstate=false])
AM_CONDITIONAL([GITSTATE], [test x$gitstate = xtrue])


# End section: Optional features/modules

AM_PROG_AR
LT_INIT([disable-shared])

AC_PROG_CXX

AC_PROG_LN_S
AC_PROG_MKDIR_P

AC_LANG(C++)

AX_CXX_COMPILE_STDCXX_11

#AX_CXXFLAGS_WARN_ALL([AM_CXXFLAGS])
AC_SUBST([WARNFLAGS])
AX_CHECK_COMPILE_FLAG([-Werror], [WARNFLAGS="${WARNFLAGS} -Werror"], [], [])
AX_CHECK_COMPILE_FLAG([-Wall], [WARNFLAGS="${WARNFLAGS} -Wall"], [], [])
AX_CHECK_COMPILE_FLAG([-Wextra], [WARNFLAGS="${WARNFLAGS} -Wextra"], [], [])
AX_CHECK_COMPILE_FLAG([-Wno-error], [BOOST_CXXFLAGS="${BOOST_CXXFLAGS} -Wno-error"], [], [])
AX_CHECK_COMPILE_FLAG([-Wdouble-promotion], [WARNFLAGS="${WARNFLAGS} -Wdouble-promotion"], [], [])
#AX_CHECK_COMPILE_FLAG([--param max-inline-insns-single=1000 -Wno-unused-function -Wuninitialized -Wredundant-decls -Wpointer-arith -Wcast-qual -Wcast-align -Wformat -Woverloaded-virtual -Wstrict-null-sentinel -Wlogical-op 	-Wmissing-declarations -Winit-self], [WARNFLAGS="${WARNFLAGS} --param max-inline-insns-single=1000 -Wno-unused-function -Wuninitialized -Wredundant-decls -Wpointer-arith -Wcast-qual -Wcast-align -Wformat -Woverloaded-virtual -Wstrict-null-sentinel -Wlogical-op 	-Wmissing-declarations -Winit-self"], [], [])

#
# LTO (link-time optimization)
#
AC_ARG_ENABLE(
    [lto],
    AC_HELP_STRING(
        [--enable-lto],
        [
												Enables or Disables Linktime Code Optimization (LTO is enabled by default)
        ]
    ),
    [
        enable_lto="$enableval"
        case $enableval in
            "no");;
            "yes");;
            *) AC_MSG_ERROR([[invalid argument --enable-lto=$disableval... stopping]]);;
        esac
    ],
    [enable_lto="no"]
)

#
# LTO Support test
#
if test "$enable_lto" != "no" ; then

	 OLD_CXXFLAGS="$CFLAGS"
	 CXXFLAGS="$CXXFLAGS -flto"
	 
	 OLD_LDFLAGS="$LDFLAGS"
	 LDFLAGS="$LDFLAGS -flto"

	 AC_MSG_CHECKING([whether $CXX supports -flto])
	 AC_RUN_IFELSE(
				[ AC_LANG_SOURCE([[
							int main(int argc, char **argv){
              return 0;
           }
				]]) ],
				[
					AC_MSG_RESULT([yes])
				],
					[
					 AC_MSG_RESULT([no])
						CXXFLAGS="$OLD_CXXFLAGS"
						LDFLAGS="$OLD_LDFLAGS"
						],
																																																												[
																																																													  AC_MSG_RESULT([guessing no])
																																																												]
																																																												)
fi



AC_TYPE_UINTMAX_T

ROOT_PATH(, [AC_DEFINE([HAVE_ROOT])dnl
foundroot=true])
AM_CONDITIONAL([HAVEROOT], [test x$foundroot = xtrue])

#BOOST_REQUIRE([1.50.0])

#AM_CONDITIONAL([NUINU],[test x$ZDRASTEVAM = xyes])

#BOOST_FOREACH
#BOOST_MATH
#BOOST_UTILITY
#BOOST_SMART_PTR
#BOOST_LAMBDA
#BOOST_PROGRAM_OPTIONS
#BOOST_STRING_ALGO
#AM_COND_IF([SAMPLESTATS],
#           [BOOST_REGEX])
#BOOST_PROGRAM_OPTIONS
#BOOST_REGEX

LIBSEQUENCE_INC_DIR=$prefix/include
LIBSEQUENCE_LIB_DIR=$prefix/lib
AC_ARG_WITH(
	[libsequence],
	AC_HELP_STRING(
		[--with-libsequence=DIR],
		[Specify the root directory for LIBSEQUENCE installation.]
		),
	[
	if ! test "$withval" = "yes" -o "$withval" = "no" ; then
		LIBSEQUENCE_INC_DIR="$withval/include"
		LIBSEQUENCE_LIB_DIR="$withval/lib"
	fi
	])

AM_CONDITIONAL([HAVELIBSEQUENCE], [test -f "$LIBSEQUENCE_LIB_DIR/libsequence.a"])
AC_SUBST(LIBSEQUENCE_INC_DIR)
AC_SUBST(LIBSEQUENCE_LIB_DIR)

AC_ARG_WITH([gtest],
  [AS_HELP_STRING([--with-gtest], [location of Google Test source])],
  [GTEST_DIR=$withval])
AC_SUBST([GTEST_DIR])
AM_CONDITIONAL([HAVE_GTEST], [test -f "$GTEST_DIR/include/gtest/gtest.h"])

AC_ARG_WITH([gmock],
  [AS_HELP_STRING([--with-gmock], [location of Google Mock source])],
  [GMOCK_DIR=$withval])
AC_SUBST([GMOCK_DIR])

#AC_MSG_CHECKING([whether CXX supports latest gcc flags])
#my_save_cxxflags="$CXXFLAGS"
#COSI_TRY_CXXFLAGS="--param max-inline-insns-single=1000 -Werror -Wall -Wextra -Wno-unused-function -Wuninitialized -Wredundant-decls -Wpointer-arith -Wcast-qual -Wcast-align -Wno-long-long -Wformat-security -Wundef -Woverloaded-virtual -Wstrict-null-s#entinel -Wshadow -Wstrict-overflow=5 -Wlogical-op -Wmissing-declarations -Wdouble-promotion -Winit-self"
#CXXFLAGS="$COSI_TRY_CXXFLAGS"
#AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
#    [AC_MSG_RESULT([yes])]
#    [COSI_CXXFLAGS=$(COSI_TRY_CXXFLAGS)],
#    [AC_MSG_RESULT([no])]
#)
#CXXFLAGS="$my_save_cxxflags"
#AC_SUBST([COSI_CXXFLAGS])

AC_SEARCH_LIBS([clock_gettime],[rt posix4])

AM_MAINTAINER_MODE([disable])

AC_CONFIG_FILES([Makefile])
#AC_CONFIG_SUBDIRS([libsigc++-2.3.1])
AC_OUTPUT
