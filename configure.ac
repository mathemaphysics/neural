#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(neural, 0.1, ryan-p-daly@uiowa.edu)
AM_INIT_AUTOMAKE(neural, 0.1)
AC_CONFIG_SRCDIR([rbmrun.c])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC

PKG_CHECK_MODULES( util, util )

AC_CONFIG_FILES( [Makefile] [neural.pc] )
AC_OUTPUT
