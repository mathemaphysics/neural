#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([neural],[0.1],[mathemaphysics@gmail.com])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([rbmrun.c])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
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

