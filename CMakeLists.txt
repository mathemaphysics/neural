cmake_minimum_required (VERSION 2.6)
project (neural)

set (tutorial_VERSION_MAJOR 1)
set (tutorial_VERSION_MINOR 0)

add_executable(rbmrun rbmrun.c)
target_link_libraries(rbmrun util)

include_directories("${PROJECT_BINARY_DIR}")

configure_file (
	"${PROJECT_SOURCE_DIR}/neural.h.in"
	"${PROJECT_BINARY_DIR}/neural.h"
	)

install (TARGETS rbmrun DESTINATION bin)
install (FILES "${PROJECT_BINARY_DIR}/neural.h" DESTINATION include)

# vim: set tabstop=4:expandtab:smarttab
