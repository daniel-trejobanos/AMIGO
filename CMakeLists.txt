cmake_minimum_required(VERSION 2.8.12)
project(AMIGO)
add_subdirectory(BayesRRcmd)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/Rcpp/cmake)

find_package(R REQUIRED)
find_package(Rcpp REQUIRED)
set(CMAKE_SHARED_LIBRARY_PREFIX "")
add_subdirectory(src)
add_custom_target(copyPackage ALL)	
set(RFILES BRoptions.R RcppExports.R)
add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/DESCRIPTION.in" "${CMAKE_CURRENT_BINARY_DIR}/DESCRIPTION")
#Copy R files
foreach(RFILE IN LISTS RFILES)
  add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/R/${RFILE}" "${CMAKE_CURRENT_BINARY_DIR}/R/${RFILE}")
	endforeach()
#Run roxygen2
#	add_custom_command(TARGET copyPackage COMMAND "${R_COMMAND}" ARGS "-e" "\"library(roxygen2);roxygen2::roxygenise()\"" WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
install(CODE "execute_process(COMMAND \"${R_COMMAND}\" CMD INSTALL . WORKING_DIRECTORY \"${PROJECT_BINARY_DIR}\")")
