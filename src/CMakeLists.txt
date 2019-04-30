#Make sure everything installed to the right place
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")



include_directories("${brr_INCLUDE_DIRS}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../BayesRRcmd/tbb/include")

set(SourceFiles  Rwrapper.cpp RcppExports.cpp)

add_library(AMIGO SHARED ${SourceFiles} )
target_link_libraries(AMIGO PRIVATE Rcpp libbrr-dev  ${TBB_IMPORTED_TARGETS})
target_include_directories(AMIGO PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../inst/include ${CMAKE_CURRENT_SOURCE_DIR})

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)

