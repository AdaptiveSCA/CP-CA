CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(CPADCA)

include_directories(include)

aux_source_directory(src/m4ri CPADCA_SRC)

add_library(CPADCA_LIB ${CPADCA_SRC})

add_executable(CPADCA test/white_box_arx.c)

target_link_libraries(CPADCA CPADCA_LIB)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O2")