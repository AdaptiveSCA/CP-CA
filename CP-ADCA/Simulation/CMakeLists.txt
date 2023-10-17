CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(CPADCA)

include_directories(include)

aux_source_directory(src/WBMatrix CPADCA_SRC)
aux_source_directory(src CPADCA_SRC)

add_library(CPADCA_LIB ${CPADCA_SRC})

add_executable(CPADCA test/attack.c)

target_link_libraries(CPADCA CPADCA_LIB m)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O2")