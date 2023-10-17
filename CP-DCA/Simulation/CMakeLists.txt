CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(CPDCA)

include_directories(include)

aux_source_directory(src/WBMatrix CPDCA_SRC)
aux_source_directory(src CPDCA_SRC)

add_library(CPDCA_LIB ${CPDCA_SRC})

add_executable(CPDCA test/attack.c)

target_link_libraries(CPDCA CPDCA_LIB m)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O2")