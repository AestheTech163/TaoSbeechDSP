cmake_minimum_required(VERSION 3.25)
project(TaoSbeechDSP C)

set(CMAKE_C_STANDARD 11)

include_directories(. common f0tracker dsp_utils aec_filters mmse_noise)

add_definitions(-DDEBUG=1)
add_definitions(-DUSE_FLOAT=1)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -O3")

add_library(common_lib STATIC
        common/arch.h
        common/vector.h
        common/vector.c
        common/numeric.h
        common/numeric.c
        common/reduce.h
        common/reduce.c
        common/io_utils.c
        common/io_utils.h)

add_library(dsp_utils STATIC
        dsp_utils/smallft.c
        dsp_utils/smallft.h
        dsp_utils/filterbank.h
        dsp_utils/filterbank.c
        dsp_utils/dsp_common_func.h
        dsp_utils/dsp_common_func.c
)

add_library(aec_filters STATIC
        aec_filters/common_states.h
        aec_filters/common_states.c
        aec_filters/aec_dual.h
        aec_filters/aec_dual.c
        aec_filters/mdf.h
        aec_filters/mdf.c
        aec_filters/aec_kalman.h
        aec_filters/aec_kalman.c
        aec_filters/aec_pnlms.h
        aec_filters/aec_pnlms.c
        )

add_library(f0tracker STATIC
        f0tracker/f0tracker.h
        f0tracker/f0tracker.c
        f0tracker/lpc.h
        f0tracker/lpc.c
        )

add_library(mmse_noise STATIC
#        dsp_utils/filterbank.h
#        dsp_utils/filterbank.c
        mmse_noise/mmse_noise.h
        mmse_noise/mmse_noise.c
        )

add_executable(demo_dualfilter demo_dualfilter.c)
target_link_libraries(demo_dualfilter common_lib dsp_utils aec_filters)

add_executable(demo_pnlms demo_pnlms.c)
target_link_libraries(demo_pnlms common_lib dsp_utils aec_filters)

add_executable(demo_kalman demo_kalman.c)
target_link_libraries(demo_kalman common_lib dsp_utils aec_filters)

add_executable(demo_f0tracker demo_f0tracker.c)
target_link_libraries(demo_f0tracker common_lib dsp_utils f0tracker )

add_executable(demo_mmsenoise demo_mmsenoise.c)
target_link_libraries(demo_mmsenoise common_lib dsp_utils f0tracker mmse_noise)

add_executable(test test.c)
target_link_libraries(test dsp_utils common_lib)