cmake_minimum_required(VERSION 3.5)
project(transport-solver)


# Function to print lists of things.
function(print_list)
    foreach(item ${ARGV})
        message("${item}")
    endforeach()
endfunction()

# Macro to get list of subdirectories within current project.
macro(subdirlist result curdir)
    file(GLOB children RELATIVE ${PROJECT_SOURCE_DIR} ${curdir}/*)
    set(dirlist "")
    foreach(child ${children})
        if(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/${child})
            list(APPEND dirlist ${child})
        endif()
    endforeach()
    set(${result} ${dirlist})
endmacro()


# Assemble list of source files, excluding mains.
set(c_files "")
set(cpp_files "")

subdirlist(top_dirs "${PROJECT_SOURCE_DIR}/src/")

print_list(${top_dirs})

foreach(top_dir ${top_dirs})
    file(GLOB_RECURSE sources RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/${top_dir}/*.c)
    list(APPEND c_files ${sources})

    file(GLOB_RECURSE sources RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/${top_dir}/*.cpp)
    list(APPEND cpp_files ${sources})
endforeach()


add_executable(solver.x "src/main_solver.cpp" ${c_files} ${cpp_files})

# Get list of subdirectories of src/
file(GLOB top_files RELATIVE "${PROJECT_SOURCE_DIR}/src/" "${PROJECT_SOURCE_DIR}/src/*")



file(GLOB_RECURSE c_files   RELATIVE "${PROJECT_SOURCE_DIR}" *.c  )
file(GLOB_RECURSE cpp_files RELATIVE "${PROJECT_SOURCE_DIR}/src/*" *.cpp)


target_include_directories(solver.x
    PRIVATE ${PROJECT_SOURCE_DIR}/src
    )

target_compile_definitions(solver.x
    PRIVATE SPACE_DIMS=1
    )
