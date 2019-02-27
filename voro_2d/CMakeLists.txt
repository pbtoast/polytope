#add_library(voro_2d common.cc cell_2d.cc container_2d.cc v_base_2d.cc v_compute_2d.cc
#            c_loops_2d.cc wall_2d.cc cell_nc_2d.cc ctr_boundary_2d.cc)
add_library(voro_2d
            ${PROJECT_SOURCE_DIR}/voro_2d/common.cc
            ${PROJECT_SOURCE_DIR}/voro_2d/cell_2d.cc
            ${PROJECT_SOURCE_DIR}/voro_2d/container_2d.cc
            ${PROJECT_SOURCE_DIR}/voro_2d/v_base_2d.cc
            ${PROJECT_SOURCE_DIR}/voro_2d/v_compute_2d.cc
            ${PROJECT_SOURCE_DIR}/voro_2d/c_loops_2d.cc
            ${PROJECT_SOURCE_DIR}/voro_2d/wall_2d.cc
            ${PROJECT_SOURCE_DIR}/voro_2d/cell_nc_2d.cc
            ${PROJECT_SOURCE_DIR}/voro_2d/ctr_boundary_2d.cc)

# Library install targets
install (TARGETS voro_2d DESTINATION ${PROJECT_BINARY_DIR}/lib)
