#-----------------------------------------------------------------------------------
# - A Python-wrapping module for CMake using PYB11Generator/pybind11
#
# Modified from the UseSWIG module from the cmake2.8 release.
#
# Defines the following macros:
#   PYB11_GENERATE_BINDINGS(module_name module_list)
#     - Generates the Python bindings for each module in the list
#
# Internal macros to do the work:
#   PYB11_REGISTER_MODULE(name)
#     - Builds list of Python source code and generated C++ code
# 
#
# Variables that must be set before calling PYB11_GENERATE_BINDINGS:
#   PYB11_DIR
#     - Location of the python source code to build the bindings
#   PYTHON_EXECUTABLE
#     - Python executable
#   PYTHON_LIB_DIR
#     - Python lib (typically obtained from the executable root)
#
# Optional variables
#   PYB11_ADDITIONAL_ARGS
#     - Additional arguments to follow the command in PYB11_GENERATE_BINDINGS
#
# To get the names of the generated source
# use: ${PYB11_GENERATED_SOURCE}
#-----------------------------------------------------------------------------------

#
# Runs "python polytopeMOD.py" to generate the bindings.
# This is where the PYB11Generator module is actually used. Calls the previous
# internal macros to generate the source file lists
#
macro(PYB11_GENERATE_BINDINGS)
  # Make a list of Python modules that our MODule depends on.
  execute_process(COMMAND python3 -c "from modulefinder import ModuleFinder as MF; finder = MF(); finder.run_script('polytopeMOD.py'); mods = [x[1].__file__ for x in finder.modules.items() if x[0] != '__main__']; for m in mods: print m + ';',"
                  OUTPUT PYB11_ЅOURCE)
  # The generated C++ bindings have the same names with a different suffix
  string(REPLACE ".py" ".cc" PYB11_GENERATED_SOURCE ${PYB11_SOURCE})
  # Append the top-level module/generated bindings.
  list(APPEND PYB11_SOURCE "${PYB11_MODULE_NAME}MOD.py")
  list(APPEND PYB11_GENERATED_SOURCE "${PYB11_MODULE_NAME}.cc")

  # Here's the rule that generates the bindings.
  add_custom_command(
    OUTPUT ${PYB11_GENERATED_SOURCE} ${PYB11_GENERATED_HEADER}
    COMMAND ${PYTHON_EXECUTABLE} 
      -c 'from PYB11Generator import * \; import sys \; sys.path.append(\"${PYB11_DIR}\") \; sys.path.append(\"${PROJECT_BINARY_DIR}/src/PYB11\") \; import ${PYB11_MODULE_NAME}MOD \; PYB11generateModule(${PYB11_MODULE_NAME}MOD, \"${PYB11_MODULE_NAME}\" ) ' 
    DEPENDS ${PYB11_SOURCE}
    )
endmacro()
