add_custom_target(doxygen_parametricbem2d
COMMAND doxygen Doxyfile
COMMENT "Generating 2DParametricBEM Documentation"
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doxygen
)
