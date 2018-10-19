file(REMOVE_RECURSE
  "libanalytical_functions.pdb"
  "libanalytical_functions.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/analytical_functions.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
