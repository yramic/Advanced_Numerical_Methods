file(REMOVE_RECURSE
  "libintegration.pdb"
  "libintegration.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/integration.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
