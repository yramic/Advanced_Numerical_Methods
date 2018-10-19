file(REMOVE_RECURSE
  "libgeometry.pdb"
  "libgeometry.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/geometry.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
