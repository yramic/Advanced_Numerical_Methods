file(REMOVE_RECURSE
  "libquadrature.pdb"
  "libquadrature.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/quadrature.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
