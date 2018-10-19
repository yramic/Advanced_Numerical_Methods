file(REMOVE_RECURSE
  "libsparse_operators.pdb"
  "libsparse_operators.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/sparse_operators.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
