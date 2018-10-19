file(REMOVE_RECURSE
  "libbem_operator.pdb"
  "libbem_operator.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/bem_operator.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
