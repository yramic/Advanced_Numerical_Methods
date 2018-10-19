file(REMOVE_RECURSE
  "libbem_symmetries.pdb"
  "libbem_symmetries.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/bem_symmetries.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
