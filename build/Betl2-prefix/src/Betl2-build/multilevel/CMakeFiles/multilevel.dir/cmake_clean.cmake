file(REMOVE_RECURSE
  "libmultilevel.pdb"
  "libmultilevel.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/multilevel.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
