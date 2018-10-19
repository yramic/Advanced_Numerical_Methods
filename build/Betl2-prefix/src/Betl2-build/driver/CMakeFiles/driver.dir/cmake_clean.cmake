file(REMOVE_RECURSE
  "libdriver.pdb"
  "libdriver.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/driver.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
