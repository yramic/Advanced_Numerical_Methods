file(REMOVE_RECURSE
  "libinput_interface.pdb"
  "libinput_interface.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/input_interface.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
