if(EXISTS "/Users/bobschreiner/ETH/AdvNumCse/Code/cmake-build-debug/HMAT/LowRankMerge/LowRankMerge_test_mysolution[1]_tests.cmake")
  include("/Users/bobschreiner/ETH/AdvNumCse/Code/cmake-build-debug/HMAT/LowRankMerge/LowRankMerge_test_mysolution[1]_tests.cmake")
else()
  add_test(LowRankMerge_test_mysolution_NOT_BUILT LowRankMerge_test_mysolution_NOT_BUILT)
endif()
