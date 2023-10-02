#!/usr/bin/env bash
# This is a slightly modified version of compile_test.sh to fit gitlab ci workflow

task() {

    if [[ $d =~ /CMakeFiles/ ]]; then #should not be checked
      return
    fi

        cd "$d"
        [[ $d =~ ./CONVQUAD/(.*)/ ]]
        cmd_test="./${BASH_REMATCH[1]}_test_mastersolution";
        cmd_mastersolution="./${BASH_REMATCH[1]}_mastersolution";

        if [[ -f "$cmd_test" ]]; then
          echo "Executing $cmd_test";
          eval $cmd_test
        else
          echo "*** WARNING: No unit tests found in $d ***";
        fi

        cd ../..
}
set -e
for d in ./CONVQUAD/*/ ;
do
  task $d #parallelizing did not work, can be modified at a later stage
done
