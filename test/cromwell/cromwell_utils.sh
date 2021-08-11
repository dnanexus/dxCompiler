#!/bin/bash
# code that might help when moving tests from cromwell.

# find all the wdl files to use in tests
#find ../cromwell_failed -name '*.wdl' | rev | cut -d'/' -f 1 | rev | tr '\n' ' '


#rename: <name>.inputs -> <name>_input.json

#for i  in **/*.inputs; do # Whitespace-safe and recursive
#  mv "$i" `echo "$i" | sed 's/\.inputs/_input.json/'`
#done

# move all >not Succeeded< statuses to different folder (only >WDL<)

#for i in _tests/*.test; do # Whitespace-safe and recursive
#  metadata_status=`cat $i | grep status`
#  if [[ $metadata_status == *"Succeeded"* ]]; then
#    continue
#  else
#    workflow_file=`cat $i | grep "workflow:"`
#    if [[ $workflow_file == *".cwl"* ]]; then
#      continue
#    elif [[ -z $workflow_file ]]; then
#      continue
#    else
#      wdl_file=`echo $workflow_file | rev | cut -d':' -f 1  | rev | cut -c2- | cut -d'/' -f 1`
#      mv $wdl_file ../cromwell_failed/$wdl_file
#    fi
#  fi
#
#done



# Create results file from test files.

#for i in _tests/*.test; do # Whitespace-safe and recursive
#  workflow_file=`cat $i | grep "workflow:"`
#  if [[ $workflow_file == *".cwl"* ]]; then
#    continue
#  elif [[ -z $workflow_file ]]; then
#    continue
#  else
#    wdl_file_folder=`echo $workflow_file | rev | cut -d':' -f 1  | rev | cut -c2- | cut -d'/' -f -1`
#    wdl_file=`echo $workflow_file | rev | cut -d':' -f 1  | rev | cut -c2-`
#    if [ -d ../cromwell_failed/$wdl_file_folder ]; then
#      ./search_outputs.py $i ../cromwell_failed/$wdl_file
#    fi
#  fi
#done