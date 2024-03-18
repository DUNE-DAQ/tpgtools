#!/bin/bash

integtest_list=( "test_tpg_implementation" )

usage() {
  declare -r script_name=$(basename "$0")
  echo """
Usage:
"${script_name}" [options(s)]

Options:
    -h, --help : prints out usage information 
    -d, --default : run default configuration 
"""
  let counter=0
  echo "List of available test:"
  echo ""
}

TEMP=`getopt -o hd --long help,stop-on-failure -- "$@"`
eval set -- "$TEMP"

while true; do
  case "$1" in 
    -h|--help) 
      usage
      exit 0;
      ;;
    -d|default)
      shift
      break
      ;;
    --stop-on-failure)
      let stop_on_failure=1
      shift
      ;;
    --)
      shift
      break
      ;;  
  esac
done

export PYTEST_DISABLE_PLUGIN_AUTOLOAD=1
tests=../tests
pattern_name=(patt_golden)
channel_number=(0 10 32 63)
clock_tick=(1 10 53 63)
threshold=(499)
num_frames=(2)
alg=(SimpleThreshold AbsRS)

#pytest --printout=more -k "test_all_params[patt_golden-32-31-64-2-AbsRS]" ../tests
#pytest --printout=more -k test_all_params[${pattern_name}-${channel_number}-${clock_tick}-${threshold}-${num_frames}-${alg}] ${tests}

# run all tests
for p in ${pattern_name[@]}
do 
  for c in ${channel_number[@]}
  do
    for k in ${clock_tick[@]}
    do
      for t in ${threshold[@]}
      do
        for n in ${num_frames[@]}
        do
          for a in ${alg[@]}
          do
          pytest --printout=more -k test_all_params[${p}-${c}-${k}-${t}-${num_frames}-${a}] ${tests}
          done
       	done
      done
    done
  done
done



