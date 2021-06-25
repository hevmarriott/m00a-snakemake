set -e
filelist=$1
pm_pop -b whamg.sh -r ${filelist}
pm_pop -b manta.sh -r ${filelist}
pm_pop -b delly.sh -r ${filelist}
pm_pop -b melt.sh -r ${filelist}
pm_stage &
sleep 60
pm_submit -n
wait
pm_submit -n
