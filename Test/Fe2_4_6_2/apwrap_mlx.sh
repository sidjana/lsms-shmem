#export VT_FILTER_SPEC=/ccs/home/sjh/tests/microbenchmarks/shmem_vttau_consttime/filter_file
#export VT_PLUGIN_CNTR_METRICS="VTRAPL_*"
export VT_METRICS=PAPI_TOT_INS,PAPI_L3_TCM,PP0_ENERGY:PACKAGE0,DRAM_ENERGY:PACKAGE0
export VT_METRICS_SEP=","
export VTRAPL_INTERVAL=10000
#export LD_PRELOAD=/ccs/home/sjh/opt/vt_install/lib/libosh_wrap.so:/ccs/home/sjh/opt/vt_install/lib/libvt-mt.so
export VT_BUFFER_SIZE="3G"
export VT_MAX_FLUSHES=15
export VT_MODE=
export VT_VERBOSE=2
export LD_LIBRARY_PATH=/ccs/proj/stf010/opt/papi/papi_install/lib/:/ccs/home/sjh/opt/vt_install/lib:/ccs/home/sjh/opt/ctools/ctool_2.12//lib:/ccs/home/sjh/opt/vt-plugin/:/ccs/proj/stf010/opt/openshmem_install/lib/

export VT_FILE_PREFIX=$1.$(/bin/hostname).$$
echo $VT_FILE_PREFIX >> $1.apwrap

#export VT_GNU_NMFILE=./nm_log

#If VampirTrace has been built with resource usage support, it 
#is able to record this information as performance counters to the trace.
#export VT_RUSAGE=all
#export VT_RUSAGE INTV=50
export VT_CPUIDTRACE=yes

#vtcc -DVTRACE_PTHREAD -DVTRACE -vt:cc oshcc 
#To enable tracing of all C-Pthread API functions include the header vt_user.h

exec $@
#exec hwloc-bind socket:0.core0 $@

