#!/bin/tcsh
### Job name job (default is name of pbs script file)
#PBS -N SIZE_vs_I_20
###set u id, check it with id command
### Number and which of nodes as well as number of CPUs per node:
#PBS -l nodes=1:ppn=16
#
### Declare job non-rerunable
#PBS -r n
### resource limits: amount of memory and CPU time ([[h:]m:]s).
#PBS -l walltime=23:00:00
#
### Output files - by default, <scriptname>.[eo]NN where NN is the
### job ID
#PBS -e /home/dscherli/SIZE_vs_I/ALQ_20_0.001.1/ERROR
#PBS -o /home/dscherli/SIZE_vs_I/ALQ_20_0.001.1/SALIDA
### Queue name (default, ib, ...)
#PBS -q gpu
#
#############
### This job\'s working directory
echo \"Working directory is $PBS_O_WORKDIR\"
cd $PBS_O_WORKDIR
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
# Run your executable
setenv LD_LIBRARY_PATH /opt/cuda/6.5/lib64:/opt/intel/2013/lib/intel64:/opt/intel/2013/mkl/lib/intel64:/lib64:/lib:/usr/lib64:/usr/lib:/usr/local/lib64:/usr/local/lib:/home/dscherli/lio/lioamber:/home/dscherli/lio/g2g
/home/lucas/jpm/transport_EPS/liosolo -i 20.in -b 631 -c 20.xyz -v > salida

