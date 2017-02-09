#/home/uriel/progs/amber/lib/:/home/uriel/progs/amber/lib64/:/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/intel64/:/opt/intel/composer_xe_2013_sp1.3.174/mkl/lib/intel64/:/usr/local/cuda-6.0/lib64:/home/uriel/progs/liogit/liotransport2/lio/lioamber/:/home/uriel/progs/liogit/liotransport2/lio/g2g/:/home/uriel/progs/g09/bsd:/home/uriel/progs/g09/private:/home/uriel/progs/g09:/usr/local/cuda-6.0/lib64:/opt/intel/composer_xe_2013_sp1.3.174/mkl/lib/intel64/:/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/intel64/:/usr/local/cuda-6.0/src/
#/home/uriel/progs/amber/lib/:/home/uriel/progs/amber/lib64/:/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/intel64/:/opt/intel/composer_xe_2013_sp1.3.174/mkl/lib/intel64/:/usr/local/cuda-6.0/lib64:/home/uriel/progs/liogit/liotransport/lioamber/:/home/uriel/progs/liogit/liotransport/g2g/:/home/uriel/progs/g09/bsd:/home/uriel/progs/g09/private:/home/uriel/progs/g09:/usr/local/cuda-6.0/lib64:/opt/intel/composer_xe_2013_sp1.3.174/mkl/lib/intel64/:/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/intel64/:/usr/local/cuda-6.0/src/
#runprog=/home/uriel/progs/liogit/lio-william/lio/liosolo/liosolo
#runprog=/home/uriel/progs/liogit/lio/liosolo/liosolo
#runprog=/home/uriel/progs/liogit/lio-cublas/liosolo/liosolo
runprog=/home/lucas/jpm/transport_EPS/liosolo/liosolo
#runprog=/home/uriel/progs/liogit/lio-rc/lio/liosolo/liosolo

${runprog} -i 20.in -c 20.xyz -ib -bs 'SBKJC' -fs 'DZVP Coulomb Fitting' -b basis -v > salida.juan
#$LIOBIN -i cl2.in -ib -bs 'SBKJC' -fs 'DZVP Coulomb Fitting' -b basis -c cl2.xyz -v > $SALIDA
