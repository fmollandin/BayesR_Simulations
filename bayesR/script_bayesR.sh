#!/bin/sh
#$ -m e
#$ -M fanny.mollandin@inrae.fr
#$ -l h_vmem=1G
#$ -q longq

#cat=upstream
#perf=perf1

#input_app=$cat'_'$perf'_app'
#input_test=$cat'_'$perf'_test'
#output_app=$cat'_'$perf'_app'
#output_test=$cat'_'$perf'_test'

scenario=$1
id=$2

app=sim_$scenario$id'_app'
test=sim_$scenario$id'_test'
app_out=sim_$scenario$id'_x0.5_app'
test_out=sim_$scenario$id'_x0.5_test'


ncat=4
varcat=0.0,0.0001,0.001,0.005
vara=0.8
it=50000
burn=30000




../bayesR/bayesR_Fanny2 -bfile $app -out $app_out -seed 14 -ndist $ncat -gpin $varcat -vara $vara -numit $it -burnin $burn 

../bayesR/bayesR_Fanny2 -bfile $test -predict -ndist $ncat -gpin $varcat -out $test_out -model $app_out.model -freq $app_out.frq -param $app_out.param 