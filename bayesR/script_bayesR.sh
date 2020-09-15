#!/bin/sh

app=sim_app
test=sim_test
app_out=sim_app
test_out=sim_test


ncat=4
varcat=0.0,0.0001,0.001,0.01
it=50000
burn=30000




./bayesR -bfile $app -out $app_out -seed 14 -ndist $ncat -gpin $varcat -vara $vara -numit $it -burnin $burn 

./bayesR -bfile $test -predict -ndist $ncat -gpin $varcat -out $test_out -model $app_out.model -freq $app_out.frq -param $app_out.param 
