Fortran code writen by Didier Boichard (INRAE) based on real genotype data, select randomly a defined number of QTL -low, medium & large- and give their an effect based on chosen part of variance.
<br />
Compile the source code:<br />
_gfortran simul_ped_v2.f -o simul_ped_<br />
and use aso so : <br />
_./simul_ped < param_sim_

The param file should include, in this order: <br />
5 0 8500                                #number of big, medium, small QTLs <br />
0.8 1                                   #heritability, part of genetic variance explained by QTLs <br />
0.03 0 0.0001                           #part of variance explained by each big, medium and small QTL (or relative variance between big, medium, low and low variances) <br />
29                                      #number of chromosome, 0 one phase file <br />
'/myfiles/ped_seq'                      #path of pedigree file <br />
'/myfiles/phase_50K_'                   #if number of chromosome=0, 'path/name' of the files, if not all the phase phase should have the name 'path/name'1 for chromosome 1, etc <br />
out_phase_                              #markers files, for each chromosome <br />
info_QTL                                #information on simulated QTLs <br />
phase                                   #genotyping file name <br />
simperf                                 #simulated phenotypes file name <br />
frq_sim                                 #SNP frequency file name <br />
0.15                                    #minimum MAF <br />
o                                       #Do we keep QTL in the markers file? o (Oui/Yes) n (Non/No)  <br />
p                                       #Input: P or p for phases files (2 lines/individual), T or t for typages (1 line/individual) <br />
p                                       #Output: P or p for phases files (2 lines/individual), T or t for typages (1 line/individual) <br />
1                                       #seed <br />
<br />
The file "simperf" contains the simulated phenotypes for each individual, and the file "info_QTL" contains informations on the simulated QTLs. <br /><br />
simperf: id / phenotype <br />
info_QTL: QTL count / chromosome / index on the chromosome / index on the genome / variance / frequency / effect
