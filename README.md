# DNAorigins
Supplementary code for paper "Bayesian inference of origin firing time distributions, origin interference and licensing probabilities from NGS data".

The repository consists of two .cpp source files essential for running an algorithm: algorithm_with_NA.cpp, algorithm_plain.cpp  

algorithm_with_NA.cpp is designed to handle the data with unsequencable region, algorithm_plain.cpp - without it.

The files are to be compiled with C++11

Instructions for running

make

It will create two executables ap and ap_NA for data without unsequencable region and with one respectively

./ap -h displays

Allowed Options:
  -h [ --help ]         Help screen
  -c [ --chr ] arg      Enter the chromosome number
  -r [ --rat ] arg (=0) Rat or not, rat==2 means human data
  -w [ --wt ] arg (=0)  If Rat, WT or not
  -l [ --left ] arg     Number of left origin (in case of human data 1 or 2)

Example of running ./ap -c 10 -l 1011

./ap_NA -h displays

Allowed Options:
  -h [ --help ]         Help screen
  -r [ --rat ] arg (=0) Rat or not
  -w [ --wt ] arg (=0)  If Rat, WT or not
  -l [ --left ] arg     Number of left origin

Example of running ./ap_NA -l 1014


Processed data files are provided. Where not stated otherwise data are obtained from paired-end NGS Okazaki Fragment data

chrX_f/r.txt

X - chromosome number, f corresponds to the data on the forward strand, r - to the data on the reverse strand

left/right indicates whether the data corresponds to the left or right triplet as described in the paper

For all chromosomes besides the human one the number to the left from the chromosome number corresponds to the name of the left origin in the triplet of origins (e.g. if this number is 1013 then the corresponding eft origin is ARS1013); for human chromosome 2 we ended left triplet as 1 and right triplet as 2; wt/rat1 corresponds to the datasets analysed from single-end sequenced data. Wt corresponds to wt, rat1 - to temperature sensitive mutant;

chr2_207 ARS207,207.5,207.8

chr7_717 ARS717, 718, 719

chr7_718 ARS718, 719, 720

chr8_813 ARS813, 815, 816

chr8_815 ARS815, 816, 818

chr10_1001 ARS1001, 1004,1005

chr10_1004 ARS1004, 1005,1006

chr10_1005 ARS1005, 1006,1007

chr10_1006 ARS1006, 1007,1007.5

chr10_1007 ARS1007, 1007.5,1008

chr10_10075 ARS1007.5, 1008,1009

chr10_1008 ARS1008, 1009,1010

chr10_1009 ARS1009, 1010,1005

chr10_1010 ARS1010, 1011, 1013

chr10_1011 ARS1011, 1013, 1014

chr10_1013 ARS1013, 1014, 1015

chr10_1014 ARS1014,1015,1018 missing values between 1098 and 1336 for paired-end sequenced data, 1108 and 1336 for wt/rat1

chr10_1015 ARS1015, 1018, 1019 missing values 588 and 826 for paired-end sequenced data, 605 and 833 for wt/rat1

chr10_1018 ARS1018, 1019, 1021

chromosome 2, human HeLa cells, 98.25-99.3 Mbp created from pile-up data from the corresponding paper:

chr2_1  

chr2_2




./ap -c 8 -r 0 -l 813

/Applications/MATLAB_R2020a.app/bin/matlab -nosplash -nodesktop -r "run('/Users/alinabazarova/Downloads/DNAorigins-master/visualise.m');"



