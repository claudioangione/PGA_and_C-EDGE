%% Runs the EDGE_k algorithm (the single-objective of EDGE_K is in the file genetic_operator.m
%%

function aux = RUN_CEDGEk(num_cores);

warning('off', 'all')
format long
addpath(genpath('C:\Program Files\MATLAB\R2014b\toolbox\cobra'));
%addpath(genpath('..\cobra'));
%initCobraToolbox

model_name = 'recon2_merged_bio_IDH1-2';

num_populations = 0;
while (exist(['populations_' model_name '/population' num2str(num_populations + 1) '.mat'], 'file') == 2)  % check the existance and counts how many solution files are in the folder
    num_populations = num_populations +1;
end 
disp([num2str(num_populations) ' populations detected.'])

num_objectives = 2;
%num_cores = 4;
addAttachedFiles(gcp,{'glpk.m','glpkcc'});  %adds these files to the independent workers of the parallel pool

expFBA_Boolean(128,384,model_name,num_objectives,num_populations,num_cores);
