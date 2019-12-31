function aux = RUN(pop,gen,num_cores);

warning('off', 'all')
format long
%addpath(genpath('C:\Program Files\MATLAB\R2014a\toolbox\cobra'));
addpath(genpath('..\cobra'));
%initCobraToolbox

model_name = 'recon2_merged_bio_PHGDH';

num_populations = 0;
while (exist(['populations_' model_name '/population' num2str(num_populations + 1) '.mat'], 'file') == 2)  % check the existance and counts how many solution files are in the folder
    num_populations = num_populations +1;
end 
disp([num2str(num_populations) ' populations detected.'])

%parallel computing toolbox matlab
parpool(num_cores);

num_objectives = 2;
%num_cores = 4;
addAttachedFiles(gcp,{'glpk.m','glpkcc'});  %adds these files to the independent workers of the parallel pool

expFBA(pop,gen,model_name,num_objectives,num_populations,num_cores);
%expFBA(12,30,model_name,num_objectives,num_populations,num_cores);