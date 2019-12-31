%run this script if the RAM memory is low, and the parallel pool will be shut
%down and re-open after every generation of a population, and the memory
%will be cleared

delete(gcp('nocreate'))

model_name = 'recon2_merged_bio_lactate';

num_populations = 0;
while (exist(['populations_' model_name '/population' num2str(num_populations + 1) '.mat'], 'file') == 2)  % check the existance and counts how many solution files are in the folder
    num_populations = num_populations +1;
end
disp([num2str(num_populations) ' populations detected.'])


for i=num_populations+1:384
    RUN(128,i,4)
    delete(gcp('nocreate'))
    clear all hidden
end