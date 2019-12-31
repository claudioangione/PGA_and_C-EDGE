num_populations = 0;
model_name = 'recon2_merged_bio_PHGDH';

M = 2; %number of objective functions
load('geni.mat');
V = numel(geni);

cd(['populations_' model_name]);
while (exist(['./population' num2str(num_populations + 1) '.mat'], 'file') == 2)  % check the existance and counts how many population files are in the folder
    num_populations = num_populations +1;
end 
disp([num2str(num_populations) ' populations detected.'])



primo_file=load('population1.mat');
dim = size(primo_file.chromosome);
num_cols = dim(2);
num_rows = dim(1);
primo_obj = num_cols-3;
secondo_obj = num_cols-2;

disp([num2str(num_rows) ' individuals in each population.'])


%matrice=[];
matrice=zeros(num_rows*num_populations, num_cols);



for i=1:num_populations
    current_file = load(strcat('population',num2str(i),'.mat'));
    chromosome = current_file.chromosome;
    matrice((i-1)*num_rows+1:i*num_rows,:) =  chromosome;
    iterations_left=num_populations-i
end

cd ..
    %matrice = chromosome(:,primo:secondo);

array_aux = zeros(num_populations * num_rows,1);
array_aux_2 = zeros(num_populations * num_rows,1);
for i = 1:num_populations 
    array_aux((i-1)*num_rows+1 : i*num_rows) = i;
end
for i = 1:num_rows 
    array_aux_2(i:num_rows:end) = i;
end
matrice_tracked = [matrice array_aux array_aux_2];  % add last columns to tell (1) from which population the row comes from, and (2) which position it occupied in the population. the number of the comlumns for the first and second objectives are the same.
clear matrice

matrice_sorted = sortrows(matrice_tracked,primo_obj); %sort matrice according to the first objective


obj_sorted = matrice_sorted(:,[end-(4+M-1):end-4 end-1 end]);  %takes the M objectives, plus the two columns (at the end of the matrix) added by me a few lines before

obj_non_redund = zeros(size(obj_sorted,1),size(obj_sorted,2));   %initialised at its maximum possible size, but later we eill eliminate all the final rows that are left with only zeros, and therefore have not been used to save the non redundant populations
j = 1;

for i=2:size(obj_non_redund,1)
    if abs(obj_sorted(i,1)-obj_sorted(i-1,1))~=0 || abs(obj_sorted(i,2)-obj_sorted(i-1,2))~=0 || abs(obj_sorted(i,3)-obj_sorted(i-1,3))~=0  %for three objectives
        obj_non_redund(j,:)=obj_sorted(i,:);
        j = j+1;
    end
    iterations_left=size(obj_non_redund,1)-i
end
       
obj_non_redund = obj_non_redund(sum(abs(obj_non_redund'))~=0,:); % selects only those rows such that the elements are not all zeros, i.e. the sum of the absolute value of the row is nonzero. THis eliminates rows where all the elements are 0, and therefore have not been used from the initialized obj_not_refund matrix



feasible = [obj_non_redund(:,1:M),zeros(size(obj_non_redund,1),1),obj_non_redund(:,end-1:end)];
non_dominated_app = non_domination(feasible,M,0);  %non_domination looks for minimisation ,so keep the objectives negative if you want to define the non_dominated points as those points that maximise those objectives

indici_non_dom = find(non_dominated_app(:,M+1)==1);
non_dominated = non_dominated_app(indici_non_dom,:);
indici_others = find(non_dominated_app(:,M+1)==0);
others = non_dominated_app(indici_others,:);

%treshold = 1;
%non_dominated = non_dominated(find((non_dominated(:,2) > treshold)==1),:);    %selects only certain rows of the non_dominated and others
%others = others(find((others(:,2) > treshold)==1),:);



addpath(genpath('./'));
non_dominated_chromosomes = zeros(size(non_dominated,1),numel(geni));
n_pop = non_dominated(:,end-1);
n_ind = non_dominated(:,end);
for i=1 : size(non_dominated,1)
    disp(['non dominated chromosomes' num2str(length(non_dominated)-i)]);
    %sol(i) = ['populations_recon2_merged_bio_IDH1-2/population' num2str(n_pop) '.mat'];
    sol=['populations_' model_name '/population' num2str(n_pop(i)) '.mat'];
    chromosome = load(sol);
    chromosome = chromosome.chromosome;
    non_dominated_chromosomes(i,:) = chromosome(n_ind(i),1:V);
end



others_chromosomes = zeros(size(others,1),numel(geni));
n_pop = others(:,end-1);
n_ind = others(:,end);
for i=1 : size(others,1)
    disp(['other chromosomes' num2str(length(others)-i)]);
    sol=['populations_' model_name '/population' num2str(n_pop(i)) '.mat'];
    chromosome = load(sol);
    chromosome=chromosome.chromosome;
    others_chromosomes(i,:) = chromosome(n_ind(i),1:V);
end

%%
% change the sign (the first or second column need to be changed according to the fact whether we
% have already changed it or not in evaluate_objective
     
% others(:,1) = - others(:,1);
% %others(:,4:end) = - others(:,4:end);
% non_dominated(:,1) = - non_dominated(:,1);
% %non_dominated(:,4:end) = - non_dominated(:,4:end);


%%

x = - non_dominated(:,1);
y = non_dominated(:,2);
z = non_dominated(:,3);
 
fig1=figure(1);
plot3(x,y,z,'*:', 'Color',[1 0 0])
grid on
 
title('');
xlabel('Biomass [h^{-1}]');
ylabel('Oncoflux [mmolh^{-1}gDW^{-1}]');
zlabel('Distance from WT gene expression');

hold on

%subplot(1,1,1)
%axes('position',[0.5 0.5 0.4 0.4]);

fig2=figure(2);
num_others = size(others,1);

%x = others(:,1);
%y = others(:,2);
%plot(x,y,'*','Color',[0 0 0])

colore=240;     % [200 200 200] is grey, almost white ([255 255 255]). Instead, [0 0 0] is black
for i=1:num_others
    x = - others(i,1);
    y = others(i,2);
    z = others(i,3);
    colore = 10 + 245*(1 - (others(i,M+2))/num_populations);  % others(i,4) contains the number (with the minus sign, though) of the population to which others(i,:) belongs. The +10 avoidt the points are too black and allows for more variability of colours
    plot3(x,y,z,'*','Color',[colore/255 colore/255 colore/255])
    hold on
end



if exist('non_dominated.mat')
    delete('non_dominated.mat')
end
if exist('others.mat')
    delete('others.mat')
end

save non_dominated.mat non_dominated;
save non_dominated_chromosomes.mat non_dominated_chromosomes;
save others.mat others;
save others_chromosome.mat others_chromosomes;


x = - non_dominated(:,1);
y = non_dominated(:,2);
z = non_dominated(:,3);
plot3(x,y,z,'*:', 'Color',[1 0 0])

grid on
hold off

%axes('position',[0 0 1 1]);
title('');
xlabel('Biomass [h^{-1}]');
ylabel('Oncoflux [mmolh^{-1}gDW^{-1}]');
zlabel('Distance from WT gene expression');

[h_m h_i]=inset(fig2,fig1);


%exportfig(gcf, 'figura.pdf', 'color', 'cmyk', 'Width', '12', 'FontMode', 'scaled', 'FontSize', '1' );
