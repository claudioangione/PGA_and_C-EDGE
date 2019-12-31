% Compute robustness using our generalised EDGE starting from a point in the objective space, and start perturbing all its genes to see on average how much that point is moving in the objective space. Considering that the only way to avoid weighting is multiplying the coordinates, I think that the correct measure in that case is evaluating the difference between two points by computing the absolute product of the difference between the two coordinates
% |(x1-x2)(y1-y2)|.
% In general, it's the volume of the hypercube whose opposite vertices are our two points on the Pareto.

load('geni.mat')
load('reaction_expression.mat')
load('recon2_merged_bio_IDH1-2.mat')
addpath(genpath('C:\Program Files\MATLAB\R2014b\toolbox\cobra'));
%addpath(genpath('..\cobra'));
%initCobraToolbox
load('non_dominated.mat');

V = numel(geni);
M = 3; %number of objectives

eps = 0.1;

%% COMPUTE PROFILES OF THE NONDOMINATED POINTS
% n_population = non_dominated(:,M+2);
% n_individual = non_dominated(:,M+3);
% 
% profiles_nondom = zeros(numel(n_population),V);
% for i = 1:numel(n_population)
%     file_name = ['populations_recon2_merged_bio_IDH1-2/population' num2str(n_population(i)) '.mat'];
%     file_pop = load(file_name);
%     chromosome = file_pop.chromosome;
%     profiles_nondom(i,:) = chromosome(n_individual(i),1:V);
% end
load('profiles_nondom.mat'); %if the profiles_nondom have been already computed before

ix_high_biomass = find(-non_dominated(:,1) > 1.5);
ix_high_oncoflux = find(-non_dominated(:,2) > 2000);
ix_low_distance = find(non_dominated(:,3) < 100);

ix_of_interest = intersect(ix_high_biomass,intersect(ix_high_oncoflux,ix_low_distance));

profiles_nondom = profiles_nondom(ix_of_interest,:);

%%

robustness_nondom_plus = zeros(numel(ix_of_interest),1);
robustness_nondom_minus = zeros(numel(ix_of_interest),1);
array_out_eps_plus = cell(numel(ix_of_interest),1);
array_out_eps_minus = cell(numel(ix_of_interest),1);

for i = 1:numel(ix_of_interest) %loop over all the nondominated individuals
    clc;
    i
    f = fopen('parfor_robustness_progress.txt', 'w');
    fprintf(f, [num2str(i) '\n']);
    fclose(f);

    profile = profiles_nondom(i,:);
    init_out = [-1 -1 1].*evaluate_objective(profile,M,V,fbarecon,geni,reaction_expression)   % we let FBA maximise oncoflux so we are sure that the oncoflux we get is the maximum possible with that particular gene expression profile
             
    out_eps_minus = cell(numel(geni),1);
    out_eps_plus = cell(numel(geni),1);
    difference_hypervolume_plus = zeros(numel(profile),1);
    difference_hypervolume_minus = zeros(numel(profile),1);
    
    for j=1:numel(geni)
        eps_profile_plus = profile; %returns back to the original profile (needed because the previous step modified it)
        eps_profile_plus(j) = profile(j) + eps;
        out_eps_plus{j} = [-1 -1 1].*evaluate_objective(eps_profile_plus,M,V,fbarecon,geni,reaction_expression)   % we let FBA maximise oncoflux so we are sure that the oncoflux we get is the maximum possible with that particular gene expression profile
         
        eps_profile_minus = profile; %returns back to the original profile (needed because the previous step modified it)
        eps_profile_minus(j) = max(profile(j) - eps,0);  %ensure the profile never goes negative
        out_eps_minus{j} = [-1 -1 1].*evaluate_objective(eps_profile_minus,M,V,fbarecon,geni,reaction_expression)   % we let FBA maximise oncoflux so we are sure that the oncoflux we get is the maximum possible with that particular gene expression profile
        
        difference_hypervolume_plus(j) = prod(abs(out_eps_plus{j} - init_out));
        difference_hypervolume_minus(j) = prod(abs(out_eps_minus{j} - init_out));
    end
    
   array_out_eps_plus{i} = out_eps_plus;
   array_out_eps_minus{i} = out_eps_minus;
   robustness_nondom_plus(i) = max(difference_hypervolume_plus);
   robustness_nondom_minus(i) = max(difference_hypervolume_minus);
end

%save('robustness_nondominated.mat',robustness_nondom);

    