addpath('C:\Program Files\SBML\libSBML-5.1.0b0-libxml2-x64\bindings\matlab\matlab');
addpath(genpath('C:\Program Files\MATLAB\R2014b\toolbox\cobra'));
initCobraToolbox



load('profiles_nondom.mat');
load('geni.mat');
load('recon2_merged_bio_IDH1-2.mat');
%cell_fbarecon = readCbModel('gastric_cancer.xml');
cell_fbarecon = load('breast_normal.mat');




%%%%NO!! la nuova procedura Partire da recon2 curated e togliere le reazioni non presenti nel
%%%%nuovo modello. Per toglierle, mettere a zero lower bound e upper bound
%%%%e basta, e poi testare con i profiles_nondom che ho gia'

% new_profile = ones(1,numel(model.genes));
% 
% for i=1:numel(profiles_nondom)
%     profile = profiles_nondom(i,:);
%     for i = 1:numel(profile)  %for every gene expression in the profiles_nondom, we search for the position of the gene and we report it in the new model
%         new_gene_index = find(ismember(model.genes,geni{i});    %the genes of the original Recon model that will not be found in the new model, will be left at 1 in the new model
%         new_profile(new_gene_index) = profile(i);
%     end
% end