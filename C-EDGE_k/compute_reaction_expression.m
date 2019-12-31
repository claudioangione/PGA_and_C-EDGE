%load genes_reactions.mat
%load geo_s_genes_reactions.mat


%load('geo_s_react.mat');
load('recon2_merged.mat');

genesets = fbarecon.grRules;

%the genesets are already written nicely and ordered by reaction, so there
%is no need to assign each geneset to a reaction, so there is no need to
%use the following code. We only need to convert "AND" into "and", "OR"
%into "or" because of how associate_genes_reactions.m is written

%% WARNING: we can remove brackets only in Recon because brackets are only used for gene names and not for associative rules between AND and OR. Strangely there is no associativity in the grRules!
%However, in general, to avoid problems I am substituting all the '((' and
%'))' (that are the useful brackets because the first '(' is part of the associative rules, while the second '(' is part of the gene name) with '[' and ']' to save them
% then I am removing all the '(' and ')', 
%and finally I am substituting back '[' and ']' with '(' and ')' respectvely. 
% This way, all the useless brackets part of the gene names are
%removed, while the useful brackets are kept
genesets = regexprep(genesets,'((','['); 
genesets = regexprep(genesets,'))',']'); 
genesets = regexprep(genesets,'(',''); 
genesets = regexprep(genesets,')',''); 
genesets = regexprep(genesets,'[','('); 
genesets = regexprep(genesets,']',')'); 

%%
genesets = regexprep(genesets,' AND ',' and '); 
genesets = regexprep(genesets,' OR ',' or '); 
%genesets = upper(genesets);  %upper('str') converts any lowercase characters in the string str to the corresponding uppercase characters and leaves all other characters unchanged. We need it because some genes in the genesets are YorXXX while in the fbamodel are YORXXX

% %Boolean_genesets = hacky_genes_reactions_conversion('geo_s_react');
%Boolean_genesets = fbamodel.grRules;  %fbamodel.grRules contains the genesets (also fbamodel.rules, but it contains only gene numbers and not their names)

% Boolean_genesets = strcat('(',Boolean_genesets,')');
% Boolean_genesets = strcat('(',Boolean_genesets,')');
% Boolean_genesets = regexprep(Boolean_genesets,', ',' and '); 
% Boolean_genesets = regexprep(Boolean_genesets,'_AND_',' and '); 
% %Boolean_genesets = regexprep(Boolean_genesets,'_OR_',' or '); 
% Boolean_genesets = regexprep(Boolean_genesets,'_OR_',') or (');  %brackets are necessary to force associative rules when these are not included in the original Boolean genesets
% 
% Boolean_genesets = regexprep(Boolean_genesets,'  ',' '); 
% Boolean_genesets = regexprep(Boolean_genesets,'  ',' '); 
% Boolean_genesets = regexprep(Boolean_genesets,'  ',' '); 
% 
% 

% full_G = full(fbamodel.G); %geneset-to-reaction matrix
% fbamodel.pts; %genesets
% 
% genesets = cell(size(full_G,2),1);   %we need to map the genesets to the reactions (all the genesets are associated to a reaction, but non every reaction in the model has a corresponding geneset (e.g. the exchange or sink reactions do not have genesets)
% genesets(:) = {''};
% for j = 1 : size(full_G,2)
%     i = find(full_G(:,j) == 1);  %find the geneset j responsible for reaction i (which means G(j,i)=1)
%     if ~isempty(i)
%         genesets(j) = Boolean_genesets(i);
%     end
% end


reaction_expression = cell(length(genesets),1);  %initializes an empty cell array, where each cell will be a string (it's the only way to create an array of strings in matlab)
reaction_expression(:) = {''};

parfor_progress(length(genesets)); %initialise
parfor i = 1:length(genesets)
    %i
    aux = associate_genes_reactions(genesets{i});
    reaction_expression{i} = aux;  
    %parfor_progress;
end

reaction_expression



geni = recon2.genes; %we take genes from recon2 because we have added some gene-to-reaction rules from recon2 that were not in reconNielsen

%% everything follows is useless because it's all in fbamodel.genes
% we also prepare the array of genes
%geni = gpr2genes(fbamodel.grRules);
%geni = unique(geni); %removes duplicates, we need unique array of genes, without duplicates
%reaction_expression = compute_geneset_expression('geo_s_react',hacky_genes_reactions_conversion('geo_s_react'));

