% extracts all the subsets found by NSGA and computes

load('non_dominated.mat');
load('others.mat');
load('recon2_merged_bio_IDH1-2.mat');   %model with objectives and positions changed. Even if it's already in the workspace, must be reloaded every time this script starts because this script changes lower and upper bounds (multiplying something X old_bounds, so as to compute the output of FBA, so when a new output is required the orginal model must be reloaded because we need the original bounds
load('reaction_expression.mat');
load('geni.mat');
load('geni_names.mat');

M=2; %number of objectives
V = numel(geni);    %the length of the input individuals (without ranking, crowding distance and outputs) is equal to the number of reactions minus 2 that is the number of biomass and synthetic objective in the model (the final 2 reactions)

chroms = [non_dominated(:,end); others(:,end)];
pops = [non_dominated(:,end-1); others(:,end-1)];

eps=0.01;

subset = cell(numel(pops),1);
interesting_subsets = cell(numel(pops),1);
interesting_subsets_names = cell(numel(pops),1);

EDGE_k_interesting = cell(numel(pops),1);
EDGE_kminus1_interesting = cell(numel(pops),1);
%EDGE_diff_interesting = cell(numel(pops),1);
EDGE_prod_diff_interesting = cell(numel(pops),1);

parfor i=1:numel(pops)
    n_pop = pops(i);
    n_chrom = chroms(i);
    sol=['population' num2str(n_pop) '.mat'];
    cd populations_recon2_merged_bio_IDH1-2
    chromosome = load(sol);
    chromosome=chromosome.chromosome;
    cd ..
    
    x = chromosome(n_chrom,1:V);
    subset_selected = find(x==1);
    
    k = numel(subset_selected);
    child_eps = ones(1,V);
    child_eps(subset_selected) = eps;
    child_ko = ones(1,V);
    child_ko(subset_selected) = 0;
    out_k_eps = [-1 -1 1].*evaluate_objective_EDGE(child_eps, M, V, fbarecon, geni, reaction_expression);
    out_k_ko= [-1 -1 1].*evaluate_objective_EDGE(child_ko, M, V, fbarecon, geni, reaction_expression);
    EDGE_k = out_k_eps(1) - out_k_ko(1); % we only consider the biomass
    
    [indices, subsets] = cSubsets(k, k-1); %finds all the subsets of 1:k having k-1 elements, and we will compute their EDGE
    EDGE_kminus1 = zeros(size(subsets,1),1);
    for ixSubs = 1:size(subsets,1)
        child_eps = ones(1,V);
        child_eps(subset_selected(subsets(ixSubs,:))) = eps;
        child_ko = ones(1,V);
        child_ko(subset_selected(subsets(ixSubs,:))) = 0;
        out_kminus1_eps = [-1 -1 1].*evaluate_objective_EDGE(child_eps, M, V, fbarecon, geni, reaction_expression);
        out_kminus1_ko = [-1 -1 1].*evaluate_objective_EDGE(child_ko, M, V, fbarecon, geni, reaction_expression);
        EDGE_kminus1(ixSubs) = out_kminus1_eps(1) - out_kminus1_ko(1);
    end
    
    subset{i} = subset_selected;
    interesting_subsets{i} = geni(subset{i});
    interesting_subsets_names{i} = geni_names(subset{i});
    EDGE_k_interesting{i} = EDGE_k;
    EDGE_kminus1_interesting{i} = EDGE_kminus1';
    
    
    %     max_EDGE_kminus1 = max(EDGE_kminus1);
    %     min_EDGE_kminus1 = min(EDGE_kminus1);
    %
    %     if isempty(max_EDGE_kminus1) max_EDGE_kminus1=0; end     %prevents the case k=1 for which there are no subsets
    %
    %     EDGE_diff = abs(EDGE_k - max_EDGE_kminus1);
    %     EDGE_diff_interesting{i} = EDGE_diff;
    %     fprintf('|EDGE_k - max{EDGE_(k-1)}| = %.15f, k = %d\n',EDGE_diff, k);
    %
    %     EDGE_sumdiff = abs(EDGE_k - max_EDGE_kminus1)+abs(EDGE_k - min_EDGE_kminus1);
    %     EDGE_sumdiff_interesting{i} = EDGE_sumdiff;
    %     fprintf('|EDGE_k - max{EDGE_(k-1)}| + |EDGE_k - min{EDGE_(k-1)}| = %.15f, k = %d\n',EDGE_diff, k);
    
    EDGE_prod_diff = 1;
    for ixSubs = 1: k
        EDGE_prod_diff = EDGE_prod_diff * abs(EDGE_k - EDGE_kminus1(ixSubs));
    end
    EDGE_prod_diff_interesting{i} = EDGE_prod_diff;
end

cd ..