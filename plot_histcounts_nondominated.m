load('non_dominated_chromosomes.mat');
load('others_chromosomes.mat');
load('others.mat')
load('non_dominated.mat');
load('geni_names.mat')
load('geni.mat')
n_bins = 5000;

all_chromosomes = [non_dominated_chromosomes; others_chromosomes];
all_solutions = [non_dominated; others];
%all_solutions = [non_dominated];
%%
all_biomass_values = all_solutions(:,1);


lim_inf = 1e-10;
%lim_sup = mean(all_biomass_values) + 2*std(all_biomass_values);
lim_sup = max(all_biomass_values) - 1e-10;

ix_low_biomass = find(all_biomass_values < lim_inf);
ix_mid_biomass = find(all_biomass_values >= lim_inf & all_biomass_values <= lim_sup);
ix_high_biomass = find(all_biomass_values > lim_sup);


%% using non_dominated.mat points


average_genes_low_biomass = mean(all_chromosomes(ix_low_biomass,:));
average_genes_mid_biomass = mean(all_chromosomes(ix_mid_biomass,:));
average_genes_high_biomass = mean(all_chromosomes(ix_high_biomass,:));
std_genes_low_biomass = std(all_chromosomes(ix_low_biomass,:));
std_genes_mid_biomass = std(all_chromosomes(ix_mid_biomass,:));
std_genes_high_biomass = std(all_chromosomes(ix_high_biomass,:));

low_biomass_chromosomes = all_chromosomes(ix_low_biomass,:);
mid_biomass_chromosomes = all_chromosomes(ix_mid_biomass,:);
high_biomass_chromosomes = all_chromosomes(ix_high_biomass,:);

tab_low_biomass = [num2cell((1:numel(geni))') geni geni_names num2cell(average_genes_low_biomass') num2cell(std_genes_low_biomass')];
tab_mid_biomass = [num2cell((1:numel(geni))') geni geni_names num2cell(average_genes_mid_biomass') num2cell(std_genes_mid_biomass')];
tab_high_biomass = [num2cell((1:numel(geni))') geni geni_names num2cell(average_genes_high_biomass') num2cell(std_genes_high_biomass')];


sorted_tab_low = sortrows(tab_low_biomass,-4);
sorted_tab_mid = sortrows(tab_mid_biomass,-4);
sorted_tab_high = sortrows(tab_high_biomass,-4);

histc1 = histcounts(low_biomass_chromosomes,n_bins);
histc2 = histcounts(mid_biomass_chromosomes,n_bins);
histc3 = histcounts(high_biomass_chromosomes,n_bins);



clear A
A.Low_Biomass = low_biomass_chromosomes;
A.Mid_Biomass = mid_biomass_chromosomes;
A.High_Biomass = high_biomass_chromosomes;


average_genes_low_biomass = mean(all_chromosomes(ix_low_biomass,:));
average_genes_mid_biomass = mean(all_chromosomes(ix_mid_biomass,:));
average_genes_high_biomass = mean(all_chromosomes(ix_high_biomass,:));
std_genes_low_biomass = std(all_chromosomes(ix_low_biomass,:));
std_genes_mid_biomass = std(all_chromosomes(ix_mid_biomass,:));
std_genes_high_biomass = std(all_chromosomes(ix_high_biomass,:));



figure
nhist(A,'linewidth',1,'box','color','sequential','binfactor',1,'proportion','maxbins',20,'samebins','smooth','location','East','xlabel','Gene expression average fold-change','ylabel','Frequency');
text(2.1,0.12,'p-value = 0.034','FontSize',10);

%%
% %% plot automatically genone-wide average expression
% figure
% bplot(low_biomass_chromosomes,'linewidth',0.3);
% title('Low Biomass - gene-wise distribution of gene expression across Pareto points');
% figure
% bplot(mid_biomass_chromosomes,'linewidth',0.3);
% title('Mid Biomass - gene-wise distribution of gene expression across Pareto points');
% figure
% bplot(high_biomass_chromosomes,'linewidth',0.3);
% title('High Biomass - gene-wise distribution of gene expression across Pareto points');
% 
% figure
% bplot(mid_biomass_chromosomes(2,:)')


% 
% %% using others.mat points (we have a better coverage of the objective space)
% ix_suboptimal = find(others(:,end-1) >= 300); %we take only the "others" points from the 300th generation onwards, because we want these points to be a result of NSGA in any case
% others_suboptimal = others(ix_suboptimal,:); 
% others_suboptimal_chromosomes = others_chromosomes(ix_suboptimal,:); 
% 
% n_points = 10;
% interv = linspace(-1e-5,1.7,n_points);
% for i = 1:n_points-1 % we have n_points-1  intervals
%     biomass_ixs{i} = find((others_suboptimal(:,1) >= interv(i)) & (others_suboptimal(:,1) <= interv(i+1)));
%     A_others.(['interval_' num2str(i)]) = others_suboptimal_chromosomes(biomass_ixs{i},:);
% end
