addpath(genpath('C:\Program Files\MATLAB\R2014b\toolbox\stats\stats'));
load('non_dominated.mat');
load('geni.mat');
load('geni_names.mat');
V = numel(geni);
M = 3; %number of objectives


%% COMPUTE PROFILES OF THE NONDOMINATED POINTS
%n_population = non_dominated(:,M+2);
%n_individual = non_dominated(:,M+3);
%
% profiles_nondom = zeros(numel(n_population),V);
% for i = 1:numel(n_population)
%     file_name = ['populations_recon2_merged_bio_IDH1-2/population' num2str(n_population(i)) '.mat'];
%     file_pop = load(file_name);
%     chromosome = file_pop.chromosome;
%     profiles_nondom(i,:) = chromosome(n_individual(i),1:V);
% end
load('profiles_nondom.mat'); %if the profiles_nondom have been already computed before
%%

ix_high_biomass = find(-non_dominated(:,1) > 1.5);
ix_high_oncoflux = find(-non_dominated(:,2) > 2000);
ix_low_distance = find(non_dominated(:,3) < 100);

ix_of_interest = intersect(ix_high_biomass,intersect(ix_high_oncoflux,ix_low_distance));

profiles = profiles_nondom(ix_of_interest,:);
genes_vs_profiles = profiles'; %we need to use profiles' and not profiles, because otherwise it would compute the correlation (and all the following measures) between profiles along all the genes, while we want the correlation between genes along the profiles
 

dist_correlation_vector = pdist_corr_1(genes_vs_profiles, 'correlation'); %pdist_corr_1 is my modified version of pdist, in which we use x - 1 instead of x - mean(x) all over the place in the definition of the correlation distance here  http://www.mathworks.co.uk/help/stats/pdist.html?refresh=true
dist_correlation_matrix = squareform(dist_correlation_vector); 
 
clusterTree = linkage(dist_correlation_vector, 'average');







%% *******************************************************************************************************************************************************************************
%% try k-means clustering with a few number of clusters: the result was that the max(silhouette) was with k=30, so we will continue with k=30
% 
% mean_silhouette = zeros(1,50);
% for NoClust = 2:50
%     [cidx, ctrs] = kmeans_corr_1(genes_vs_profiles, num_clusters, 'dist','corr', 'rep',5, 'disp','final'); %kmeans_corr_1 is my modified version of pdist, in which we use x - 1 instead of x - mean(x) all over the place in the definition of the correlation distance here  http://www.mathworks.co.uk/help/stats/pdist.html?refresh=true
%     %figure
%
%     % We now compute silhouette values. Large silhouette value, greater than 0.6, indicating that the cluster is somewhat separated from neighboring clusters. However, the third cluster contains many points with low silhouette values, and the first contains a few points with negative values, indicating that those two clusters are not well separated.
%     % The average s(i) over all data of a cluster is a measure of how tightly grouped all the data in the cluster are. 
%     figure;
%     [silh5,h] = silhouette(genes_vs_profiles,cidx,'corr');
%     mean_silhouette(NoClust) = mean(silh5);
%     h = gca;
%     h.Children.EdgeColor = [.8 .8 1];
%     xlabel 'Silhouette Value';
%     ylabel 'Cluster';
%     
% end
% 
% plot(mean_silhouette);

%%
num_clusters = 30;


%% k-means clustering

[cidx, ctrs] = kmeans_corr_1(genes_vs_profiles, num_clusters, 'dist','corr', 'rep',5, 'disp','final'); %kmeans_corr_1 is my modified version of pdist, in which we use x - 1 instead of x - mean(x) all over the place in the definition of the correlation distance here  http://www.mathworks.co.uk/help/stats/pdist.html?refresh=true
figure
for c = 1:num_clusters
    subplot(5,ceil(num_clusters/5),c);
    plot(genes_vs_profiles((cidx == c),:)');
    %axis tight
end
suptitle('K-Means Clustering of Genes');

% We now compute silhouette values. Large silhouette value, greater than 0.6, indicating that the cluster is somewhat separated from neighboring clusters. However, the third cluster contains many points with low silhouette values, and the first contains a few points with negative values, indicating that those two clusters are not well separated.
% The average s(i) over all data of a cluster is a measure of how tightly grouped all the data in the cluster are.
figure;
[silh5,h] = silhouette_corr_1(genes_vs_profiles,cidx,'corr');
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value';
ylabel 'Cluster';


% PLOT AVERAGE IN EVERY CLUSTER
figure
for c = 1:num_clusters
    subplot(5,ceil(num_clusters/5),c);
    plot(ctrs(c,:)');
    %axis tight
    axis off    % turn off the axis
end
suptitle('K-Means Clustering of Profiles');
%     
%clustergram(genes_vs_profiles(:,2:end));










%% *******************************************************************************************************************************************************************************
%% try hierarchical clustering with a few number of clusters: the result was that the max(silhouette) was with k=12, so we will continue with k=12
% 
% mean_silhouette = zeros(1,50);
% for NoClust = 2:50
%     clusters = cluster(clusterTree, 'maxclust', NoClust);
% 
%     % We now compute silhouette values. Large silhouette value, greater than 0.6, indicating that the cluster is somewhat separated from neighboring clusters. However, the third cluster contains many points with low silhouette values, and the first contains a few points with negative values, indicating that those two clusters are not well separated.
%     % The average s(i) over all data of a cluster is a measure of how tightly grouped all the data in the cluster are. 
%     figure;
%     [silh5,h] = silhouette_corr_1(genes_vs_profiles,clusters,'corr');
%     mean_silhouette(NoClust) = mean(silh5);
%     h = gca;
%     h.Children.EdgeColor = [.8 .8 1];
%     xlabel 'Silhouette Value';
%     ylabel 'Cluster';
%     
% end

%%
num_clusters = 12;

%% Hierarchical clustering


clusters = cluster(clusterTree, 'maxclust', num_clusters);
numel_clusters = accumarray(clusters,1);   %cardinality of clusters, i.e. number of genes (elements) in each cluster
gene_names_in_cluster = cell(num_clusters,1);

figure
for c = 1:num_clusters
    subplot(3,ceil(num_clusters/3),c);
    plot(genes_vs_profiles((clusters == c),:)');
    genes_in_cluster = find(clusters==c);
    gene_names_in_cluster{c} = strjoin(geni_names(genes_in_cluster)',', ');
    title({   [num2str(numel_clusters(c)) ':  ' strjoin(geni_names(genes_in_cluster(1:min(numel_clusters(c),3)))',', ')]     ,     strjoin(geni_names(genes_in_cluster(4:min(numel_clusters(c),6)))',', ')     }, 'FontSize',7);
    %axis tight
end
suptitle('Hierarchical Clustering of Genes');



% We now compute silhouette values. Large silhouette value, greater than 0.6, indicating that the cluster is somewhat separated from neighboring clusters. However, the third cluster contains many points with low silhouette values, and the first contains a few points with negative values, indicating that those two clusters are not well separated.
% The average s(i) over all data of a cluster is a measure of how tightly grouped all the data in the cluster are.
figure;
[silh5,h] = silhouette_corr_1(genes_vs_profiles,clusters,'corr');
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value';
ylabel 'Cluster';






%% PLOT DISTANCE BETWEN GENES THROUG MULTI-DIMENSIONAL SCALING APPLIED TO THE CORRELATIO MATRIX TO MAKE It BECOME A DISTANCE MATRIX

[Y,stress] = mdscale(dist_correlation_vector,2,'criterion','metricstress');
figure
plot(Y(:,1),Y(:,2),'.','LineWidth',2);
C = clusters; %colour according to hierarchical clustering
colormap(jet(256))
scatter(Y(:,1),Y(:,2),200,C,'.');
title('Hierarchical clustering (k=12)');
colorbar;
%gname(geni);


figure
C = cidx;  %colour according to kmeans clustering
colormap(jet(256))
scatter(Y(:,1),Y(:,2),200,C,'.');
title('K-Means Clustering (k=30)');
colorbar;


load('EDGE_results.mat');
figure
low_edge = find(edge_bio < 0.1);
high_edge = find(edge_bio >= 0.1);
C = edge_bio(low_edge);  %colour according to the low edge score
colormap(jet(256))
scatter(Y(low_edge,1),Y(low_edge,2),200,C,'.');
title('EDGE score');
colorbar;
hold on
scatter(Y(high_edge,1),Y(high_edge,2),200,'black','X');
%gname(geni_names);
disp('Genes with highest EDGE:');
disp(geni_names(high_edge))



cluster_edge_table = [geni geni_names num2cell(edge_bio) num2cell(clusters) num2cell(cidx) ];
cluster_edge_table = sortrows(cluster_edge_table,3);
distance_from_zero = norm(Y);

low_distance_genes = find(distance_from_zero < 0.03);
for i = 1:numel(low_distance_genes)
    occurrences_low_distance_genes{i} = find(~cellfun('isempty',strfind(fbarecon.grRules,geni{low_distance_genes(i)})));
end


%exportfig(gcf, 'figura.pdf', 'color', 'cmyk', 'Width', '1000', 'Height', '600', 'FontMode', 'scaled', 'FontSize', '1' );