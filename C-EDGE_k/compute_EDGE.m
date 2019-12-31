load('geni.mat')
load('geni_names.mat');
load('reaction_expression.mat')
load('recon2_merged_bio_IDH1-2.mat')
%addpath(genpath('C:\Program Files\MATLAB\R2014b\toolbox\cobra'));
addpath(genpath('..\cobra'));
%initCobraToolbox

M=3; %number of objectives
V = numel(geni);
%load('population124.mat');  %from non_dominated, we have chosen the 6th individual in population 124
%x = chromosome(6,1:V);
x = ones(1,V); %we start from the all-one configuration, than we impose only the expression of that gene to be epsilon (exactly as in EDGE), while all the others are kept 1.

init_out = evaluate_objective(x,M,V,fbarecon,geni,reaction_expression)   % we let FBA maximise oncoflux so we are sure that the oncoflux we get is the maximum possible with that particular gene expression profile
init_out(1) = -init_out(1);
init_out(2) = -init_out(2);

eps = 0.1;


% %% single gene epsilon perturbation (useless because we find the same result on the diagonal of the pairwise correlation matrix)
% edge_bio = zeros(numel(geni),1);
% edge_onco = zeros(numel(geni),1);
% edge_bio_onco_hyperv = zeros(numel(geni),1);
% out_eps_bio = zeros(numel(geni),1);
% out_eps_oncoflux = zeros(numel(geni),1);
% out_ko_bio = zeros(numel(geni),1);
% out_ko_oncoflux = zeros(numel(geni),1);
% 
% parfor i = 1:numel(geni)
%     %eps = zeros(1,numel(geni));
%     %eps(i) = 0.1;
%     
%     x_eps = x; x_eps(i) = eps;
%     out_eps = evaluate_objective(x_eps,M,V,fbarecon,geni,reaction_expression)   % we let FBA maximise oncoflux so we are sure that the oncoflux we get is the maximum possible with that particular gene expression profile
%     out_eps(1) = -out_eps(1);
%     out_eps(2) = -out_eps(2);
% 
%     x_ko = x; x_ko(i)=0;
%     out_ko = evaluate_objective(x_ko,M,V,fbarecon,geni,reaction_expression)   % we let FBA maximise oncoflux so we are sure that the oncoflux we get is the maximum possible with that particular gene expression profile
%     out_ko(1) = -out_ko(1);
%     out_ko(2) = -out_ko(2);
% 
%     out_eps_bio(i) = out_eps(1);
%     out_eps_oncoflux(i) = out_eps(2);
%     out_ko_bio(i) = out_ko(1);
%     out_ko_oncoflux(i) = out_ko(2);
%     
%     edge_bio(i) = out_eps(1)-out_ko(1);
%     edge_onco(i) = out_eps(2)-out_ko(2);
%     edge_bio_onco_hyperv(i) = out_eps(1)*out_eps(2) - out_ko(1)*out_ko(2);
% 
%     
% %     disp([edge_bio(i) edge_onco(i) edge_bio_onco_hyperv(i)]); 
%     
% end



%% pairwise epsilon perturbation

starting_row = input('Insert starting row: ','s');
starting_row = str2num(starting_row);

if starting_row > 1 
%     load('edge_bio_pairwise.mat');
%     load('edge_oncoflux_pairwise.mat');
%     load('edge_bio_onco_hyperv_pairwise.mat');
    load('out_eps_bio_pairwise.mat');
    load('out_ko_bio_pairwise.mat');
    load('out_eps_oncoflux_pairwise.mat');
    load('out_ko_oncoflux_pairwise.mat');
else
    edge_bio_pairwise = zeros(numel(geni),numel(geni));
    edge_oncoflux_pairwise = zeros(numel(geni),numel(geni));
    edge_bio_onco_hyperv_pairwise = zeros(numel(geni),numel(geni));
    out_eps_bio_pairwise = zeros(numel(geni),numel(geni));
    out_ko_bio_pairwise = zeros(numel(geni),numel(geni));
    out_eps_oncoflux_pairwise = zeros(numel(geni),numel(geni));
    out_ko_oncoflux_pairwise = zeros(numel(geni),numel(geni));
end

for  i = starting_row:numel(geni)
    clc;
    i
    f = fopen('parfor_progress.txt', 'w');
    fprintf(f, [num2str(i) '\n']);
    fclose(f);
%     aux_edge_bio = zeros(1,numel(geni));
%     aux_edge_oncoflux = zeros(1,numel(geni));
%     aux_bio_onco_hyperv = zeros(1,numel(geni));
    aux_out_eps_bio = zeros(1,numel(geni));
    aux_out_eps_oncoflux = zeros(1,numel(geni));
    aux_out_ko_bio = zeros(1,numel(geni));
    aux_out_ko_oncoflux = zeros(1,numel(geni));
    
    parfor j = i:numel(geni)
        disp(j);
        x_eps = x; x_eps([i j]) = eps;
        out_eps = evaluate_objective(x_eps,M,V,fbarecon,geni,reaction_expression);   % we let FBA maximise oncoflux so we are sure that the oncoflux we get is the maximum possible with that particular gene expression profile
        aux_out_eps_bio(j) = -out_eps(1);
        aux_out_eps_oncoflux(j) = -out_eps(2);
    
        x_ko = x; x_ko([i j])=0;
        out_ko = evaluate_objective(x_ko,M,V,fbarecon,geni,reaction_expression);   % we let FBA maximise oncoflux so we are sure that the oncoflux we get is the maximum possible with that particular gene expression profile
        aux_out_ko_bio(j) = -out_ko(1);
        aux_out_ko_oncoflux(j) = -out_ko(2);
    
%         aux_edge_bio(j) = out_eps(1)-out_ko(1);
%         aux_edge_oncoflux(j) = out_eps(2)-out_ko(2);
%         aux_bio_onco_hyperv(j) = out_eps(1)*out_eps(2) - out_ko(1)*out_ko(2);
%        disp([aux_edge_bio(j) aux_edge_oncoflux(j) aux_bio_onco_hyperv(j)]); 
        
    end
    
   
    out_eps_bio_pairwise(i,:) = aux_out_eps_bio;
    out_ko_bio_pairwise(i,:) = aux_out_ko_bio;
    out_eps_oncoflux_pairwise(i,:) = aux_out_eps_oncoflux;    
    out_ko_oncoflux_pairwise(i,:) = aux_out_ko_oncoflux;
    
    save('out_eps_bio_pairwise.mat','out_eps_bio_pairwise');
    save('out_ko_bio_pairwise.mat','out_ko_bio_pairwise');
    save('out_eps_oncoflux_pairwise.mat','out_eps_oncoflux_pairwise');
    save('out_ko_oncoflux_pairwise.mat','out_ko_oncoflux_pairwise');
end

edge_bio_pairwise = out_eps_bio_pairwise - out_ko_bio_pairwise;
edge_oncoflux_pairwise = out_eps_oncoflux_pairwise - out_ko_oncoflux_pairwise;
edge_bio_onco_hyperv_pairwise = out_eps_bio_pairwise.*out_eps_oncoflux_pairwise - out_ko_bio_pairwise.*out_ko_oncoflux_pairwise;

    
% we need now to put in the lower triangular part the same nubers calculated for the upper triangular part
edge_bio_pairwise = edge_bio_pairwise + edge_bio_pairwise.' - eye(size(edge_bio_pairwise)).*edge_bio_pairwise;
edge_oncoflux_pairwise = edge_oncoflux_pairwise + edge_oncoflux_pairwise.' - eye(size(edge_oncoflux_pairwise)).*edge_oncoflux_pairwise;
edge_bio_onco_hyperv_pairwise = edge_bio_onco_hyperv_pairwise + edge_bio_onco_hyperv_pairwise.' - eye(size(edge_bio_onco_hyperv_pairwise)).*edge_bio_onco_hyperv_pairwise;
out_eps_bio_pairwise = out_eps_bio_pairwise + out_eps_bio_pairwise.' - eye(size(out_eps_bio_pairwise)).*out_eps_bio_pairwise;
out_ko_bio_pairwise = out_ko_bio_pairwise + out_ko_bio_pairwise.' - eye(size(out_ko_bio_pairwise)).*out_ko_bio_pairwise;
out_eps_oncoflux_pairwise = out_eps_oncoflux_pairwise + out_eps_oncoflux_pairwise.' - eye(size(out_eps_oncoflux_pairwise)).*out_eps_oncoflux_pairwise;
out_ko_oncoflux_pairwise = out_ko_oncoflux_pairwise + out_ko_oncoflux_pairwise.' - eye(size(out_ko_oncoflux_pairwise)).*out_ko_oncoflux_pairwise;

save('edge_bio_pairwise.mat','edge_bio_pairwise');
save('edge_oncoflux_pairwise.mat','edge_oncoflux_pairwise');
save('edge_bio_onco_hyperv_pairwise.mat','edge_bio_onco_hyperv_pairwise');




%% CHECK OCCURRENCES OF HIGH-EDGE GENES
high_edge_genes = find(edge_bio > 1);
for i = 1:numel(high_edge_genes)
    occurrences_high_edge_genes{i} = find(~cellfun('isempty',strfind(fbarecon.grRules,geni{high_edge_genes(i)})));
end

%% PLOT RESULTS
% SEE statistics_on_genes for further plotting results combining EDGE with
% clustering

high_edge = find(edge_bio>1.5);
geni_names(high_edge)

edge_table = [geni geni_names num2cell(edge_bio)];
edge_table = sortrows(edge_table,3);


% PLOT EDGE ONLY TAKING INTO ACCOUNT LOW VALUES
M1 = edge_bio_pairwise;
figure
imagesc(M1);
title('EDGE')
[rows,cols] = ind2sub(size(M1),find(M1>1.5));
cb = colorbar;
set(cb,'XTickLabel',cellstr(sprintf('%.10f\n',yt))); %for Matlab r2014b. See plot_and_export_color for the alternative for older versions of Matlab


M1_low = edge_bio_pairwise;
%M1_low(find(M1_low>0.1))=0;
figure
imagesc(M1_low);
colorbar;
set(cb,'XTickLabel',cellstr(sprintf('%.10f\n',yt))); %for Matlab r2014b. See plot_and_export_color for the alternative for older versions of Matlab
interesting_part = find(M1_low<0.1);
caxis([min(M1_low(interesting_part)) max(M1_low(interesting_part))]); %Chang colormap scaling: maps to the min all the values less than the min, and to the max all the values greater than the max.


M2 = edge_oncoflux_pairwise;
figure
interesting_part = find(M2>0.5);
non_interesting_part = setdiff(1:(size(M2,2)^2),interesting_part); %set to NaN all the non interesting indices of the matrix, so they will not be plotted
imagesc(M2);
[rows,cols] = ind2sub(size(M2),interesting_part);
colorbar;
set(cb,'XTickLabel',cellstr(sprintf('%.10f\n',yt))); %for Matlab r2014b. See plot_and_export_color for the alternative for older versions of Matlab


M2_low = edge_oncoflux_pairwise;
interesting_part = find(M2_low<0.1);
non_interesting_part = setdiff(1:(size(M2_low,2)^2),interesting_part); %set to NaN all the non interesting indices of the matrix, so they will not be plotted
figure
imagesc(M2_low);
colorbar;
set(cb,'XTickLabel',cellstr(sprintf('%.10f\n',yt))); %for Matlab r2014b. See plot_and_export_color for the alternative for older versions of Matlab
caxis([min(M2_low(interesting_part)) max(M2_low(interesting_part))]); %Chang colormap scaling: maps to the min all the values less than the min, and to the max all the values greater than the max.




%% COMPUTE PAIRS OF GENES THAT ARE TOXIC TOGETHER BUT NOT TOXIC ALONE, AND FLAG THESE PAIRS IN A TABLE
list=[];
value=[];
table_flag = zeros(numel(geni),numel(geni));

[rows_pairwise_toxic,cols_pairwise_toxic] = ind2sub(size(edge_bio_pairwise),find(edge_bio_pairwise < 0));
[~,idx] = unique(sort([rows_pairwise_toxic cols_pairwise_toxic],2),'rows','stable'); %remove the useless permutation in the pair
rows_pairwise_toxic = rows_pairwise_toxic(idx);
cols_pairwise_toxic = cols_pairwise_toxic(idx);

for i = 1:numel(rows_pairwise_toxic)
    if (edge_bio(rows_pairwise_toxic(i)) >= 0) && (edge_bio(cols_pairwise_toxic(i)) >= 0) %this means they are both non toxic alone, although they are toxic together!
        list(end+1) = i;
        value(end+1) = abs(edge_bio_pairwise(rows_pairwise_toxic(i),cols_pairwise_toxic(i)) - edge_bio(rows_pairwise_toxic(i))) + abs(edge_bio_pairwise(rows_pairwise_toxic(i),cols_pairwise_toxic(i)) - edge_bio(cols_pairwise_toxic(i)));
        table_flag(rows_pairwise_toxic(i), cols_pairwise_toxic(i)) = 1; %flag the fact that that pair is intresting!
        table_flag(cols_pairwise_toxic(i), rows_pairwise_toxic(i)) = 1; %flag the fact that that pair is intresting!
    end
end
interesting_pairs = [geni_names(rows_pairwise_toxic(list)) geni_names(cols_pairwise_toxic(list))];
%interesting_pairs = [geni(rows_pairwise_toxic(list)) geni(cols_pairwise_toxic(list))];


%% COMPUTE PAIRS OF GENES THAT ARE NON-TOXIC TOGETHER BUT TOXIC ALONE, AND FLAG THESE PAIRS IN A TABLE
list2=[];
value2=[];
table_flag2 = zeros(numel(geni),numel(geni));

[rows_pairwise_toxic2,cols_pairwise_toxic2] = ind2sub(size(edge_bio_pairwise),find(edge_bio_pairwise > 0));
[~,idx] = unique(sort([rows_pairwise_toxic2 cols_pairwise_toxic2],2),'rows','stable'); %remove the useless permutation in the pair
rows_pairwise_toxic2 = rows_pairwise_toxic2(idx);
cols_pairwise_toxic2 = cols_pairwise_toxic2(idx);

for i = 1:numel(rows_pairwise_toxic2)
    if (edge_bio(rows_pairwise_toxic2(i)) <= 0) && (edge_bio(cols_pairwise_toxic2(i)) <= 0) %this means they are both non toxic alone, although they are toxic together!
        list2(end+1) = i;
        value2(end+1) = abs(edge_bio_pairwise(rows_pairwise_toxic2(i),cols_pairwise_toxic2(i)) - edge_bio(rows_pairwise_toxic2(i))) + abs(edge_bio_pairwise(rows_pairwise_toxic2(i),cols_pairwise_toxic2(i)) - edge_bio(cols_pairwise_toxic2(i)));
        table_flag2(rows_pairwise_toxic2(i), cols_pairwise_toxic2(i)) = 1; %flag the fact that that pair is intresting!
        table_flag2(cols_pairwise_toxic2(i), rows_pairwise_toxic2(i)) = 1; %flag the fact that that pair is intresting!
    end
end


interesting_pairs2 = [geni_names(rows_pairwise_toxic2(list2)) geni_names(cols_pairwise_toxic2(list2))];
%interesting_pairs2 = [geni(rows_pairwise_toxic2(list2)) geni(cols_pairwise_toxic2(list2))];


%% EVALUATE A VERY INTERESTING PARTICULAR CASE FOR THE PREVIOUS ANALYSIS
i = 627;
j = 628;
x_eps = x; 
range = 0:0.0001:0.0249; 
count=0;
out_eps_var = zeros(numel(range),M);
for eps = range
    eps
    count = count+1;
    x_eps([i j]) = eps;
    out_eps_var(count,:) = evaluate_objective(x_eps,M,V,fbarecon,geni,reaction_expression);   % we let FBA maximise oncoflux so we are sure that the oncoflux we get is the maximum possible with that particular gene expression profile
end
%out_eps_var(:,1:2) = -out_eps_var(:,1:2);
plot(range,out_eps_var(:,1));
xlabel('\epsilon');
ylabel('Biomass = EDGE-pairwise');
title('Genes (8560.2, 8560.1) yield 0 biomass only if both KO, and maximum biomass if only one of them is KO');
