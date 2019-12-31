% % Compute controllability using our generalised EDGE starting from a point in the objective space, and start perturbing all its genes to see on average how much that point is moving in the objective space. Considering that the only way to avoid weighting is multiplying the coordinates, I think that the correct measure in that case is evaluating the difference between two points by computing the absolute product of the difference between the two coordinates
% % |(x1-x2)(y1-y2)|.
% % In general, it's the volume of the hypercube whose opposite vertices are our two points on the Pareto.
%


load('geni.mat')
load('reaction_expression.mat')
load('recon2_merged_bio_PHGDH.mat')
load('pos_genes_in_react_expr.mat');
load('ixs_geni_sorted_by_length.mat');
% addpath(genpath('C:\Program Files\MATLAB\R2014b\toolbox\cobra'));
% %addpath(genpath('..\cobra'));
% %initCobraToolbox

load('non_dominated.mat');
load('others.mat');
load('non_dominated_chromosomes.mat'); %if the profiles_nondom have been already computed before
load('others_chromosomes.mat');

V = numel(geni);
M = 2;
eps = 0.01;

all_chromosomes = [non_dominated_chromosomes; others_chromosomes];
all_points = [non_dominated; others];

all_points(:,2) = -all_points(:,2); %because we maximised for the second objective, therefore the second column has negative values (due to how NSGA-II works)

ix_of_interest = find(all_points(:,2) > mean(all_points(:,2))+2*std(all_points(:,2)));  %points whose phgdh is outlier, i.e. above mean+2*std_dev
chromosomes_of_interest = all_chromosomes(ix_of_interest,:);
points_of_interest = all_points(ix_of_interest,:);

init_out = points_of_interest(:,1:2); %we already have the output for the points-of_interest points, so we don't need the previous instruction


% %% CONTROLLABILITY ANALYSIS
% %%
%
% controllability_nondom_plus = zeros(numel(ix_of_interest),1);
% controllability_nondom_minus = zeros(numel(ix_of_interest),1);
% controllability_nondom = zeros(numel(ix_of_interest),1);
% array_out_eps_plus = cell(numel(ix_of_interest),1);
% array_out_eps_minus = cell(numel(ix_of_interest),1);
% 
% parfor i = 1:numel(ix_of_interest) %loop over all the "of interest" individuals
%     
%     i
%     ix = ix_of_interest(i);
%     
%     profile = all_chromosomes(ix,1:V);
%     init_out = [1 -1].*evaluate_objective(profile,M,V,fbarecon,geni,reaction_expression,pos_genes_in_react_expr, ixs_geni_sorted_by_length)   % put [-1 +1] exactly the opposite of the -1 and +1 put at the end of evaluate objective, in order to cancel their effect. We let FBA maximise oncoflux so we are sure that the oncoflux we get is the maximum possible with that particular gene expression profile
%     
%     
%     %     %% PERTURBATION ONE GENE AT A TIME
%     %     out_eps_minus = cell(numel(geni),1);
%     %     out_eps_plus = cell(numel(geni),1);
%     %     difference_hypervolume_plus = zeros(numel(profile),1);
%     %     difference_hypervolume_minus = zeros(numel(profile),1);
%     %
%     %     parfor j=1:numel(geni)
%     %         eps_profile_plus = profile; %returns back to the original profile (needed because the previous step modified it)
%     %         eps_profile_plus(j) = profile(j) + eps;
%     %         out_eps_plus{j} = [1 -1].*evaluate_objective(eps_profile_plus,M,V,fbarecon,geni,reaction_expression)
%     %
%     %         eps_profile_minus = profile; %returns back to the original profile (needed because the previous step modified it)
%     %         eps_profile_minus(j) = max(profile(j) - eps,0);  %ensure the profile never goes negative
%     %         out_eps_minus{j} = [1 -1].*evaluate_objective(eps_profile_minus,M,V,fbarecon,geni,reaction_expression)
%     %
%     %         difference_hypervolume_plus(j) = prod(abs(out_eps_plus{j} - init_out));
%     %         difference_hypervolume_minus(j) = prod(abs(out_eps_minus{j} - init_out));
%     %     end
%     %
%     %     array_out_eps_plus{i} = out_eps_plus;
%     %     array_out_eps_minus{i} = out_eps_minus;
%     %     controllability_nondom_plus(i) = max(difference_hypervolume_plus);
%     %     controllability_nondom_minus(i) = max(difference_hypervolume_minus);
%     %     controllability_nondom(i) = max(controllability_nondom_plus(i),controllability_nondom_minus(i))
%     
%     
%     %% PERTURBATION ALL GENES AT A TIME (GLOBAL)
%     eps_array = eps*ones(1,numel(profile));
%     eps_profile_plus = profile + eps_array;
%     out_eps_plus = [1 -1].*evaluate_objective(eps_profile_plus,M,V,fbarecon,geni,reaction_expression,pos_genes_in_react_expr, ixs_geni_sorted_by_length)
%     
%     eps_profile_minus = max(profile - eps_array,0);  %ensure the profile never goes negative
%     %    eps_profile_minus = bsxfun(@max,profile - eps_array,zeros(1,numel(profile)));
%     out_eps_minus = [1 -1].*evaluate_objective(eps_profile_minus,M,V,fbarecon,geni,reaction_expression,pos_genes_in_react_expr, ixs_geni_sorted_by_length)
%     
%     
%     
%     array_out_eps_plus{i} = out_eps_plus;
%     array_out_eps_minus{i} = out_eps_minus;
%     
%     difference_hypervolume_plus = prod(abs(array_out_eps_plus{i}  - init_out));
%     difference_hypervolume_minus = prod(abs(array_out_eps_minus{i} - init_out));
%     
%     % we finally compute the controllability: high controllability means the value is easy to change, low robustness
%     controllability_nondom_plus(i) = difference_hypervolume_plus/prod(abs(init_out)); %we normalize by the original value of the hypervolume (it's like we have finally computed a derivative of the hypervolume)
%     controllability_nondom_minus(i) = difference_hypervolume_minus/prod(abs(init_out));
%     controllability_nondom(i) = max(controllability_nondom_plus(i),controllability_nondom_minus(i));
%     %%
% end




% save('controllability_nondom_plus.mat',controllability_nondom_plus);
% save('controllability_nondom_minus.mat',controllability_nondom_minus);
% save('controllability_nondom.mat',controllability_nondom);

% %% PLOT controllability LOW BIOMASS
%
% load('non_dominated.mat');
% M=2;
% figure
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.04], 0.1, [0.1 0.01]);
% if ~make_it_tight,  clear subplot;  end
%
% %axis tight
%
% subplot(4,1,2:4);
% %subaxis(plot_rows,plot_columns,ix_plot, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
% x = non_dominated(ix_low_biomass(1:end-8),1); %here and in the following, we do not consider the last 6 points of low_biomass because they are actually medium biomass, and the last two because they have Infinite roustness index (because they are with 0 denominator in the computation)
% y = -non_dominated(ix_low_biomass(1:end-8),2);
% %z = non_dominated(:,3);
% %C = non_dominated(:,M+2);
% scatter(x,y,40,'red','*');
% hold on
% plot(x,y,'Color','black','LineWidth',0.3,'LineStyle','--'); %draws a black line passing through the Pareto points
%
% xlim([min(x) max(x)]);
% aux_xlim = xlim;
% xlabel('Biomass [h^{-1}]');
% ylabel('Oncoflux IDH1+IDH2 [mmolh^{-1}gDW^{-1}]');
%
%

%% l'altezza dallo stem e' la controlability stessa (prima serviva togliere 1 perche' in realta' calcolavo la robustness, non la controlability)
load('results_controllability.mat');
figure;
scatter(init_out(:,1),controllability_nondom);
figure;
ixs_not_inf = find(controllability_nondom < 1e20);
scatter(init_out(ixs_not_inf,1),log(controllability_nondom(ixs_not_inf)));

figure
stem_value = zeros(numel(controllability_nondom),1);
stem_value(ixs_not_inf) = log(controllability_nondom(ixs_not_inf));

cmap = parula(256);
color_index = ceil((stem_value-min(stem_value))/(max(stem_value)-min(stem_value))*(size(cmap,1)-1)) + 1;  %maps all stem_values into the interval [0+1; 255+1] = [1;256]

for i=1:numel(stem_value)
    if stem_value(i) ~= 0 
        stem(i,stem_value(i),'LineWidth',1,'Marker','square','Color',cmap(color_index(i),:));
    end
    hold on   %the stem are long if the controllability is high, and therefore the controllability_nondom_low_biomass(i) is high, because the late indicates actually fragility
end

% xlim(aux_xlim);
set(gca,'YTick',[]) %removes ticks
set(gca,'XTick',[])
x0=100;
y0=100;
width=550;
height=150;
set(gcf,'units','points','position',[x0,y0,width,height]);


% % %instructions for horizontal colorbar
% % ax=gca;
% % pos=get(gca,'pos');
% % set(gca,'pos',[pos(1) pos(2) pos(3) pos(4)*0.95]);
% % pos=get(gca,'pos');
% % hc=colorbar('location','northoutside','position',[pos(1) pos(2)+pos(4)+0.01 pos(3) 0.02]);
% % set(hc,'xaxisloc','top');
% % set(hc,'YTick',[])
% % text_leg = 'High controllability';
% % t = text(-1.98E-7,2.6E4,text_leg);
% % text_leg = 'Low controllability';
% % t = text(-0.385E-7,2.6E4,text_leg);


%instructions for colorbar on the left
ax=gca;
pos=get(gca,'pos');
set(gca,'pos',[pos(1) pos(2) pos(3) pos(4)*0.95]);
pos=get(gca,'pos');
hc=colorbar('location','westoutside','position',[pos(1)-0.03 pos(2) 0.02 pos(4)]);
set(hc,'YTick',[]) %removes ticks
set(hc,'XTick',[])
text_leg = 'HC';
t = text(-185,19,text_leg);
text_leg = 'LC';
t = text(-185,-37,text_leg);
text_leg = '0';
t = text(-160,0,text_leg);
box on


