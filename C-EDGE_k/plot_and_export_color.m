function plot_and_export_color
%      
% if (matlabpool('size') == 0)  %opens only if it is closed 
%    matlabpool('open','local',4)
% end

load('others.mat');
load('non_dominated.mat');

M = 2; %number of objective functions
     
% the following instructions are needed now before non_domination will always look
% for the minimisation of the objectives!
others(:,1)= - others(:,1);
others(:,2)= - others(:,2);
%others(:,4:end)=-others(:,4:end);
non_dominated(:,1) = - non_dominated(:,1);
non_dominated(:,2) = - non_dominated(:,2);
%non_dominated(:,4:end)=-non_dominated(:,4:end);


num_populations = max(abs(non_dominated(:,M+2)));
size_pop = max(abs(others(:,M+3)));



num_others = size(others,1);

others = sortrows(others,M+2);   %sorting the table according to the generation in which they have been found
%others=flipud(others);  %reverses the order. commenting this instruction makes the late generations hidden behind the early generations, and therefore they become visible only if that area of the Pareto front was not explored already by the early generations. 

x = others(:,1); 
y = others(:,2);




%colore(:,1) = 180*(1 - (others(:,4)-min(others(:,4)))/(max(others(:,4))-min(others(:,4))));  % others(i,4) contains the number (with the minus sign, though) of the population to which others(i,:) belongs. The +10 avoidt the points are too black and allows for more variability of colours
%colore(:,2) = 30+200*(1 - (others(:,4)-min(others(:,4)))/(max(others(:,4))-min(others(:,4))));
%colore(:,3) = 255;


%a = makeColorMap([180/255 230/255 255/255],[0/255 30/255 255/255],max(others(:,4))-min(others(:,4)));
%colormap(a);

%parfor_progress(num_others); % Initialize

% for i=1:num_others/100
%     %parfor_progress;
%     disp(num_others-i);
%     plot(x(i),y(i),'*','Color',colore(i,:)./255);
%     hold on
% end



%% PLOT OTHERS

C = others(:,M+2);
fig3=figure;
scatter(x,y,40,C,'*');

%set(gcf,'renderer','OpenGL')  %lets my laptop handle the graphics, rather than leaving it to MATLAB. This will be changed again later to export high quality pictures

colormap(jet(256)) 
%colorbar;
%colormap(flipud(colormap)) %reverses the colormap

% instructions for horizontal colorbars
ax=gca;
pos=get(gca,'pos');
set(gca,'pos',[pos(1) pos(2) pos(3) pos(4)*0.95]);
pos=get(gca,'pos');
hc=colorbar('location','northoutside','position',[pos(1) pos(2)+pos(4)+0.03 pos(3) 0.02]);
set(hc,'xaxisloc','top');
set(hc,'YTick', 10:10:max(C));




%% PLOT NON-DOMINATED

hold on
x = non_dominated(:,1);
y = non_dominated(:,2);
C = non_dominated(:,M+2);
plot(x,y,'*:', 'Color',[1 0 0]);
scatter(x,y,40,C,'*');

% K = convhull([x y]);  %convex hull of all the Pareto-optimal points, so we can plot the Pareto optimal surface trisurf
% trs = trisurf(K,x,y,C);
% set(trs, 'EdgeColor', 'none', 'FaceColor', 'black'); %Properties of the trisurf: facealpha sets transparency
% set(trs,'FaceAlpha',0.1);
grid on
hold on

set(fig3,'units','normalized','outerposition',[0 0 1 1]); %we need to force full screen, otherwise the xtick that we set a few instructions below will be only 4, and then if I make the picture full-screen through the maximize-window button of Windows, the picture will have now 8 ticks, but those 4 ticks will be repeated for the other 4 ticks, finally generating wrong ticks for the x axis
% aux = xlim;
% xlim([min(x) max(x)]);
xlab = get(gca,'xtick');

%%set label's decimal digits
% vers = version; %matlab version
% if ~isempty(strfind(vers, 'R2014b'))
%     set(gca,'xticklabel',cellstr(sprintf('%.12f\n',xlab))); %for Matlab r2014b
% else
%     set(gca,'XTickLabel',sprintf('%.12f|',xlab)); %for Matlab r2014a and previous versions
% end


 
title('');
xlabel('|EDGE_k - max{EDGE_{k-1}}| [h^{-1}]');
ylabel('k');




hold off

%% PLOT NON-DOMINATED ONLY
fig_nondom = figure;

x = non_dominated(:,1);
y = non_dominated(:,2);
C = non_dominated(:,M+2);
scatter(x,y,40,[0.5 0 0],'*');
grid on
hold on

set(fig_nondom,'units','normalized','outerposition',[0 0 1 1]); %we need to force full screen, otherwise the xtick that we set a few instructions below will be only 4, and then if I make the picture full-screen through the maximize-window button of Windows, the picture will have now 8 ticks, but those 4 ticks will be repeated for the other 4 ticks, finally generating wrong ticks for the x axis
aux = xlim;
xlim([min(x) max(x)]);
xlab = get(gca,'xtick');
vers = version; %matlab version
if ~isempty(strfind(vers, 'R2014b'))
    set(gca,'xticklabel',cellstr(sprintf('%.12f\n',xlab))); %for Matlab r2014b
else
    set(gca,'XTickLabel',sprintf('%.12f|',xlab)); %for Matlab r2014a and previous versions
end


title('Non-dominated points');
xlabel('|EDGE_k - max{EDGE_{k-1}}| [h^{-1}]');
ylabel('k');




%fill3(x',y',z',0);
%L = plot(x,y,'*:', 'Color',[0 0 100]./255);
%uistack(L, 'bottom');


% x = - non_dominated(:,1);
% y = non_dominated(:,2);
% fig2=figure(2);
% plot(x,y,'*:', 'Color','Black')
% grid on
% hold off
% title('');
% xlabel('Acetate [mmol h^{-1} gDW^{-1}]');
% ylabel('Biomass [h^{-1}]');
% 
% [h_m h_i]=inset(fig1,fig2);


set(gcf,'renderer','painters')  %to export high quality eps or pdf
% colormapeditor
% aux = xlim;
% xlim([0 aux(2)]);
% colorbar('YTick', 100:100:1500);

% exportfig(gcf, 'figura.pdf', 'color', 'cmyk', 'Width', '20', 'Height', '12', 'FontMode', 'scaled', 'FontSize', '1.2' );
