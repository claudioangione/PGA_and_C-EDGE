function f = evaluate_objective(x, M, V, fbarecon, geni, reaction_expression)

% gene_importance decide l'intensita' che il numero di gene expression (normalizzato
% a 1) ha sui lower e upper bound dei flussi.


%% function f = evaluate_objective(x, M, V)
% Function to evaluate the objective functions for the given input vector
% x. x is an array of decision variables and f(1), f(2), etc are the
% objective functions. The algorithm always minimizes the objective
% function hence if you would like to maximize the function then multiply
% the function by negative one. M is the numebr of objective functions and
% V is the number of decision variables. 
%
% This functions is basically written by the user who defines his/her own
% objective function. Make sure that the M and V matches your initial user
% input.

yt=x';      % x' is the transpose of x, that is the gene expression array associated with the child (two codon usages)
yt=yt(1:V);     % needed because sometimes, and especially with child_3 in genetic_operator, there is a bug and all the child is passed to this function here, including the final rank, crowding distance and objective functions (instead of passing only the V decision variables)
%fbarecon.present=ones(fbarecon.nrxn,1);

%fbamodel.present = ~(fbamodel.G' * yt);   %istruzione vecchio approccio knockout

%NEW CODE: DATA UN GENE EXPRESSION ARRAY, CALCOLA AUTOMATICAMENTE GLI OUTPUT DELL'FBA

%load('geni.mat');
%load('reaction_expression.mat');

%geni( cellfun(@isempty, reaction_expression) ) = {'1'};  %replaces all the empty names of genes with the name '1'
%reaction_expression( cellfun(@isempty, reaction_expression) ) = {'1'};  %replaces all the empty cells of gene expression (e.g. exchange reactions) with '1', i.e. gene expressed nomally

%reaction_expression = cellfun(@(c) c{1},reaction_expression);
%geni = cellfun(@(c) c{1},geni);

eval_reaction_expression = cell(numel(reaction_expression),1);
ixs = find(yt~=1);

for ix=1:numel(ixs) %loop over the array of the non-1 gene expressions, in order to replace the names of genes in geni_reazioni.mat with their values. All the gene set expressions with only 1s as gene values , will be left empty and at the end of this loop everything empty will be substituted with 1 anyway. This avoids looping over all the genes yt, which is very expensive
     i =ixs(ix); 
     matches = strfind(reaction_expression,['/' geni{i} '/']);     %this and the following instruction find the locations of the gene 'bXXXX' in the array reaction_expression
     posizioni_gene = find(~cellfun('isempty', matches));
     for j=1:length(posizioni_gene)      %for each string of reaction_expression, we replace the substring 'bXXXX' with the number representing its gene expression
         if isempty(eval_reaction_expression{posizioni_gene(j)})
             eval_reaction_expression{posizioni_gene(j)} = reaction_expression{posizioni_gene(j)};
         end
         eval_reaction_expression{posizioni_gene(j)} = strrep(eval_reaction_expression{posizioni_gene(j)}, ['/' geni{i} '/'], num2str(yt(i),'%.15f'));  %Matlab strangely truncates decimal digits when using num2str. Addimg %.12f at least ensures that 12 decimal digits are included in the number converted into string
     end
end

for i=1:numel(eval_reaction_expression)
    if isempty(eval_reaction_expression{i})
        eval_reaction_expression{i} = '1';
    else
        for j=1:numel(geni) %in the eval_reaction_expression in which already all the non-1 genes have been substituted with their values, we now have to substitute all the remaining genes with the number 1
            eval_reaction_expression{i} = strrep(eval_reaction_expression{i}, ['/' geni{j} '/'], '1.0');  %we do this job only for the non-1 reaction expressions, so we save a lot of computational time. Note that we need '1.0' because the regecpr below is looking for NUM.NUM strings, and therefore putting only '1' will give an infinite loop
        end
    end
end


%eval_reaction_expression( cellfun(@isempty, eval_reaction_expression) ) = {'1'};  %replaces all the empty cells of gene expression (e.g. exchange reactions or reactions whose genes have all gene expression 1) with 1, i.e. gene expressed nomally

num_reaction_expression = zeros(1,length(eval_reaction_expression));
gamma = zeros(1,length(reaction_expression));

for i=1:length(num_reaction_expression)
    str = eval_reaction_expression{i};
    
    while (numel(strfind(str,')')) > 32) %if there are more than 32 parentheses, matlab is unable to run EVAL. So we need to reduce these parentheses manually by starting to eval smaller pieces of the string 
        to_replace = 'min.\d*+\.+\d*,\d*+\.+\d*.|max.\d*+\.+\d*,\d*+\.+\d*.';  %searches for all the strings of kind min(NUM.NUM,NUM.NUM) or max(NUM.NUM,NUM.NUM)
        substrings_to_replace = regexp(str, to_replace, 'match');
        for j = 1:numel(substrings_to_replace)
            ss_rep = substrings_to_replace{j};
            str = strrep(str,ss_rep,num2str(eval(ss_rep),'%.15f'));
        end              
    end
    
    num_reaction_expression(i) = eval(str);   %evaluates the cells like they are numerical expressions (so as to compute min and max of gene expressions)

%     match_geneset = find(yt == num_reaction_expression(i));     %this instruction finds the location of the resultant geneset expression in the gene expression array, so as to understand the gene responsible for this geneset expression
%     
%     if length(match_geneset) > 1
%         gamma(i) = 1; %in this case, there are multiple matches, which meansthe expression was exactly 1, which means almost surely that the gene in the model was not present in the probe genes, and therefore its expression and its variance must be set to 1
%     else
%         gamma(i) = gene_importance(match_geneset);
%     end
end

gamma = 3.*ones(1,length(reaction_expression));


for i=1:length(num_reaction_expression)   %loop over the array of the geneset expressions
    
        fbarecon.lb(i) = fbarecon.lb(i)*num_reaction_expression(i)^gamma(i);
        fbarecon.ub(i) = fbarecon.ub(i)*num_reaction_expression(i)^gamma(i);
        
%         if num_reaction_expression(i)>=1
%             fbarecon.lb(i) = fbarecon.lb(i)*(1+gamma(i)*log(num_reaction_expression(i)));
%             fbarecon.ub(i) = fbarecon.ub(i)*(1+gamma(i)*log(num_reaction_expression(i)));
%         else
%             fbarecon.lb(i) = fbarecon.lb(i)/(1+gamma(i)*abs(log(num_reaction_expression(i))));
%             fbarecon.ub(i) = fbarecon.ub(i)/(1+gamma(i)*abs(log(num_reaction_expression(i))));
%         end
end
%END NEW CODE

%[v1, fmax, fmax_max] = flux_balance_trilevel(fbarecon,true);
[v1, fmax, fmin] = flux_balance(fbarecon,true);

% objective functions number M is 2
f(1) = fbarecon.f' * v1; % Biomass
f(2) = fmax; % max of the 1st synthetic obj (we want to minimise the max flux in the IDH reactions)
f(3) = sum(abs(round(yt*1e15)/1e15-ones(numel(yt),1))); % !!!! change also in expFBA at the beginning for the initial population !!!!! (the 1e15 rounding allows to truncate at the 15th decimal digit so as to avoid nonzero very small values due to the loss of precision of matlab when subtracting two arrays)
format longG; format compact;
%f 
f(1) = -f(1); %put - if we maximise, + if we minimise (because nsga is by default minimising the objectives, so the mminus sign means it will minimise the negative objective, which means it will maximise the objective
f(2) = -f(2);
f(3) = f(3); % !!!! change also in expFBA at the beginning for the initial population !!!!! (NOT NEEDE IN THE NEW VERSION WHERE I EVALUATE ALSO THE INITIAL INDIVIDUALS USING EVALUATE_OBJECTIVE
f

%[concentrationMatrix,excRxnNames,timeVec,biomassVec] = my_dynamicFBA_antibiotic(fbamodel,substrateRxns,initConcentrations,initBiomass,timeStep,nSteps,plotRxns);