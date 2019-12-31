function f  = genetic_operator(parent_chromosome, M, V, fbamodel, geni, reaction_expression, k_edge)

%% function f  = genetic_operator(parent_chromosome, M, V)
%
% This function is utilized to produce offsprings from parent chromosomes.
% The genetic operators corssover and mutation which are carried out with
% slight modifications from the original design. For more information read
% the document enclosed.
%
% parent_chromosome - the set of selected chromosomes.
% M - number of objective functions
% V - number of decision varaiables
% mu - distribution index for crossover (read the enlcosed pdf file)
% mum - distribution index for mutation (read the enclosed pdf file)
% l_limit - a vector of lower limit for the corresponding decsion variables
% u_limit - a vector of upper limit for the corresponding decsion variables
%
% The genetic operation is performed only on the decision variables, that
% is the first V elements in the chromosome vector.

%  Copyright (c) 2009, Aravind Seshadri
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are
%  met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%  POSSIBILITY OF SUCH DAMAGE.


[N,m] = size(parent_chromosome);
clear m
eps = 0.01;
value= k_edge ; %number of mutations allowed with respect to the parents
C = k_edge; % max number of knockouts (this equals k, i.e. the length of the subsets of genes for which the combined EDGE is evaluated)
NN = V+M;

%parfor
parfor i=1:N
    child=parent_chromosome(i,1:NN);
    j=1;
    while j<=value
        h_j=randi([1 V],1);
        if j==1
            child(j,h_j)=not(child(1,h_j));
        else
            child(j,h_j)=not(child(j-1,h_j));
        end
        
        while sum(child(j,1:V)) > C %reverts some KO so that their total number is C
            from = sum(child(j,1:V))-C;
            to = sum(child(j,1:V));
            n_knockin=randi([from to],1);
            geni_knockout=find(child(j,1:V)==1);
            
            knockin = randi([1 length(geni_knockout)],n_knockin,1);
            child(j,geni_knockout(knockin))=0;
        end
        
        while sum(child(j,1:V)) < C %adds some KO so that their total number is C
            geni_not_knockout=find(child(j,1:V)==0);
            additional_ko = randi([1 numel(geni_not_knockout)],1);
            child(j,geni_not_knockout(additional_ko)) = 1;
        end
        
        subset_selected = find(child(j,1:V) == 1); %the ones indicate the subset of k-genes selected for EDGE evaluation
        k = numel(subset_selected);
        child_eps = ones(1,V);
        child_eps(subset_selected) = eps;
        child_ko = ones(1,V);
        child_ko(subset_selected) = 0;
        out_k_eps = [-1 -1 1].*evaluate_objective_EDGE(child_eps, M, V, fbamodel, geni, reaction_expression);
        out_k_ko= [-1 -1 1].*evaluate_objective_EDGE(child_ko, M, V, fbamodel, geni, reaction_expression);
        EDGE_k = out_k_eps(1) - out_k_ko(1); % we only consider the biomass
        
        [indices, subsets] = cSubsets(k, k-1); %finds all the subsets of 1:k having k-1 elements, and we will compute their EDGE
        EDGE_kminus1 = zeros(size(subsets,1),1);
        for ixSubs = 1:size(subsets,1)
            child_eps = ones(1,V);
            child_eps(subset_selected(subsets(ixSubs,:))) = eps;
            child_ko = ones(1,V);
            child_ko(subset_selected(subsets(ixSubs,:))) = 0;
            out_kminus1_eps = [-1 -1 1].*evaluate_objective_EDGE(child_eps, M, V, fbamodel, geni, reaction_expression);
            out_kminus1_ko = [-1 -1 1].*evaluate_objective_EDGE(child_ko, M, V, fbamodel, geni, reaction_expression);
            EDGE_kminus1(ixSubs) = out_kminus1_eps(1) - out_kminus1_ko(1);
        end
        
        %         EDGE_kminus1 = max(EDGE_kminus1);
        %         if isempty(EDGE_kminus1) EDGE_kminus1=0; end     %prevents the case k=1 for which there are no subsets
        %
        %         EDGE_diff = abs(EDGE_k - EDGE_kminus1);
        
        EDGE_diff = 1;
        for ixSubs = 1: k
            EDGE_diff = EDGE_diff * abs(EDGE_k - EDGE_kminus1(ixSubs));
        end
        
        
        if EDGE_diff>0.00001
            disp('a');
        end
        
        child(j,V+1:NN) = [-EDGE_diff -EDGE_diff];  %we need to maximise the EDGE_diff (single objective optimization)
        fprintf('|EDGE_k - EDGE_(k-1)_1| * ... * |EDGE_k - EDGE_(k-1)_k|  = %.15f, k = %d\n',EDGE_diff, k);
         j=j+1;
    end% while
    child = non_domination_sort_mod(child, M, V);
    child_new(i,:) = replace_chromosome(child, M, V, 1);
end
f = child_new(:,1:NN);
