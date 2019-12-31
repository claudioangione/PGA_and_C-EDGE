% FLUX_BALANCE Flux-balance analysis of FBA model
%    [V, FMAX, FMIN] = FLUX_BALANCE(FBAMODEL) performs a basic flux-balance
%    analysis of the FBA model FBAMODEL and returns a biomass-maximizing
%    flux distribution in the vector V.  The maximimum and minimum 
%    synthetic objective possible in a biomass-maximizing flux distribution
%    is given in FMAX and FMIN, respectively.
%
%    [V, FMAX, FMIN] = FLUX_BALANCE(FBAMODEL, QUIET) performs the analysis
%    and supresses screen output if QUIET is set to true.

function [v, fmax, fmax_max, fmax_min, fmin, fmin_max, fmin_min] = flux_balance_trilevel(fbamodel, quiet)

if nargin < 2
    quiet = false;
end

param.tmlim  = -1;
param.msglev = 1;
param.save   = 0;

nrxn   = numel(fbamodel.rxns);
nmetab = numel(fbamodel.mets);

yt = ones(nrxn,1); %the old yt = fbamodel.present, meaning that all the reactions must be considered as active in the model

A = [ fbamodel.S; 
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn) ];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1) ];
ctype = char('S' * ones(1, nmetab + nnz(~yt)));
vartype = char('C' * ones(1, nrxn));
[v, vbiomass] = glpk(fbamodel.f, A, b, fbamodel.lb, fbamodel.ub, ctype, vartype, -1, param);


%min and max 2nd objective fbamodel.g
A = [ fbamodel.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      fbamodel.f' ];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass ];
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 1));     %An array of characters containing the sense of each constraint in the constraint matrix.  Each element of the array may be one of the following values %           'F' Free (unbounded) variable (the constraint is ignored).  %           'U' Variable with upper bound ( A(i,:)*x <= b(i)).  %           'S' Fixed Variable (A(i,:)*x = b(i)).   %           'L' Variable with lower bound (A(i,:)*x >= b(i)).   %           'D' Double-bounded variable (A(i,:)*x >= -b(i) and A(i,:)*x <= b(i)).
[v, fmin] = glpk(fbamodel.g, A, b, fbamodel.lb, fbamodel.ub, ctype, vartype, 1, param);
[v, fmax] = glpk(fbamodel.g, A, b, fbamodel.lb, fbamodel.ub, ctype, vartype, -1, param);


%% CODICE AGGIUNTO DA ME PER FARE TRILEVEL

%min and max 3rd objective fbamodel.h with minimum 2nd objective fbamodel.g
A = [ fbamodel.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      fbamodel.f';
      fbamodel.g'];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass;
      fmin];
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 2));     %one 'S' more than the case before because there is another constraint
[v, fmin_min] = glpk(fbamodel.h, A, b, fbamodel.lb, fbamodel.ub, ctype, vartype, 1, param);
[v, fmin_max] = glpk(fbamodel.h, A, b, fbamodel.lb, fbamodel.ub, ctype, vartype, -1, param);


%min and max 3rd objective fbamodel.h with maximum 2nd objective fbamodel.g
A = [ fbamodel.S;
      sparse(1:nnz(~yt), find(~yt), ones(nnz(~yt), 1), nnz(~yt), nrxn);
      fbamodel.f';
      fbamodel.g'];
b = [ zeros(nmetab, 1);
      zeros(nnz(~yt), 1);
      vbiomass;
      fmax];
ctype = char('S' * ones(1, nmetab + nnz(~yt) + 2)); 
[v, fmax_min] = glpk(fbamodel.h, A, b, fbamodel.lb, fbamodel.ub, ctype, vartype, 1, param);
[v, fmax_max] = glpk(fbamodel.h, A, b, fbamodel.lb, fbamodel.ub, ctype, vartype, -1, param);



if ~quiet
    fprintf('Biomass flux          %s:    %f\n',fbamodel.rxns{find(fbamodel.f==1)}, fbamodel.f' * v);
    fprintf('1st Synthetic flux    %s:  [fmin = %.15f, fmax = %.15f]\n', fbamodel.rxns{find(fbamodel.g==1)}, fmin, fmax);
    fprintf('2nd Synthetic flux    %s:  [fmin_min = %.15f, fmin_max = %.15f, fmax_min = %.15f, fmax_max = %.15f]\n', fbamodel.rxns{find(fbamodel.h==1)}, fmin_min, fmin_max, fmax_min, fmax_max);
end