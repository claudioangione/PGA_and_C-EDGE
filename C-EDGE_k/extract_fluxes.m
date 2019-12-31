function v1 = extract_fluxes(n_sol,n_gen,experiment_name)
%extracts the flux distribution for the n_sol population and the n_gen
%generation of that population

load('recon2_merged_bio_IDH1-2.mat');   %model with objectives and positions changed. Even if it's already in the workspace, must be reloaded every time this script starts because this script changes lower and upper bounds (multiplying something X old_bounds, so as to compute the output of FBA, so when a new output is required the orginal model must be reloaded because we need the original bounds

M=3;
%V=fbamodel.nbin;
V = numel(fbarecon.grRules);    %the length of the input individuals (without ranking, crowding distance and outputs) is equal to the number of reactions minus 2 that is the number of biomass and synthetic objective in the model (the final 2 reactions)

cd populations_recon2_merged_bio_IDH1-2
sol=['population' num2str(n_sol) '.mat'];
load(sol);
cd ..
x = chromosome(n_gen,1:V);



evaluate_objective(x, M, V, fbarecon, geni, reaction_expression)
