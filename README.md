# C-EDGE and PGA
Code for the paper: "Discovering code Essential Multiple Gene Effects through Large Scale Optimization: an Application to Human Cancer Metabolism"

--------------------
  Initial settings
--------------------

Download COBRA toolbox for MATLAB from http://opencobra.github.io/ and set the local COBRA folder with the instruction
```matlab
addpath(genpath('path_to_COBRA_toolbox'));
```

The following code is fully parallelised and requires the Parallel Toolbox in Matlab.


---------------------------------------------------
  Multi-objective optimization of gene expression
---------------------------------------------------

Load the model (e.g. the one with biomass and phosphoglycerate dehydrogenase set as objectives)
>> load('recon2_merged_bio_PHGDH.mat')

To start the optimization, run
```matlab
NUMBER_OF_CORES = 4  % please change depending on the number of cores available
RUN(128,384,NUMBER_OF_CORES)
```

Default population of 128 individuals; 384 populations will be generated by default. We suggest keeping this 1:3 proportion.
This will run the optimization with biomass as first objective and PHGDH as second objective. 
To change the objective reactions modify the following variables:

a) fbarecon.f selects the first objective (default: biomass)
b) fbarecon.g selects the second objective (default: phosphoglycerate dehydrogenase)

After the optimization, append_and_plot_solutions.m computes the Pareto front.
The file non_dominated.mat contains all the Pareto optimal points, while others.mat are all the other points. 
In the first two columns there are the two objective functions, while the 4th column is the number of the generation (that is, the file solutionX) in which that solution has been found, and the 5th column is the position of that solution in that generation.

Finally, plot_and_export_color.m plots the final version of the Pareto front, and statistics_on_genes.m generates the clustering and multidimensional scaling plots.


To find the indices of the reactions, and then change lower and upper bounds (fbamodel.lb and fbamodel.ub), please type
```matlab
n = find(ismember(fbarecon.rxns, 'EX_succ(e)')==1) (e.g. for succinate)
fbarecon.lb(n) = NEW_LOWER_BOUND
fbarecon.ub(n) = NEW_UPPER_BOUND
```


---------------------------------------------------
  C-EDGE algorithm
---------------------------------------------------

To run C-EDGE_1 and C-EDGE_2, run the compute_EDGE.m Matlab script.
This will load the recon2_merged_bio_PHGDH.mat metabolic model and run the C-EDGE algorithm.
The script will save the resulting vector as "c-edge_scores.mat" reporting the C-EDGE scores for each gene in the model.

To run C-EDGE_k in general, run the RUN_CEDGEk.m, matlab script in the subfolder C-EDGE_k

```matlab
RUN_CEDGEk(4) % please change depending on the number of cores available
```

In this case, the single objctive PGA (soPGA) will run. We are interested in finding the largest subsets possible of k genes such that if we compute the edge for those subsets we get interesting results, i.e. the EDGE(k) is different from all the EDGE(k-1) of the subsets of the k genes with k-1 elements. That is why, in genetic_operator.m, we maximise the objective
EDGE_diff = abs(EDGE_k - EDGE_kminus1);
