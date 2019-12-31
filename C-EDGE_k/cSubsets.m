%% The function below displays the indices of all subsets having c elements where length of the base array is N
function [indices, subsets] = cSubsets(N, c)
indices = dec2bin(1:(2 ^ N - 1)) - '0';   % indices of all subsets
counter = sum(indices, 2);                % number of elements in subsets

subsets = logical(indices(find(counter == c), :));  % logical indices appropriate subsets indices
end