function f = non_domination(x, M, V)

[N, m] = size(x);

% Number of individuals that dominate this individual
individual = zeros(N,1);
tic
for i = 1 : N
   i
     for j = 1 : N
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1 : M   %loops over the M objecties
            if x(i,V + k) < x(j,V + k)
                dom_less = dom_less + 1;
            elseif (x(i,V + k) == x(j,V + k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more + 1;
            end
        end
        if dom_less == 0 && dom_equal ~= M  %which means the i-th individual is dominated by the j-th individual in at least one objective function
            individual(i) = individual(i)+ 1;
            break;
        end
    end   
    if individual(i) == 0
        x(i,M + V + 1) = 1; %if no individuals dominate the i-th individual, than put 1 in his M+V+1-th column
        
    end
   
end
toc
f=x;
end