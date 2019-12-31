 for i=1:size(non_dominated)
for j = i:size(non_dominated)
a=non_dominated(i,1); b=non_dominated(i,2); c=non_dominated(i,3);
a1=non_dominated(j,1); b1=non_dominated(j,2); c1=non_dominated(j,3);
if ( a>a1 && b>b1 && c>c1)
disp('a')
end
end
end