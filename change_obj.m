function fbarecon = change_obj(fbarecon, num_obj, name_reaction)

index_new = find(strcmp(fbarecon.rxns,name_reaction)==1);
if isempty(index_new)
    tmp=sprintf('Reaction %s is not in the model', name_reaction);
    error (tmp);
    return;

end

switch num_obj
    case 1
        index_old = find(fbarecon.f==1);
        fbarecon.f(index_old) = 0;
        fbarecon.f(index_new) = 1;
    case 2
        index_old = find(fbarecon.g==1);
        fbarecon.g(index_old) = 0;
        fbarecon.g(index_new) = 1;
    case 3
        index_old = find(fbarecon.h==1);
        fbarecon.h(index_old) = 0;
        fbarecon.h(index_new) = 1;
end

