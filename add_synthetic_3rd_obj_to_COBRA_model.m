
% this script modifies the structure fbamodel coming out of the cobra
% routine readCbModel. The modification is that instead of a single natural
% objective array fbamodel.c, after running this script there are three arrays:
% fbamodel.f for the natural objective (biomass/growth), by default it's set
% to the last reaction
% fbamodel.g and fbamodel.h for the synthetic objective by default are set
% to the last but one and the last but two reactions

function fbamodel = add_synthetic_3rd_obj_to_COBRA_model(fbamodel)

f = fieldnames(fbamodel);
v = struct2cell(fbamodel);

posiz_c = strmatch('c',f,'exact');
num_fields = length(f);
f{posiz_c} = 'f';

f1 = cell(num_fields+2,1); %new fieldnames array with two locations more
f1(1:posiz_c) = f(1:posiz_c);
f1{posiz_c+1} = 'g'; %new fieldname 'g' in position 8
f1{posiz_c+2} = 'h'; %new fieldname 'g' in position 9
f1(posiz_c+3:num_fields+2)=f(posiz_c+1:num_fields);


v1 = cell(num_fields+2,1); %new field content array  with two locations more
v1(1:posiz_c) =v (1:posiz_c);
v1(posiz_c+1) = v1(posiz_c);  %temporary, the fbamodel.g array is equal to the fbamodel.f. This MUST be changed before running any FBA
v1(posiz_c+2) = v1(posiz_c);  %temporary, the fbamodel.h array is equal to the fbamodel.f. This MUST be changed before running any FBA
v1(posiz_c+3:num_fields+2) = v(posiz_c+1:num_fields);

fbamodel = cell2struct(v1,f1);


%% PLEASE SET OBJECTIVES USING CHANGE_OBJ
%fbamodel.g(end) = 0;
%fbamodel.g(end-1) = 1;

%set the natural objective to the last reaction (check biomass is always at
%the end if the biomass is needed as natural objective, that is the first one that gets maximisedminimised in the bilevel problem)

%fbamodel.f(find(fbamodel.f==1)) = 0;
%fbamodel.f(end) = 1;



