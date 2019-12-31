load('recon2_Nielsen.mat');
load('recon2.v03.mat');
fbarecon = reconNielsen;

vuoti_fbarecon = length(find(ismember(fbarecon.grRules,'')));
vuoti_recon2 = length(find(ismember(recon2.grRules,'')));

n_grRules = length(fbarecon.grRules);
for i = 1:n_grRules
    if strcmp('',fbarecon.grRules(i))
        ix_gr = find(strcmp(fbarecon.rxns{i},recon2.rxns));
        if ~strcmp('',recon2.grRules(ix_gr))  %skipped if also recon2.grRules(ix_gr) is empty
            fbarecon.grRules(i) = recon2.grRules(ix_gr);
        end
    end
end