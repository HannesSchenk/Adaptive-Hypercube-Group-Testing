function [P_sus] = HC_find_suspect_fast(Groups,Groups_pos,Groups_neg,C_label)
D = ndims(C_label);
switch numel(Groups_pos)
    case 0 % zero infected
        P_sus = [];
    case D % one infected
        P_sus = Groups(Groups_pos(1),:);
        for g = 2:numel(Groups_pos)
            P_sus = intersect(P_sus,Groups(Groups_pos(g),:));
        end
    otherwise % two or more infected
        if ~isempty(Groups_neg)
            P_save = rmmissing(Groups(Groups_neg(1),:));
            for g = 2:numel(Groups_neg) % des versog, wenn gruppen mit lei nan sein
                P_save = unique([P_save, Groups(Groups_neg(g),:)]);
            end
            N_tot_vec = C_label(:);
            N_tot_vec = (N_tot_vec(~isnan(N_tot_vec)));
            P_sus = N_tot_vec(~ismember(N_tot_vec,P_save));
        else
            P_sus = rmmissing(C_label(:));
        end
end
end