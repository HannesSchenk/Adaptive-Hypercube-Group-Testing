function [Groups,Groups_pos,Groups_neg] = GroupsPosNeg_fast(C_label,C_res,L)
M = [];
k = 1;
D = ndims(C_label);
Groups = nan(L*D,L^(D-1));
temp = true(L^D/L,1);
shift = L^D/L;
for i = 1:D
    for j = 1:L
        if i == 1
            b = [temp; repmat(~temp,L-1,1)];
            M = [M, circshift(b,(j-1)*shift)];
        else
            b = repmat([temp; repmat(~temp,L-1,1)],L^(i-1),1);
            M = [M, circshift(b,(j-1)*shift)];
        end
        A = C_label(logical(reshape(M(:,k),repmat(L,1,D))));
        A = (A(~isnan(A)));
        Groups(k,1:numel(A)) = A;
        k = k+1;
    end
    if size(temp,1) > 1
        temp = temp(1:end/L);
        shift = shift/L;
    end
end
 
Groups_tot = 1:L*D;
Groups_pos = [];
for g = 1:D*L % identify positive groups
    if any(C_res(logical(reshape(M(:,g),repmat(L,1,D)))))
        Groups_pos = [Groups_pos, g];
    end
end
Groups_neg = setxor(Groups_tot,Groups_pos);
% Groups_neg = setxor(Groups_neg,find(structfun(@isempty, Groups))); % deleting empty groups


end