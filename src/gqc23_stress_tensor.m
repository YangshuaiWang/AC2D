% function STh = gqc23_stress_tensor(STh, S, ufs, divFC)
% Program: compute the corrected stress for gqc23 coupling method by 
%         applying divergence free tensor field.
% Author : M.Liao
% Version: Fist release     Mar-08-2016
function STh = gqc23_stress_tensor(STh, S, ufs, divFC)

for i = 1:size(S,2)
%     s = S(2:5,i);
    iT = S(1,i);
    it = find(ufs(1,:) == iT);
    ie = ufs(2,it);
    for j = 1:numel(ie)
        ic(j) = find(divFC(1,:)==ie(j));
    end
    dfc = zeros(2,2);
    for j = 1:numel(ie)
        dfc = dfc + divFC(2:3,ic(j))*ufs(3:4, it(j))';
    end
    STh(:,:,iT) = STh(:,:,iT) + dfc;
end

end