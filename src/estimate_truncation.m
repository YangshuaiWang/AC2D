% function Tr = estimate_truncation(geom, STh, dW)
% program: computing truncation error upper bound as sqrt(\s(iT_Bdry)^2)
% author : M.Liao
% vertion: first release    Apr-28-2016
function Tr = estimate_truncation(geom, STh, model)
[~, dW] = model.Wfun(model, model.F0);
dW = 2*dW/sqrt(3);
dW = reshape(dW, 4, 1);
N = geom.N;
iX = geom.iBdry;
iT = geom_2dtri_aIXT(geom, iX);
S = STh(:,:,iT);
S = reshape(S, 4, numel(iT));
% h = geom.maxH(iT);
volT = geom.volT*sqrt(3)/2;
vol = volT(iT);
vol = vol/sum(vol)*1/3*sum(volT);
% disp(num2str(sum(vol)));
ds = bsxfun(@minus, S, dW);
% ds = bsxfun(@times, ds, vol);
% vol = vol/sum(vol)*.75*sum(volT);
Tr = sqrt(sum(vol.*sum(ds.^2)));
% Tr = sqrt(sum(S.^2));
% Tr = 3*sqrt(3)/2*(2*N-1)/sum(vol)*sqrt(sum(vol.*sum(ds.^2)));
% Tr = sqrt(3*sqrt(3)/2*(2*N-1)/sum(vol)*sum(vol.*sum(ds.^2)));
% Tr = sqrt(sum(sum(ds.^2)));
end