function jab = jacob2(s, Data)
%JACOB Summary of this function goes here
%   All random coefficients are normally distributed

iT = Data.iT;
N = size(Data.v,3);
T = max(iT);
% cdindex = [find(diff(iT)>0); length(iT)];

jab = zeros(size(Data.Xrc));

% indexing order: product - random coef - person (j,k,i)

s = permute(s, [1 3 2]);
sxv = bsxfun(@times, Data.XrcV, s);
% sumsxv = cumsum(sxv);
% sumsxv = sumsxv(cdindex,:,:);
% sumsxv = [sumsxv(1,:,:); diff(sumsxv)];
% ssxv = bsxfun(@times, sumsxv(iT,:,:), s);
% derShareTheta = mean(sxv - ssxv, 3);
% 
s = squeeze(s);

for t = 1:T
    index = iT == t;
    si = s(index,:);
    derShareDelta = diag(mean(si,2)) - si*si'/N;
    
    si = permute(si, [1 3 2]);
    sixv = sxv(index,:,:);
    sumsxv = sum(sixv);
    ssxv = bsxfun(@times, si, sumsxv);
    derShareTheta2 = mean(sixv - ssxv, 3);   
    
    jab(index,:)  = -derShareDelta\derShareTheta2;
end

end

