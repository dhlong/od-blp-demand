function [ e ] = elas( s, alphai, iT )
%ELAS Summary of this function goes here
%   Detailed explanation goes here

T = max(iT);
J = length(iT);
N = size(s,2);

sa = bsxfun(@times, s, alphai);
e = zeros(J);

for t = 1:T;
    index = iT == t;
    e(index, index) = diag(mean(sa(index,:),2)) - sa(index,:)*s(index,:)'/N;
end

end

