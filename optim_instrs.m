function [ instrs ] = optim_instrs(theta, beta, Data)
%OPTIM_INSTRS Summary of this function goes here
%   Detailed explanation goes here

pindex = 1;
prcindex = 1;

price = Data.X(:,pindex);
Z = Data.Z;

% 1st stage fitted price
phat = Z*((Z'*Z)\(Z'*price));

% if price has rc, replace with phat

Data.X(:,pindex) = phat;

if prcindex > 0
    Data.Xrc(:,prcindex) = phat;
    Data.XrcV = bsxfun(@times, Data.Xrc, Data.v);
end

% evaluate delta at xi = 0
delta = Data.X*beta;
mu = calmu(theta,Data);
emu = exp(mu);

s = calshare(delta, emu, Data.iT);
jab = jacob2(s, Data);

instrs = [Data.X jab];


end

