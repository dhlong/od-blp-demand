function [obj, deriv] = gmm(theta, Data)
%SSQ Summary of this function goes here
%   Detailed explanation goes here

[delta, s] = invertshare(theta, Data);

if any(isnan(delta))
    fprintf(' *** invert delta overflow\n');
    obj = 1e30;
    deriv = 1e30*ones(size(theta));
    return
end

X = Data.X;
Z = Data.Z;

XZ = X'*Z;
ZZ = Z'*Z;
XZZZ = XZ/ZZ;
beta = (XZZZ*XZ')\(XZZZ*Z'*delta);

xi = delta - X*beta;
xiZ = xi'*Z;
obj = xiZ/ZZ*xiZ';

if nargout > 1
	jab = jacob2(s, Data);
	dxi = jab - X*((XZZZ*XZ')\(XZZZ*Z'*jab));
	deriv = 2*xiZ/ZZ*Z'*dxi;
end

end

