function [delta, s] = invertshare(theta, Data)

% settings
toler = 1e-13;
display_iter = 1000;
max_iter_1 = 1000;
max_iter_2 = 10000;

logshare = log(Data.share);
global lastdelta

% starting value
if any((lastdelta > 1e30) | isnan(lastdelta))
    delta = logshare - log(Data.outshr);
else
    delta = lastdelta;
end

emu = exp(calmu(theta, Data));

%%
%
% Use Squared Polynomial Extrapolation Methods to speed up convergence,
% see the link for details
% <https://lirias.kuleuven.be/bitstream/123456789/482522/1/Discussion_Paper_35.pdf>
%
% Iterate: $\delta^{h+1} = \delta^h - 2\alpha^h + {\alpha^h}^2 v^h$ with
%
% $$g(\delta) = \delta + \ln S - \ln s(\delta)$$
%
% $$r^h = g(\delta^h) - \delta^h = \ln S - \ln s(\delta^h)$$
%
% $$v^h = g(g(\delta^h)) - 2g(\delta^h) + \delta^h = \ln S -
% \ln(\delta^h+r^h) - r^h$$
%
% $$\alpha^h = \frac{v^h \cdot r^h}{v^h \cdot v^h}$$
% or 
% $$\alpha^h = - sign(r^h \cdot v^h) \frac{||r^h||}{||v^h||}$$
%

% choose one of the stepsize functions
% stepsize = @(r,v) (v'*r)/(v'*v);
stepsize = @(r,v) -norm(r)/norm(v);

iter = 0;
converged = false;
while ~converged   
    s = calshare(delta, emu, Data.iT);
    r = logshare - log(mean(s,2));
    
    if iter > max_iter_1
        delta = delta + r;
    else    
        s2 = mean(calshare(delta + r, emu, Data.iT), 2);
        r2 = logshare - log(s2);

        v = r2 - r;
        a = stepsize(r,v);
        d = a^2*v - 2*a*r;
    
        delta = delta + d;
    end
    
    distance = max(abs(r(:)));
    converged = distance < toler;
    
    % display
    
    iter = iter + 1;
    if mod(iter, display_iter) == 0
        fprintf(' Iteration #%d, distance = %f\n', iter, distance);
    end
    
    % error handling
    
    if any((delta > 1e30) | isnan(delta))
        delta = nan(size(delta));
        fprintf(' *** Terminate invertshare; numerical overflow\n');
        return;
    end
    
    if iter == max_iter_1
        fprintf(' *** SPEM failed to converge; try regular contraction\n');
    elseif iter > max_iter_2
        delta = nan(size(delta));
        fprintf(' *** Terminate invertshare; max iterations reached; distance = %f\n', distance);
        return;
    end
end

lastdelta = delta;

end