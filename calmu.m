function mu = calmu(theta, Data )
%CALMU Summary of this function goes here
%   Detailed explanation goes here

%%
%
% $$\mu_{jti} = \delta_{jt} + \sum_k \sigma_k X_{jk} \nu_{jki} +
% \alpha \exp(\sigma^p \nu^p_{ki} price_{jt}) +
% \lambda \exp(\sigma^e \nu^e_{ki} dpm_{jt})$$
%

% N = size(Data.v, 3);
% J = size(Data.Xrc, 1);

params = getParams(theta);
sigmaXrc = bsxfun(@times, params.sigma', Data.Xrc);
mu = squeeze(sum(bsxfun(@times, sigmaXrc, Data.v),2));

% for i=1:N
%     alphai = params.alpha*exp(params.sigmap*Data.vprice(:,i));
%     lambdai = params.lambda*exp(params.sigmae*Data.ve(:,i));
%     mu(:,i) = sum(sigmaXrc.*Data.v(:,:,i),2) + ...
%         + alphai.*Data.price + lambdai.*Data.dpm;    
% end


end

