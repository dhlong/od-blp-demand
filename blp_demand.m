clear;
df = importdata('..\..\data\od\od-annual-6.csv');
% df = importdata('.\blp-simulate.csv');
data = df.data;

%%
global lastdelta count outshr

cdid        = data(:,1);
firmid      = data(:,2);


const       = ones(size(cdid));
share       = data(:,4);
outshr      = data(:,5);
price       = data(:,6);

% dpm         = data(:,7);
% X1          = data(:,8);
% X2          = data(:,9);
% X3          = data(:,10);

mpd         = data(:,7);
dpm         = data(:,8);
mpg         = data(:,9);
hp          = data(:,10)/100;
weight      = data(:,11)/1000;
space       = data(:,12);
pgreal      = data(:,13);

suv         = data(:,14);
truck       = data(:,15);
van         = data(:,16);
minivan     = data(:,17);

comply      = data(:,18);
cafestd     = data(:,19)/1.4;


gpm         = 1./mpg*100; % gallons per 100 miles
dpm         = gpm.*pgreal;
hpwt        = hp./weight;
trend       = cdid-1;

count = 0;

rng('default');

% number of random draws per market
N = 500;

[T, ~, iT] = unique(cdid);
[F, ~, iF] = unique([cdid, firmid], 'rows');

nT = max(iT);

%% price rc and dpm rc will be dealt with separately

% random coefficients
Xrc_lb = 'price dpm pgreal const hpwt weight';
Xrc = [price dpm pgreal const hpwt weight]; % suv truck van minivan];
% Xrc = [price dpm const X1];
Krc = size(Xrc,2);

% mean utility coefficients
X_lb = 'price dpm pgreal const hpwt weight suv truck van minivan';
X = [(price) (dpm) pgreal const hpwt weight suv truck van minivan];
K = size(X,2);

Xz = [const pgreal (dpm) hpwt weight space suv truck van minivan];
for k = 1:size(Xz,2)
    sum_firm(:,k) = accumarray(iF, Xz(:,k));
    sum_total(:,k) = accumarray(iT, Xz(:,k));
end

count_firm = sum_firm(iF,1);
count_total = sum_total(iT,1);

Z_firm = bsxfun(@rdivide, sum_firm(iF,3:end) - Xz(:,3:end), max(count_firm - 1,1));
Z_rival = bsxfun(@rdivide, sum_total(iT,3:end) - sum_firm(iF,3:end), count_total - count_firm);
Z = [Xz count_firm Z_firm count_total Z_rival];

ZZ = Z'*Z; % GMM weighted matrix
XZ = X'*Z;
XZZZ = XZ/ZZ;
A = (XZZZ*XZ')\XZZZ*Z';


%% random draws
halton_dim  = Krc*nT; % last two are for price and dpm (log-normal)
halton_skip = 1000;
halton_leap = 100;
halton_scramble    = 'RR2';
tempDraw    = haltonset(halton_dim, 'Skip', halton_skip, 'Leap', halton_leap);
tempDraw    = scramble(tempDraw, halton_scramble);

% Make uniform draws
draws       = net(tempDraw, N);
v           = reshape(draws, [N nT (Krc)]);
v           = permute(v, [2 3 1]);
v           = norminv(v);

v = bsxfun(@minus, v, mean(v,3));
v = v(iT,:,:);
% load draws.mat;
% N = size(v,3);

%% combined data
Data.X = X;
Data.Z = Z;
Data.Xrc = Xrc;
Data.iT = iT;
Data.iF = iF;
Data.price = price;
Data.share = share;
Data.outshr = outshr;
% Data.v = v(iT,1:end-2,:);
Data.v = v;
Data.XrcV = bsxfun(@times, Data.Xrc, Data.v);
% Data.vprice = squeeze(v(iT,end-1,:)); % random draws for price
% Data.ve = squeeze(v(iT,end,:)); % random draws for dpm
Data.A = A;
Data.ZZ = ZZ;
% Data.gpm = gpm;
Data.dpm = dpm;
% Data.pgreal = pgreal;
% Data.comply = comply;
% Data.cafestd = cafestd;

%% Starting values and bounds

% IV logit delta
delta = log(share) - log(outshr);
beta = A*delta;
display(beta);


%%
theta0 =    [
    %sigma
    %     0.3362
    %     0.0258
    %     0.4837
    %     0.6828
    %     0
    0.5
    0.5
    0.5
    0.6
    0.5
    0.5
    %     0.5
    
    %
    %     1.5757 %sigma cartyve
    %     1.0700
    %     2.0370
    %    -0.4915
    
    % alpha lambda
    %    -0.1
    %    -0.1
    %     0.5
    %     1
    ];

thetalb = 0*ones(size(theta0));
thetaub = 100*ones(size(theta0));


lastdelta = delta;

%% Optimization options

solveropts = nloptset('algorithm', 'LN_BOBYQA');
optObj = optiset('display', 'iter', 'tolrfun', 1e-6, 'tolafun', 1e-6,...
    'maxtime', 1e5, 'solver', 'NLOPT', 'solverOpts', solveropts);

optObjIP = optiset('display', 'iter', 'tolrfun', 1e-6, 'tolafun', 1e-6,...
    'maxtime', 1e5, 'solver', 'NLOPT');

%% Optimization routine

OptIP = opti('fun', @(x) gmm(x, Data), 'x0', theta0, ...
    'lb', thetalb, 'ub', thetaub, 'options', optObj);

% [theta, fval] = solve(OptIP);
% [delta, s] = invertshare(theta, Data);

%% Matlab optim routine
options = optimset('Display', 'iter', 'TolFun', 1e-8, 'TolX', 1e-6,...
    'MaxIter', 2000, 'GradObj', 'on', 'DerivativeCheck', 'off');
[theta, fval] = fmincon(@(x) gmm(x,Data), theta0, [], [], [], [], thetalb, thetaub, [], options);

%%
Data2 = Data;
for iter = 1:5
    delta = invertshare(theta, Data);
    beta = A*delta;
    
    Z = optim_instrs(theta, beta, Data);
    ZZ = Z'*Z; % GMM weighted matrix
    XZ = X'*Z;
    XZZZ = XZ/ZZ;
    A = (XZZZ*XZ')\XZZZ*Z';
    Data.A = A;
    Data.ZZ = ZZ;
    Data.Z = Z;
    
    [theta, fval] = fmincon(@(x) gmm(x,Data), theta0, [], [], [], [], thetalb, thetaub, [], options);
    OptIP = opti('fun', @(x) gmm(x, Data), 'x0', theta, ...
        'lb', thetalb, 'ub', thetaub, 'options', optObjIP);
    
    [theta, fval] = solve(OptIP);
    Data = Data2;
    
end

%%
theta2 = theta;
[delta,s] = invertshare(theta2, Data);
beta2 = A*delta;

alphai = beta2(1) + theta2(1)*squeeze(v(:,1,:));
% alphai = bsxfun(@rdivide, alphai, price);
e = elas(s, alphai, Data.iT);
owne = diag(e).*price./share;

margin = calmargin(s, alphai, Data.iF);
markup = margin./price;
V = cov_demand(theta2, beta2, Data);
se = sqrt(diag(V));
% se = sqrt((diag(V)));
% display([theta, se(1:numel(theta))]);
printmat([theta2, se(1:numel(theta))], 'sigma', Xrc_lb, 'sigma se');
printmat([beta2, se(numel(theta)+1:end)], 'beta', X_lb, 'beta se');