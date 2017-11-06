function [d, k, f, c] = setuprandom(varargin)
%SETUPRANDOM Random demand, as well as capacity and costs of generators.
%   SETUPRANDOM('T',1,'N',randi([3 10]),'fval',NaN)
%   version 0.02

ip = inputParser;
addOptional(ip, 'T', 1, @isnumeric);
addOptional(ip, 'N', randi([3 10]), @isnumeric);
addOptional(ip,'fval', NaN, @isnumeric);
parse(ip,varargin{:})
T = ip.Results.T;
N = ip.Results.N;
fval = ip.Results.fval;
%T = 1; N = randi([3 10]); fval = NaN;


k = randi([1 5], 1, N);
k = floor(k/min(k)); % Rescale so min(k_n) always equal 1.

% Assumption #4: 
% No one generator is "pivitoal" supplier / necessary capacitiy to meet 
% demand without any one supplier.
d = randi([2 sum(k)-max(k)], 1, T);

% Assumption #3:
% c is strictly decreasing in f and k is increasing in f.
[~, idx_k] = sort(k);
if ~isnan(fval)
    f = fval * ones(1,N);
else
    f = randi([0 25],1,N); % fixed startup cost is random integer between 0 and 250.
    f(idx_k) = sort(f); % k is increasing in f.
end

%c = zeros(1,N);
strictly = false;
while ~strictly
    c = randi([1 10],1,N); % marginal cost is random integer between 1 and 100.
    c(idx_k) = sort(c,'descend'); % c is strictly decreasing in f.

    if isnan(fval)
        strictly = all( sign(diff(c)) == -sign(diff(f)) ); % strictly decreasing
    else
        strictly = true;
    end
end

[tempk, finalorder] = sortrows([k' c'],[1 2],{'descend' 'ascend'}); %sort(k,'descend');
k = tempk(:,1)';
f = f(finalorder);
c = c(finalorder);

end

