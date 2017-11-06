% SIMULATION 01
% Comparison of efficient dispatch, hourly-vertical auction and energy
% block auction. With:
% - Assumption #4 -- no pivotal generator: d <= sum(k)-max(m)
% - Assumption #3 -- marginal costs decreasing in capacity.
% - Assuming no fixed start up.

sim.size = 10000;
sim.seed = 20171103;
rng(sim.seed);

%% Generate random examples
for si = 1:sim.size
    [simulation(si).d, simulation(si).k, simulation(si).f, simulation(si).c] =  setuprandom('T',randi([1 2]),'fval',0);
end

save(sprintf('data/simulation_%i_%i.mat',sim.seed,sim.size), 'simulation')

%clear
%load('data/simulation_20171103_10000.mat')
%load('data/simulation_20171103_5.mat')

%% Calculate costs and payments
cap = 25+10*5*2 + 1;
results = table;
for si = 1:length(simulation)

    T = length(simulation(si).d);
    N = length(simulation(si).k);

    % Efficient dispatch
    [ED_q, ED_costs] = efficientdispatch(simulation(si).d, simulation(si).k, simulation(si).f, simulation(si).c);

    % Hourly-vertical auction
    [VA_q, p, rev, cost, ~] = eqhourlycomp(ED_q, simulation(si).d, simulation(si).k, simulation(si).f, simulation(si).c, cap);
    VA_costs = sum(cost(:));
    VA_payments = sum(rev(:));
    VA_cap = any(p == cap);

    % Partition (extreme competition)
    L = partitionx(simulation(si).d, simulation(si).k);

    % Find (the first) block auction equilibrium
    [HAx_q, ~, p, ~, cost, ~] = eqblockfirst(L, simulation(si).k, simulation(si).f, simulation(si).c, cap);
    HAx_costs = sum(cost(:));
    HAx_payments = sum(p);
    HAx_cap = any(p == cap);

    % Sanity check
    if ~isequal(ED_q, HAx_q) && ED_costs ~= HAx_costs
        warning('Simulation #%i -- Allocation is not equal and neither is costs.',si);
    end

    % summary variables
    dmin = min(simulation(si).d);
    dmax = max(simulation(si).d);

    ksum = sum(simulation(si).k);
    kmean = mean(simulation(si).k);
    kmedian = median(simulation(si).k);
    kmax = max(simulation(si).k);
    kstd = std(simulation(si).k);

    cmean = mean(simulation(si).c);
    cmedian = median(simulation(si).c);
    cmin = min(simulation(si).c);
    cmax = max(simulation(si).c);
    cstd = std(simulation(si).c);

    results = [results; ... 
               table(si, T, N, ...
                     ED_costs, VA_costs, VA_payments, VA_cap, HAx_costs, HAx_payments, HAx_cap, ...
                     dmin, dmax, ksum, kmean, kmedian, kmax, kstd, cmean, cmedian, cmin, cmax, cstd)];
end
%results.Properties.VariableNames = {'i' 'T' 'N' 'ED_costs' 'VA_costs' 'VA_payments' 'HAx_costs' 'HAx_payments'};

writetable(results, sprintf('data/simulation_%i_%i.csv',sim.seed,sim.size))
