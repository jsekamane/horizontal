%% Setup
%seed = 1;
%rng(seed);
T = 2;
N = 3;

[d, k, f, c] = setuprandom()
[d, k, f, c] = setuprandom('T',2,'N',3,'fval',0)
cap = 25+10*5*2;


%% Efficient dispatch
[ED_q, ED_costs] = efficientdispatch(d, k, f, c)

%% Hourly-vertical auction
% We only focus on the `competitive` equilibrium, ie. we use that fact that
% we know the q_n,t in this equilibrium before hand (thus who wins and 
% production costs), and so we only need to solve for (winning) bids.

[VA_q,p,rev,cost,profit] = eqhourlycomp(ED_q, d, k, f, c, cap);

%VA_q;
VA_costs = sum(cost(:))
%VA_profit = sum(profit(:))
VA_payments = sum(rev(:))



%% Partition (extreme competition)
% We split demand into as many lots as is needed for every supplier to be 
% able to participate in all lots. This is an extreme case with
% "full"/"perfect" competition. The total payments will be minimized, but
% there will (most likely) be multiple (identical) equilibrium. We choose this 
% approach, because we have yet to determine a method to optimally 
% partition demand (i.e. a partition that minimizes total payments and
% which also has an unqiue (pure-strategy) equilibrium.

L = partitionx(d,k);


%% Find (the first) block auction equilibrium
% i.e. assume that the cheapest supplier wins the first x lots, where x
% corresponds to this suppliers capacity. And so on for the second cheapest,
% third cheapest, etc.
[HAx_q, lw, p, rev, cost, profit] = eqblockfirst(L, k, f, c, cap);

%HAx_q
HAx_costs = sum(cost(:))
%HAx_profit = sum(profit(:))
HAx_payments = sum(p) % = sum(rev(:))

if ~isequal(ED_q, HAx_q)
    warning("Allocation is not equal!")
    if(ED_costs ~= HAx_costs)
        warning("... and neither is costs!!")
    end
end
