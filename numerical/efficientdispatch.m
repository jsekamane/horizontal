function [xout,objfval] = efficientdispatch(d,k,f,c)
%EFFICIENTDISPATCH minimize overall generation costs
%   EFFICIENTDISPATCH(d,k,f,c)
%   version 0.02

T = size(d,2);

% The vector of variables looks as follows
% x = [u_1 u_2 ... u_N, q_1,1, q_1,2, q_2,1, q_2,2, ..., q_N,1, q_N,2]
% i.e. the first N elements are the binary varible u_n, while the 
% following 2*N elements are q_n,t (sorted by n and then t).

objf = [f repelem(c,T)]; % Objective function

% shorthand for length/size of various elements
sx = length(objf);
sn = length(f);
%sq = sx-sn;

intcon = 1:sx; % all variables are integers

% Capacity constraint; q_n,t less or equal u_n * k_n for all n,t
A = [-repelem(diag(k),T,1) eye(T*sn)];
b = zeros(T*sn,1);

% Demand constraint; Sum of q_n,t must equal d_t for all t.
Aeq = [zeros(T,sn) repmat(eye(T),1,sn)];
beq = d';

% Technical lower/upper bound constraints on variables u_n and q_n,t.
lb = zeros(sx,1); % Neither u_n nor q_n,t can be negative
ub = [ones(sn,1); repelem(k,T)']; % The ones() together with `lb` enforces that u_n is binary.

% Minimize objective function
options = optimoptions('intlinprog','Display','off');
[x,objfval,exitflag,output] = intlinprog(objf,intcon,A,b,Aeq,beq,lb,ub,options);

if(exitflag ~= 1)
    error("intlinprog did not converged to the solution x.");
    % Print intlinprog output:
    x
    objfval
    exitflag
    output    
end

xout = reshape(x(sn+1:end),T,sn)';

%[x(1:sn) reshape(x(sn+1:end),T,[])']
%testcost = x.*objf';
%sum([testcost(1:sn) reshape(testcost(sn+1:end),T,[])'], 2)

end

