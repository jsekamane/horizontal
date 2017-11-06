function [ED_q,p,rev,cost,profit] = eqhourlycomp(ED_q, d, k, f, c, cap)
%EQHOURLYCOMP We only focus on the `competitive` equilibrium, ie. we use that 
%   fact that we know the q_n,t in this equilibrium before hand (thus who 
%   wins and production costs), and so we only need to solve for (winning) 
%   bids.
%   EQHOURLYCOMP(ED_q, D, k, f, c, cap)
%   version 0.03
    
    % TO-DO: extend to cover f>0.

    % minimum bid increment  
    minbidincrement = 1;
    [sn, T] = size(ED_q);
    
    p = zeros(1,T);
    
    % Residual capacity
    rk = repmat(k',1,T)-ED_q;

    % Assume that inframarginal bid zero, and that
    %             non-dispatch bid startup cost + marginal cost.
    b = sign(rk) .* repmat((f+c)',1,T);

    % Find the marginal bidder
    % Sort by marginal costs
    [~, idx_c] = sort(c); %sortrows([c' k']); %sortrows([c' k'],[1 2],{'ascend' 'descend'});
    % cumulative capacity (sorted by marginal cost)
    kcs = cumsum( k(idx_c) );
    % Find the marginal bidder (ie. the ones who's capacity means that demand is met).
    [I,J] = find(d-repmat(kcs',1,T)<=0);
    [~,idx_mm] = unique(J, 'first'); % Find the first that makes the remaing demand zero or negative
    m = idx_c( I(idx_mm) );

    
    for t=1:T
        % Assume that marginal bidder bids 
        %   - the lowest bid of the non-dispatched (if the marginal bidder's cost is lower)
        %   - the lowest bid of the non-dispatched minus minimum bid increment (if cost are higher or equal)
        %   - the price cap if there are not any non-dispatched.
        %next = find(b(:,t) - b(m(t),t)>0,1); % the next non-dispatched
        idx_next = find(ED_q(:,t)==0);
        [~, idx_nexts] = sort(b(idx_next,t));
        next = idx_next(idx_nexts(1));
        if ~isempty(next)
            b(m(t),t) = b(next,t) - (b(next,t) > c(m(t)))  * (c(m(t)) >= c(next)) * minbidincrement;
        else
            b(m(t),t) = cap;
        end
        p(t) = b(m(t),t);
    end

    %b

    rev = (ED_q .* p)';
    cost = ( any(ED_q,2).*f' + ED_q.*repmat(c',1,T) )';
    profit = rev-cost;

end

