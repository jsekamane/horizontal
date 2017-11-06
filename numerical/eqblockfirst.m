function [q,lw,p,rev,cost,profit] = eqblockfirst(L, k, f, c, cap)
%EQBLOCKFIRST Find (the first) block auction equilibrium
%   i.e assume that the cheapest supplier wins the first x lots, where x
%   corresponds to this suppliers capacity. And so on for the second 
%   cheapest, third cheapest, etc.
%   EQBLOCKFIRST(L, k, f, c, cap)
%   version 0.01
    
    T = size(L,2);
    ssL = size(L,1);
    sn = length(k);
    
    % Prepare matix with quantity, bids and final price.
    q = zeros(sn,T);
    %ql = zeros(sn,ssL);
    lw = zeros(1,ssL);
    %b = zeros(sn,ssL);
    p = zeros(1,ssL);
    rev = zeros(1,sn);
    cost = zeros(1,sn);

    % Total size of each lot
    Ltot = sum(L,2);
    % Lot capacitiy requirement
    Lcap = max(L,[],2);
    
    %for t=1:size(L,2)
    % Remaning suppleirs, and thier capacity and costs
    rk = repmat(k,T,1)';
    rf = f;
    rc = c;

    % Loop through each lot
    for l = 1:ssL

        % Suppliers with sufficient remaning capacity:
        srk = all(L(l,:) <= rk, 2)';
        
        % cost of producing lot (for each supplier)
        cL = (rf + rc*Ltot(l)) .* srk + ...
             9999 * (1-srk) ;
        
        % sort the costs of
        [~, idx_cLs] = sort(cL);

        % winner of lot
        w = idx_cLs(1);
        lw(l) = w;
        q(w,:) = q(w,:) + L(l,:);
        
        % winner bids runner up's minimum bid (i.e. costs) and 
        % bids cap if there is no runner up.
        if( sum(srk)>1 )
            p(l) = cL(idx_cLs(2));
        else
            p(l) = cap;
        end
        
        % Revenue and costs
        rev(w) = rev(w) + p(l);
        cost(w) = cost(w) + cL(w);
        
        % reduce winners remaning capacity
        rk(w,:) = rk(w,:) - L(l,:);

        % Remove suppliers with zero remaing capacity
        %remove = find(k==0);
        %rk(remove) = [];
        %rf(remove) = [];
        %rc(remove) = [];

        %(Lcap(l) <= rk) .* cl
        %[~, scl] = sort(  )

    end
    
    profit = rev-cost;
    
    %lw
    %L
    %q
    %p
%    end

end

