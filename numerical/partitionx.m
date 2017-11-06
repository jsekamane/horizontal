function [L] = partitionx(d,k)
%PARTITIONX Partition with extreme competition. 
%   We split demand into as many lots as is needed for every supplier to be 
%   able to participate in all lots. This is an extreme case with
%   "full"/"perfect" competition. The total payments will be minimized, but
%   there will (most likely) be multiple (identical) equilibrium. We choose this 
%   approach, because we have yet to determine a method to optimally 
%   partition demand (i.e. a partition that minimizes total payments and
%   which also has an unqiue (pure-strategy) equilibrium.

    T = length(d);

    % High demand period
    %[~, tu] = max(d);
    sL = ceil( d / min(k) );
    % Number of lots, each with (at most) size; min(k)
    %sL(tu)
    % Lots (the last contains any residual demand, i.e. the size of it may be less than min(k))
    L = [ min(k)*ones(1,sL(1)-1) d(1)-(sL(1)-1)*min(k) zeros(1,max(sL)-sL(1)) ]'; % t = 1   
    if T>1
        L = [L'; min(k)*ones(1,sL(2)-1) d(2)-(sL(2)-1)*min(k) zeros(1,max(sL)-sL(2)) ]'; % t = 2
    end
    % TO-DO: double check how lots work when T=2, and demand is not integer (ie. who residuals are grouped across time).
    %        This is especially important for when f>0.


end

