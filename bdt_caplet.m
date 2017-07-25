function [f] = bdt_caplet(sig)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    global pai r_new d_new dt sqrt_dt U ;
    global K  caplet_market;
    global c;

    N = 1;
    k = 3*c+4;
    dt1 = 0.25;
    
    if c == 1
        a = 2;
    else
        a = k-3;
    end;

    for i = a:(k-1)
        for j = 1:i
            if j == 1
                pai(i,j) = 0.5 * pai(i-1,j) * d_new(i-1,j);
            elseif j == i
                pai(i,j) = 0.5 * pai(i-1,j-1) * d_new(i-1,j-1);
            else
                pai(i,j) = 0.5 * pai(i-1,j) * d_new(i-1,j) + 0.5 * pai(i-1,j-1) * d_new(i-1,j-1);
            end;
        end;
        
        for j = 1:i
            r_new(i,j) = U(i) * exp((2*(j-1)-(i-1))*sig*sqrt_dt);
        end;
        
        for j = 1:i
            d_new(i,j) = exp(-r_new(i,j)*dt);
        end;
    end;
    
    p = zeros(k,k);
    p(k,1:k) = ones(1,k);
    
    for i = (k-1):(-1):(k-3)
        for j = 1:i
            p(i,j) = (0.5 * p(i+1,j) + 0.5 * p(i+1,j+1)) * d_new(i,j);
        end;
    end;
     
    L = (1./p((k-3),1:(k-3)) - 1) * (1/dt1);
    payoff = N * dt1 * max(L-K(c),0).* p((k-3),1:(k-3));
    f = (sum(pai((k-3),1:(k-3)) .* payoff) - caplet_market(c));

end

