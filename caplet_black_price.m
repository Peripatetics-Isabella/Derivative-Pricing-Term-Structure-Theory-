function [f] = caplet_black_price(sig)

    global forward_rate  K  tao  ZCB  cap  previous_caplet  caplet_market  caplet_black;
    global c;
    
    N_caplet = 1;
    if c ~= 1
            j = c;
            d1 = (log(forward_rate(j)/K(c)) + 0.5 * (sig^2) * (tao*(j))) / (sig * sqrt(tao*(j)));
            d2 = (log(forward_rate(j)/K(c)) - 0.5 * (sig^2) * (tao*(j))) / (sig * sqrt(tao*(j)));
        
            previous_caplet(j,c) = N_caplet * tao * ZCB(3*j+4) * (forward_rate(j) * normcdf(d1) - K(c) * normcdf(d2)); 
        
        
    end;
    
    previous_caplet_sum = sum(previous_caplet(:,c));
    caplet_market(c) = cap(c) - previous_caplet_sum;
    
    d1 = (log(forward_rate(c)/K(c)) + 0.5 * (sig^2) * (tao*(c))) / (sig * sqrt(tao*(c)));
    d2 = (log(forward_rate(c)/K(c)) - 0.5 * (sig^2) * (tao*(c))) / (sig * sqrt(tao*(c)));
        
    caplet_black(c) = N_caplet * tao * ZCB(3*c+4) * (forward_rate(c) * normcdf(d1) - K(c) * normcdf(d2));
    
    f = (caplet_market(c) - caplet_black(c))^2;
    

end

