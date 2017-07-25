function [f] = lmm_initial_test_calibration(x)

    global sigma  rho  k;
    
    dt = 0.25;
    
    beta = x(1,1);
    a = x(1,2);
    bit = x(1,3);
    c = x(1,4);
    d = x(1,5);
    
    fai = x(:,6);
    
    % caplet_volatility
    sigma = zeros(1,k);
    for T = 1:k
        sigma(T) = fai(T) * ((a+bit) * (T*dt) * exp(-c * (T*dt)) + d);
    end;
    
    % correlation
    for i=1:k
        for j=1:k
            rho(i,j) = exp(-abs(i-j)*dt*beta);
        end;
    end;
    
    f = 0;
    
 end

