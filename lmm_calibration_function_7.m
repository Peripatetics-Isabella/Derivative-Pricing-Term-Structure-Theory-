function [f] = lmm_calibration_function_7(x)

    global black_vol  w  sigma  rho  swaption_sigma  swaption_volatility  k  forward_rate  swap_rate_lmm;
    
    dt = 0.25;
    alpha = x(1,1);
    beta = x(1,2);  % not fix beta
    % beta = 0.1;     % fix beta = 0.1
    gamma = x(1,3);
    a = x(1,4);
    b7 = x(1,5);
    c = x(1,6);
    d = x(1,7);
    
    fai_square = x(:,8);
    
    % caplet_volatility
    sigma = zeros(1,k);
    for T = 1:k
        integral = 0;
        for t = 1:T
            integral = integral + ((b7*(T*dt-t*dt)+a)*exp(-c*(T*dt-t*dt))+d)^2 * dt;
        end;
        sigma(T) = fai_square(T) / (T*dt) * integral;
    end;
        
    A = sum((sigma-black_vol).^2);
    
    % correlation
    for i=1:k
        for j=1:k
            rho(i,j) = alpha + (1-alpha)*exp(-abs(i-j)*dt*beta*exp(gamma*dt*min(i,j)));
        end;
    end;
    
    % swaption_volatility
    swaption_sigma = zeros(8);
    for tn = 2:8
        for t0 = 1:(tn-1)
            for i = 1:k*dt
                for j = 1:k*dt
                    I = 0;
                    for l = 0.25:dt:t0
                        I = I + rho(i/dt,j/dt) * sigma(i/dt) * sigma(j/dt) * dt;
                    end;
                    swaption_sigma(t0,tn) = swaption_sigma(t0,tn) + (w(i,tn)*w(j,tn)*forward_rate(4*i)*forward_rate(4*j)/(swap_rate_lmm(t0,tn)^2))*I;
                end;
            end;
            swaption_sigma(t0,tn) = sqrt( 1/t0 * swaption_sigma(t0,tn) );
            
        end;
    end;
    
    B = sum(sum((swaption_sigma-swaption_volatility).^2));
        
    f = A + B;

end

