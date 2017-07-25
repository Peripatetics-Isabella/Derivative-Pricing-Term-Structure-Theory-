%generate drift term
k = 32;
fr0 = forward_rate_list(2:2+k-1)/100;
phi = diag(normrnd(0,1,[1,k]));
for i=1:k
    store = 0;
    for j=(i+1):k
        store = store+fr0(j)*tao*rho(i,j)*sigma(j)/(1+fr0(j)*tao);
    end;
    mu(i) = -sigma(i)*store;
end;

par=(sum(b*sqrt(dt)*phi,2))';

fr1_hat = transpose(fr0).*exp(mu*tao-0.5*sigma.^2*tao+sigma.*par);

for i=1:k
    s = 0;
    for j=(i+1):k
        s = s+fr1_hat(j)*tao*rho(i,j)*sigma(j)/(1+fr1_hat(j)*tao);
    end
    mu_hat(i) = -sigma(i)*s;
end;

mu_new = 0.5*(mu+mu_hat);

% use all the above outputs to generate all F(T0);
fr1 = transpose(fr0).*exp(mu_new*tao-0.5*sigma.^2*tao+sigma.*par);

% construct the long jump forward rates matrix
long_jump_fr = [];

long_jump_fr = [long_jump_fr;fr1'];

for T=2:k
    
    % grasp the sigma(t)
    sigma_loop = sigma(T:k);    
    
    % generate new correlation matrix
    alpha_loop = 0.05;
    beta_loop = 0.1;
    gamma_loop = 0.05;
    for i=1:(k-T+1)
        for j=1:(k-T+1)
            rho_loop(i,j) = alpha_loop + (1-alpha_loop)*exp(-abs(i-j)*tao*beta_loop*exp(gamma_loop*tao*min(i,j)));
        end;
    end;
    
    % generate square root of the correlation matrix
    b_loop = sqrtm(rho_loop);
    
    % generate new drift term
    fr0_loop = long_jump_fr(T:k,(T-1));
    phi_loop = diag(normrnd(0,1,[1,k-T+1]));
    for i=1:(k-T+1)
        store_loop = 0;
        for j=(i+1):(k-T+1)
            store_loop = store_loop+fr0_loop(j)*tao*rho_loop(i,j)*sigma_loop(j)/(1+fr0_loop(j)*tao);
        end;
        mu_loop(i) = -sigma_loop(i)*store_loop;
    end;

    par_loop=(sum(b_loop*sqrt(dt)*phi_loop,2))';

    fr1_hat_loop = transpose(fr0_loop).*exp(mu_loop*tao-0.5*sigma_loop.^2*tao+sigma_loop.*par_loop);

    for i=1:(k-T+1)
        s_loop = 0;
        for j=(i+1):k
            s_loop = s_loop+fr1_hat(j)*tao*rho(i,j)*sigma(j)/(1+fr1_hat(j)*tao);
        end
        mu_hat_loop(i) = -sigma_loop(i)*s_loop;
    end;

    mu_new_loop = 0.5*(mu_loop+mu_hat_loop);
    
    % use the above outputs to generate F(T)
    fr1_loop = transpose(fr0_loop).*exp(mu_new_loop*tao-0.5*sigma_loop.^2*tao+sigma_loop.*par_loop);
    
    % append new F(T) to the long jump forward rates matrix
    long_jump_fr(T:k,T) = fr1_loop;
    
    clear rho_loop;
    clear mu_loop;
    clear mu_hat_loop;
    clear mu_new_loop;
end;
