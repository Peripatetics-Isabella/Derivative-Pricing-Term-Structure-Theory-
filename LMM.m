% read market data
expiry = xlsread('black_vol.xlsx','A:A');
market_vol = xlsread('black_vol.xlsx','B:B'); % cap volatility
forward_rate_list = xlsread('ForwardRates.xlsx','C:C'); % Market Forward Rates

% set the 'jumps'
n = length(expiry); % given product length
tao = 1/4; % given tao
dt = 1/4; % given long jump
m = 1/dt;
global k;
k = m*n; % total market cap vols needed

% use linear interpolation to generate cap vols needed for generation of
% piecewise-constant volatilities
global black_vol;
for i=1:k
   if i<=m
       black_vol(i) = market_vol(1);
   else
       low = floor(i/m);
       up = ceil(i/m);
       black_vol(i) = (i/m-low)*(market_vol(up)-market_vol(low))+market_vol(low);
   end;  
end;
black_vol = black_vol/100;


% Calibration

% swap_rate_lmm
global w;
w = zeros(8);
for tn = 2:8
    for t0 = 1:(tn-1)
        denominator = 0;
        for j = (t0+tao):tao:tn;
            denominator = denominator + tao * ZCB(12*j);
        end;
        w(t0,tn) = tao * ZCB(12*t0) /denominator;
    end;
end;

global swap_rate_lmm;
swap_rate_lmm = zeros(8);
for tn = 2:8
    for t0 = 1:(tn-1)
        for l = t0:tn
            swap_rate_lmm(t0,tn) = swap_rate_lmm(t0,tn) + w(l,tn) * forward_rate(4*l);
        end;
    end;
end;



% black_vol(2:k) = iv;
% black_vol = iv';
global rho swaption_volatility;
swaption_volatility = xlsread('swaptions black vols_10-27-16_Product 4.xlsx','B23:I30')/100;

sigma = zeros(1,k);
rho = zeros(k);


% % function 2 for caplet volatility
% yeta = 0.11:0.01:0.42;
% alpha = 0.05;
% beta = 0.1;
% gamma = 0.05;
% 
% x = zeros(k,4);
% x(1,1) = alpha;
% x(1,2) = beta;
% x(1,3) = gamma;
% x(:,4) = yeta;
% 
% solution = fsolve(@lmm_calibration, x);
% lmm_calibration(solution);


% function 7 for caplet volatility
alpha = 0.05;
beta = 0.1;
gamma = 0.05;
a = 1;
b7 = 1;
c = 1;
d = 1;
fai_square = 0.11:0.01:0.42;
x = zeros(k,8);
x(1,1) = alpha;
x(1,2) = beta;
x(1,3) = gamma;
x(1,4) = a;
x(1,5) = b7;
x(1,6) = c;
x(1,7) = d;
x(:,8) = fai_square;

solution = fsolve(@lmm_calibration_function_7, x);

global swaption_sigma
swaption_sigma = [];

lmm_calibration_function_7(solution);


% Test Rebunato's Accuracy
error_swaption_volatility = abs(swaption_sigma - swaption_volatility);


%generate square root of corrlation matrix
global b;
b = sqrtm(rho);

% % generate piecewise-constant volatilities (here we use the easy way
% % because tao is constant)
% 
% % sigam = zeros(1,k);
% global sigma;
% sigma(1) = black_vol(1);
% for i=2:k
%     sigma(i)=sqrt(i*tao*black_vol(i)^2-(i-1)*tao*black_vol(i-1)^2);
% end;

% % generate corraltion matrix
% % here, tried the fifth method, but during the loop, it got stuck for the
% % last 2 steps, mainly due to the limitation of the model and we don't
% % wanna hard-code the last 2 correlation matrice.
% alpha = 0.05;
% beta = 0.1;
% gamma = 0.05;
% for i=1:k
%     for j=1:k
%         rho(i,j) = alpha + (1-alpha)*exp(-abs(i-j)*tao*beta*exp(gamma*tao*min(i,j)));
%     end;
% end;



%generate drift term
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

% output the long jump forward rates matrix
csvwrite('Long_Jump_Forward_Rates.csv',long_jump_fr);






