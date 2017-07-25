% read market data
expiry = xlsread('black_vol.xlsx','A:A');
forward_rate_list = xlsread('ForwardRates.xlsx','C:C'); % Market Forward Rates
% set the 'jumps'
n = length(expiry); % given product length
tao = 1/4; % given tao
dt = 1/4; % given long jump
m = 1/dt;
k = m*n; % total market cap vols needed
    
% calculate zero-coupon bond prices using market forward rates 
forward_rate = forward_rate_list(1:k)/100;

zcb(1) = 1;
for i = 2:(k+1)
    zcb(i) = zcb(i-1) / (1 + forward_rate(i-1) * tao);
end;

% calculate swap rates using market forward rates
cumulative_sum = 0;
for i = 1:(k-1)
    cumulative_sum = cumulative_sum + tao * zcb(i+1);
    swap_rate(i) = (zcb(2)-zcb(i+2)) / cumulative_sum;  
end;
    
% caplet pricing
principle = 1;
iteration = 10000;
caplet_output = [];
for sim = 1:iteration
    Monte_Carlo
    FRM = long_jump_fr; %FRM = csvread('Long_Jump_Forward_Rates.csv');
        
    for i= 1:(k-1)
        pro = 1;
        for j=(i+1):k
            pro = pro*(1+tao*FRM(j,i));
        end;
        caplet_value(i) = zcb(k+1)*principle*max((FRM(i,i)-swap_rate(i)),0)*pro;
    end;
    caplet_output = [caplet_output;caplet_value];
end;

caplet_price = sum(caplet_output)/iteration;

error_caplet_price_lmm = abs(caplet_price-caplet_market);
















