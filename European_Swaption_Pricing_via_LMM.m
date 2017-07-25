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
clear swap_rate;
cumulative_sum = 0;
for i = 1:(k-1)
    cumulative_sum = cumulative_sum + tao * zcb(i+1);
    swap_rate(i) = (zcb(2)-zcb(i+2)) / cumulative_sum;  
end;

% European Swaption pricing
principle = 1;
iteration = 5000;
maturity_of_underlying_swap = [1,2,3,4,5,6,7]; % modify here to set maturity of underlying swap
maturity_of_swaption = 8; % modify here to set maturity of swaption

start_point = maturity_of_underlying_swap/tao;
stop_point = maturity_of_swaption/tao;

for sim = 1:iteration
    Monte_Carlo
    FRM = long_jump_fr; %FRM = csvread('Long_Jump_Forward_Rates.csv');
    
    for im = 1:length(maturity_of_underlying_swap)
        % value the underlying swap 
        swap_pro = 1;
        swap_payment = zeros(1,stop_point-start_point(im)+1);
        for i = 1:(stop_point-start_point(im)+1)
            swap_pro = swap_pro/(1+FRM((i+start_point(im)-1),start_point(im))*tao);
            swap_payment(i) = swap_pro;
        end;
        
        European_Swaption_value(sim,im) = zcb(stop_point+1)*principle*max((1-swap_payment(stop_point-start_point(im)+1)-strike_euro(start_point(im)*tao)*tao*sum(swap_payment)),0)/swap_payment(stop_point-start_point(im)+1);
        
    end;
end;

European_Swaption_price = sum(European_Swaption_value)/iteration;
European_Swaption_price = European_Swaption_price';

error_european_swaption_LMM = abs(European_Swaption_price - swaption_black);














