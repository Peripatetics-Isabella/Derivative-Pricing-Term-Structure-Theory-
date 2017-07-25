% % Project_BDT

% % Interval = 1M; Tenor = 3M
% global forward_rate;
% forward_rate = xlsread('forward_rate_30_360_1m_3m.xlsx', 'C2:C98'); % without '%'
% forward_rate = forward_rate / 100;
% M = length(forward_rate);

global tao 
tao = 0.25;

term = 8;

% Interval = 3M; Tenor = 3M
original_forward_rate  = xlsread('Forward Rate_10-27-16_30-360-3M-LIBOR_Product 4.xlsx', 'C:C'); % without '%'
global forward_rate;
forward_rate = original_forward_rate(2:term/tao+2); % without '%'
forward_rate = forward_rate / 100;
M = 97;

% Market Zero-Coupon Bond Price (ZCB)

global ZCB
ZCB = zeros(M,1);
ZCB(1) = 1;
for i = 2:length(forward_rate)
    ZCB(3*i-2) = ZCB(3*i-5) / (1 + forward_rate(i-1) * tao);
end;


% Interpolation for ZCB
for i = 2:3:95
    ZCB(i) = exp(2/3 * log(ZCB(i-1)) + 1/3 * log(ZCB(i+2)));
    ZCB(i+1) = exp(1/3 * log(ZCB(i-1)) + 2/3 * log(ZCB(i+2)));
end


% % Swap Rate / K
% swap_rate = zeros((M-1)/3-1,1);
% cumulative_sum = 0;
% 
% for i = 1:length(swap_rate)
%     cumulative_sum = cumulative_sum + tao * ZCB(3*i+4);
%     swap_rate(i) = (ZCB(4)-ZCB(3*i+4)) / cumulative_sum;  
% end;
% 
% K = swap_rate;


% % Swap Rate / K
% swap_rate = zeros(M-1,1);
% cumulative_sum = 0;
% 
% for i = 3:length(swap_rate)
%     cumulative_sum = cumulative_sum + tao * ZCB(i+1);
%     swap_rate(i) = (ZCB(4)-ZCB(i+1)) / cumulative_sum;  
% end;
% 
% global K
% K = zeros((M-1)/3-1,1);
% for i = 1:length(K)
%     K(i) = swap_rate(3*i+3);
% end;


% Swap Rate / K
swap_rate = zeros(length(forward_rate)-2,1);
cumulative_sum = 0;

for i = 1:length(swap_rate)
    cumulative_sum = cumulative_sum + tao * ZCB(3*i+4);
    swap_rate(i) = (ZCB(4)-ZCB(3*i+4)) / cumulative_sum;  
end;

global K
K = swap_rate;


% Cap Black Volatility
original_cap_vol = xlsread('cap_volatility_10.27.2016.xlsx', 'B:B');  % without '%'
cap_vol = original_cap_vol(2:term+1);  % without '%'
cap_vol = cap_vol / 100;


% Caplet Volatility with Interplation
global sigma
sigma = zeros(M,1);
for i = 1:length(cap_vol)
    sigma(12*i+1) = cap_vol(i);
end;

for i = 1:12
    sigma(i) = sigma(13);
end;

for i = 19:12:(M-6)
    sigma(i-3) = 0.75 * sigma(i-6) + 0.25 * sigma(i+6);
    sigma(i) = 0.5 * sigma(i-6) + 0.5 * sigma(i+6);
    sigma(i+3) = 0.25 * sigma(i-6) + 0.75 * sigma(i+6);
end;

for i = 14:3:(M-2)
    sigma(i) = sigma(i+2);
    sigma(i+1) = sigma(i+2);
end;


% Market Cap Price
N_caplet = 1;
global caplet;
caplet = zeros(length(K), length(K));
for i = 1:(length(K))
    for j = 1:i
        d1 = (log(forward_rate(j)/K(i)) + 0.5 * (sigma(3*i+1)^2) * (tao*(j))) / (sigma(3*i+1) * sqrt(tao*(j)));
        d2 = (log(forward_rate(j)/K(i)) - 0.5 * (sigma(3*i+1)^2) * (tao*(j))) / (sigma(3*i+1) * sqrt(tao*(j)));
        
        caplet(j,i) = N_caplet * tao * ZCB(3*j+4) * (forward_rate(j) * normcdf(d1) - K(i) * normcdf(d2));
        
    end;
end;

global cap;
cap = sum(caplet);


% Market Caplet Price (Adjusting Implied Volatility)
global previous_caplet;
previous_caplet = zeros(length(cap), length(cap));

global caplet_black;
caplet_black = zeros(1, length(cap));

iv = zeros(length(cap),1);
global c;

for i = 1:length(cap)
    c = i;
    
    % iv(i) = fminbnd(@caplet_black_price,0,1);
    % iv(i) = fsolve(@caplet_black_price,0.3);
    iv(i) = fminsearch(@caplet_black_price,0.5);
    
    for j = i:length(cap)
                
        d1 = (log(forward_rate(i)/K(j)) + 0.5 * (iv(i)^2) * (tao*(i))) / (iv(i) * sqrt(tao*(i)));
        d2 = (log(forward_rate(i)/K(j)) - 0.5 * (iv(i)^2) * (tao*(i))) / (iv(i) * sqrt(tao*(i)));
        
        previous_caplet(i,j) = N_caplet * tao * ZCB(3*i+4) * (forward_rate(i) * normcdf(d1) - K(j) * normcdf(d2)); 
       
    end;
        
    d1 = (log(forward_rate(i)/K(i)) + 0.5 * (iv(i)^2) * (tao*(i))) / (iv(i) * sqrt(tao*(i)));
    d2 = (log(forward_rate(i)/K(i)) - 0.5 * (iv(i)^2) * (tao*(i))) / (iv(i) * sqrt(tao*(i)));
        
    caplet_black(i) = N_caplet * tao * ZCB(3*i+4) * (forward_rate(i) * normcdf(d1) - K(i) * normcdf(d2));
    
end;

global caplet_market;
caplet_market = caplet_black;


% Sigma for BDT
sigma(1) = iv(1);
sigma(2) = iv(1);
sigma(3) = iv(1);
sigma(4) = iv(1);

for i = 1:length(iv)
    sigma(3*i+4) = iv(i);
    sigma(3*i+3) = iv(i);
    sigma(3*i+2) = iv(i);
end;


% Calibration to Zero-Coupon Bond Price
global U 
U = zeros(M,1);
global r 
r = zeros(M,M);
global d 
d = zeros(M,M);
global pi 
pi = zeros(M,M);

global dt 
dt = 1/12;
global sqrt_dt 
sqrt_dt = sqrt(dt);

pi(1,1) = 1;

U(1) = -log(ZCB(2))/dt;
r(1,1) = U(1);
d(1,1) = exp(-r(1,1)*dt);

for i = 2:(M-1)
    for j = 1:i
        if j == 1
            pi(i,j) = 0.5 * pi(i-1,j) * d(i-1,j);
        elseif j == i
            pi(i,j) = 0.5 * pi(i-1,j-1) * d(i-1,j-1);
        else
            pi(i,j) = 0.5 * pi(i-1,j) * d(i-1,j) + 0.5 * pi(i-1,j-1) * d(i-1,j-1);
        end;
    end;
    
    syms u;
    syms r_temp;
    syms d_temp;
    syms f;
    for j = 1:i
        r_temp(j) = u * exp((2*(j-1)-(i-1))*sigma(i)*sqrt_dt);
    end;
    
    for j = 1:i
        d_temp(j) = exp(-r_temp(j)*dt);
    end;
    
    f = sum((pi(i,1:i)).*d_temp)-ZCB(i+1);
    
    U(i) = solve(f);
    
    
    for j = 1:i
        r(i,j) = U(i) * exp((2*(j-1)-(i-1))*sigma(i)*sqrt_dt);
    end;

    for j = 1:i
        d(i,j) = exp(-r(i,j)*dt);
    end;
    
end;


% Value Caplets without adjusting volatilities
% First Caplet is (0,3M), payoff at 4M
caplet_tree_original = zeros(length(caplet_market)-1,1);
global N 
N = 1;
for i = 1:(length(caplet_market))
    payoff = N * 0.25 * max(r(3*i+1,1:3*i+1)-K(i),0);
    caplet_tree_original(i) = sum(payoff .* pi(3*i+1,1:3*i+1));
end;

error_caplet = abs(caplet_tree_original - caplet_market');


% Calibration to Caplet Price
global r_new
r_new = zeros(M,M);
global d_new 
d_new = zeros(M,M);
global pai 
pai = zeros(M,M);

pai(1,1) = 1;
r_new(1,1) = U(1);
d_new(1,1) = exp(-r_new(1,1)*dt);

sigma_new = zeros(M,1);
iv_bdt = zeros(length(caplet_market),1);

% Solve sigma_new for each Caplet
for l = 1:length(caplet_market)
    
    c = l;
    k = 3 * c + 4;

    % iv_bdt(c) = fminbnd(@bdt_caplet,0,1);
    % iv_bdt(c) = fminsearch(@bdt_caplet,0.5);
    iv_bdt(c) = fsolve(@bdt_caplet,0.3);
    
    if c == 1;
        a = 2;
    else
        a = k-3;
    end;
    
    for i = a:(k-1)
        sigma_new(i) = iv_bdt(c);
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
            r_new(i,j) = U(i) * exp((2*(j-1)-(i-1))*sigma_new(i)*sqrt_dt);
        end;
    
        for j = 1:i
            d_new(i,j) = exp(-r_new(i,j)*dt);
        end;
        
    end;
    
end;


% Error of ZCB
ZCB_new = sum(pai,2);               % Calculate the sum of each row to get the ZCB_new
error_ZCB = abs(ZCB - ZCB_new);


% % European Swaption

% Black price for European Swaption
% European Swaption ((1y,2y,3y,4y,5y,6y,7y)--8y)

sigma_black = zeros(7,1);
sigma_black(1) = xlsread('volatility_ATM_swaptions_10.27.2016.xlsx', 'Quotes', 'H6') / 100;    % Black volatility (in %) obtained from Bloomberg Terminal
sigma_black(2) = xlsread('volatility_ATM_swaptions_10.27.2016.xlsx', 'Quotes', 'G7') / 100;
sigma_black(3) = xlsread('volatility_ATM_swaptions_10.27.2016.xlsx', 'Quotes', 'F8') / 100;
sigma_black(4) = xlsread('volatility_ATM_swaptions_10.27.2016.xlsx', 'Quotes', 'E9') / 100;
sigma_black(5) = xlsread('volatility_ATM_swaptions_10.27.2016.xlsx', 'Quotes', 'D10') / 100;
sigma_black(6) = xlsread('volatility_ATM_swaptions_10.27.2016.xlsx', 'Quotes', 'C11') / 100;
sigma_black(7) = xlsread('volatility_ATM_swaptions_10.27.2016.xlsx', 'Quotes', 'B12') / 100;

T0 = [13,25,37,49,61,73,85];
Tn = 8 * 12 + 1;




% sigma_black = zeros(6,1);
% sigma_black(1) = xlsread('volatility_ATM_swaptions_10.27.2016.xlsx', 'Quotes', 'G6') / 100;
% sigma_black(2) = xlsread('volatility_ATM_swaptions_10.27.2016.xlsx', 'Quotes', 'F7') / 100;
% sigma_black(3) = xlsread('volatility_ATM_swaptions_10.27.2016.xlsx', 'Quotes', 'E8') / 100;
% sigma_black(4) = xlsread('volatility_ATM_swaptions_10.27.2016.xlsx', 'Quotes', 'D9') / 100;
% sigma_black(5) = xlsread('volatility_ATM_swaptions_10.27.2016.xlsx', 'Quotes', 'C10') / 100;
% sigma_black(6) = xlsread('volatility_ATM_swaptions_10.27.2016.xlsx', 'Quotes', 'B11') / 100;
% 
% T0 = [13,25,37,49,61,73];
% Tn = 7 * 12 + 1;

principal = 1;
tao = 0.25;


swap_rate = zeros(length(sigma_black),1);
strike_euro = swap_rate;
swaption_black = zeros(length(sigma_black),1);

for l = 1:length(sigma_black)
    
    s = 0;
    for i = (T0(l)+3):3:Tn
        s = s + tao * ZCB(i);
    end;

    swap_rate(l) = (ZCB(T0(l))-ZCB(Tn)) / s;
    strike_euro(l) = swap_rate(l);         % ATM

    T = (T0(l)-1)/12;

    d1 = (log(swap_rate(l)/strike_euro(l)) + 0.5 * sigma_black(l) * sigma_black(l) * T) / (sigma_black(l)*sqrt(T));
    d2 = (log(swap_rate(l)/strike_euro(l)) - 0.5 * sigma_black(l) * sigma_black(l) * T) / (sigma_black(l)*sqrt(T));

    swaption_black(l) = principal * s * (swap_rate(l) * normcdf(d1) - strike_euro(l) * normcdf(d2));
    
end;

% BDT Price for European Swaption

swaption_BDT = zeros(length(sigma_black),1);
coupon = principal * strike_euro * tao;

for l = 1:length(sigma_black)
    v = zeros(Tn,Tn);
    v(Tn,:) = principal + coupon(l);

    for i = (Tn-1):(-1):T0(l)
        if mod((i-1),3) == 0
            for j = 1:i
                v(i,j) = (0.5 * v(i+1,j) + 0.5 * v(i+1,j+1)) * d_new(i,j) + coupon(l);
            end;
        else
            for j = 1:i
                v(i,j) = (0.5 * v(i+1,j) + 0.5 * v(i+1,j+1)) * d_new(i,j);
            end;
        end;
    end;

    payment = max(principal-(v(T0(l),1:T0(l))-coupon(l)), 0);
    swaption_BDT(l) = sum(pai(T0(l),1:T0(l)) .* payment);
end;

% Error of European Swaption

error_european_swaption_BDT = abs(swaption_BDT-swaption_black);
    

% % Bermudan Swaption (2-year Bermudan option on a 5-year receiver swap)

% Coupon rate = swap_rate(T0 = 2y, Tn = 7y)
N = 1;
T0_b = 2 * 12 + 1;
Tn_b = 7 * 12 + 1;
tenor = Tn_b - T0_b;
tao = 0.25;

s = 0;
for i = (T0_b+3):3:Tn_b
    s = s + tao * ZCB(i);
end;
strike = (ZCB(T0_b)-ZCB(Tn_b)) / s;
coupon_b = N * strike * tao;

% Coupon Bond

cb = zeros(Tn_b);
cb(Tn_b,:) = N + coupon_b;

coupon_bond = zeros(T0_b,Tn_b);

for i = (Tn_b-1):(-1):1
    if mod((i-1),3) == 0 && i ~= 1
        for j = 1:i
            cb(i,j) = (0.5 * cb(i+1,j) + 0.5 * cb(i+1,j+1)) * d_new(i,j) + coupon_b;
        end;
    else
        for j = 1:i
            cb(i,j) = (0.5 * cb(i+1,j) + 0.5 * cb(i+1,j+1)) * d_new(i,j);
        end;
    end;
end;

for i = (tao*12+1):(tao*12):T0_b
    coupon_bond(i,:) = cb(i,:) - coupon_b;
end;
coupon_bond(coupon_bond==-coupon_b) = 0;

% % 5-year Coupon Bonds
% coupon_bond = zeros(T0_b,T0_b);
% 
% for l = 1:((T0_b-1)/3)
%     t0 = l*tao/dt+1;
%     tm = t0 + tenor;
%     cb = zeros(tm);
%     cb(tm,:) = N + coupon_b;
%     for i = (tm-1):(-1):t0
%         if mod((i-1),3) == 0
%             for j = 1:i
%                 cb(i,j) = (0.5 * cb(i+1,j) + 0.5 * cb(i+1,j+1)) * d_new(i,j) + coupon_b;
%             end;
%         else
%             for j = 1:i
%                 cb(i,j) = (0.5 * cb(i+1,j) + 0.5 * cb(i+1,j+1)) * d_new(i,j);
%             end;
%         end;
%     end;
%     
%     coupon_bond(t0,1:t0) = cb(t0,1:t0)-coupon_b;
%     
% end;

% Early Exercise Value (EE)
ee = coupon_bond;
ee = ee - N;
ee(ee==-N) = 0;

% Bermudan Swaption
bs = zeros(Tn_b);
bs(T0_b,:) = max(ee(T0_b,:),0);

for i = (T0_b-1):(-1):1
    if mod((i-1),3) == 0 && i ~= 1
        for j = 1:i
            bs(i,j) = max(((0.5 * bs(i+1,j) + 0.5 * bs(i+1,j+1)) * d_new(i,j)),ee(i,j));
        end;
    else
        for j = 1:i
            bs(i,j) = (0.5 * bs(i+1,j) + 0.5 * bs(i+1,j+1)) * d_new(i,j);
        end;
    end;
end;

bermudan_swaption = bs(1,1);


% % Real Product (Product 4)

% Add credit spread to BDT rates
% tr_rate = xlsread('treasury_note_rate.xlsx', 'E2:J2'); % without '%'
% tr_rate = tr_rate / 100;
% treasury_rate = zeros(10,1);
% treasury_rate(1) = tr_rate(1);
% treasury_rate(2) = tr_rate(2);
% treasury_rate(3) = tr_rate(3);
% treasury_rate(5) = tr_rate(4);
% treasury_rate(7) = tr_rate(5);
% treasury_rate(10) = tr_rate(6);
% 
% treasury_rate(4) = 0.5 * treasury_rate(3) + 0.5 * treasury_rate(5);
% treasury_rate(6) = 0.5 * treasury_rate(5) + 0.5 * treasury_rate(7);
% treasury_rate(8) = 2/3 * treasury_rate(7) + 1/3 * treasury_rate(10);
% treasury_rate(9) = 1/3 * treasury_rate(7) + 2/3 * treasury_rate(10);
 
jm_coupon_rate = zeros(8,1);
jm_coupon_rate(1:4) = 0.02;
jm_coupon_rate(5) = 0.0225;
jm_coupon_rate(6) = 0.025;
jm_coupon_rate(7) = 0.04;
jm_coupon_rate(8) = 0.06;

% credit_spread = jm_coupon_rate - treasury_rate(1:8);

credit_spread = 50 / 10000;              % obtained from 'http://www.boursorama.com/bourse/taux/cds-CDS_5A_JPMORGAN_CHASE_&_CO-3xJPM'

r_cs = r_new;

% for i = 1:length(r_cs)
%     for j = 1:length(credit_spread)
%         if (i-1)/12 >= (j-1) && (i-1)/12 < j
%             r_cs(i,1:i) = r_cs(i,1:i) + credit_spread(j);
%             
%         end; 
%     end;    
% end;

r_cs = r_cs + credit_spread;
d_cs = exp(-r_cs*dt);


% % Without adding credit spread
% r_cs = r_new;
% d_cs = exp(-r_cs*dt);


% Valuation of product 4 through BDT tree

N_rp = 1000;
frequency = 0.5;
redeem_start_t = 4;
tn = length(r_cs);

jm_coupon = N_rp * frequency * jm_coupon_rate;

call_value = zeros(tn);

for i = (redeem_start_t*12+1):6:tn
    x = (i-1)/12;
    for j = 1:length(jm_coupon)
        if x >= (j-1) && x < j
            call_value(i,1:i) = N_rp + jm_coupon(j);
        end;
    end;
end;

rp = zeros(tn);
rp(tn,:) = N_rp + jm_coupon(length(jm_coupon));

rp_coupon = zeros(tn);
for i = 7:6:(tn-1)
    rp_coupon(i,1:i) = jm_coupon(floor((i-1)/12)+1);
end;

for i = (tn-1):(-1):1
    for j = 1:i
        if i >= (redeem_start_t*12+1) && mod((i-1),6) == 0
            rp(i,j) = min((0.5*rp(i+1,j)+0.5*rp(i+1,j+1))*d_cs(i,j)+rp_coupon(i,j),call_value(i,j));
        elseif mod((i-1),6) == 0 && i ~= 1
            rp(i,j) = (0.5*rp(i+1,j)+0.5*rp(i+1,j+1))*d_cs(i,j)+rp_coupon(i,j);
        else
            rp(i,j) = (0.5*rp(i+1,j)+0.5*rp(i+1,j+1))*d_cs(i,j);
        end;
    end;
end;

real_product_BDT = rp(1,1);








