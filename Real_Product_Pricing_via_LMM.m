% % Real Product (Product 4)

N_rp = 1000;
frequency = 0.5;

jm_coupon_rate_lmm = zeros(k,1);
jm_coupon_rate_lmm(2:2:14) = 0.02;
jm_coupon_rate_lmm(16:2:18) = 0.0225;
jm_coupon_rate_lmm(20:2:22) = 0.025;
jm_coupon_rate_lmm(24:2:26) = 0.04;
jm_coupon_rate_lmm(28:2:32) = 0.06;

jm_coupon_lmm = N_rp * frequency * jm_coupon_rate_lmm;

credit_spread = 50 / 10000;              % obtained from 'http://www.boursorama.com/bourse/taux/cds-CDS_5A_JPMORGAN_CHASE_&_CO-3xJPM'

call_value_lmm = zeros(k,1);
call_value_lmm = N_rp + jm_coupon_lmm;
call_value_lmm(call_value_lmm == N_rp) = 0;

redeem_start_t = 4;

A = sum(jm_coupon_lmm((frequency/tao):(frequency/tao):((redeem_start_t-frequency)/tao)) .* ZCB(7:6:43));

iteration = 10000;

global FRM;
FRM = cell(iteration,1);
for iter = 1:iteration
    Monte_Carlo
    FRM{iter} = long_jump_fr; 
end;
for iter = 1:iteration
    for frm1 = 1:k
        for frm2 = 1:k
            if FRM{iter}(frm1,frm2) ~= 0
                FRM{iter}(frm1,frm2) = FRM{iter}(frm1,frm2) + credit_spread;
            end;
        end;
    end;
end;

Qi = zeros(iteration,k);
cont_value_lmm = zeros(iteration,k);
va = zeros(iteration,1);

for iter = 1:iteration
   Qi(iter,k) = call_value_lmm(k);
end;

for j = (k-frequency/tao):(-frequency/tao):(redeem_start_t/tao)
    X1 = zeros(iteration,1);
    X2 = zeros(iteration,1);
    
    for iter = 1:iteration
        cont_value_lmm(iter,j) = P(iter,j,k) * Qi(iter,j+frequency/tao) + jm_coupon_lmm(j);
        X1(iter) = FRM{iter}(j,j);
        X2(iter) = S(iter,j);
    end;
    
    Y = cont_value_lmm(:,j);
    ONE = ones(length(Y),1);
    X3 = X1.^2;
    X4 = X2.^2;
    X5 = X1.*X2;
    X = [ONE,X1,X2,X3,X4,X5];
    b_coefficient = regress(Y,X);
    Y_hat = (b_coefficient' * X')';
    
    for iter = 1:iteration
        if call_value_lmm(j) < Y_hat(iter,1)
            Qi(iter,j) = call_value_lmm(j) / P(iter,j,k);
        else
            Qi(iter,j) = cont_value_lmm(iter,j) / P(iter,j,k);
        end;
    end;
    
end;

for iter = 1:iteration
    va(iter) = Qi(iter,redeem_start_t/tao) * P(iter,1,k);
    va(iter) = va(iter) / (1 + tao*(forward_rate_list(1)/100));
end;

real_product_lmm = sum(va) / iteration + A;











