% % Bermudan Swaption (2y-5y) TM = 7, which is fixed.

N = 1;
iteration = 10000;

global FRM;
FRM = cell(iteration,1);


for iter = 1:iteration
    Monte_Carlo
    FRM{iter} = long_jump_fr; 
end;

tao = 0.25;
T0 = 0.25*(1/tao);
Tn = 2*(1/tao);
TM = 7*(1/tao);

bs_lmm = zeros(iteration,Tn);
Q = zeros(iteration,Tn);
ee_lmm = zeros(iteration,Tn);
cont_value_lmm = zeros(iteration,Tn);


for iter = 1:iteration
    s = 0;
    for j = Tn:(TM-1)
        s = s + N * tao * (FRM{iter}(j,Tn)-strike) * P(iter,Tn,j+1);
    end;
    bs_lmm(iter,Tn) = max(s,0);
    Q(iter,Tn) = bs_lmm(iter,Tn) / P(iter,Tn,TM);
end;


for ti = (Tn-1):(-1):T0
    X1 = zeros(iteration,1);
    X2 = zeros(iteration,1);
    
    for iter = 1:iteration
        s = 0;
        for j = (ti+1):(TM-1)
            s = s + N * tao * (FRM{iter}(j,ti)-strike) * P(iter,ti,j+1);
        end;
        ee_lmm(iter,ti) = s;
        cont_value_lmm(iter,ti) = Q(iter,ti+1) * P(iter,ti,TM);
        
        X1(iter) = FRM{iter}(ti,ti);
        X2(iter) = S_bs(iter,ti);
    end;
    
    Y = cont_value_lmm(:,ti);
    ONE = ones(length(Y),1);
    X3 = X1.^2;
    X4 = X2.^2;
    X5 = X1.*X2;
    X = [ONE,X1,X2,X3,X4,X5];
    b_coefficient = regress(Y,X);
    Y_hat = (b_coefficient' * X')';
    
    for iter = 1:iteration
        if ee_lmm(iter,ti) > Y_hat(iter,1)
            bs_lmm(iter,ti) = ee_lmm(iter,ti);
        else
            bs_lmm(iter,ti) = cont_value_lmm(iter,ti);
        end;
    
        Q(iter,ti) = bs_lmm(iter,ti) / P(iter,ti,TM);
    
    end;
    
end;

bermudan_swaption_lmm = max(1/iteration * ZCB(TM*3+1) * sum(Q(:,1)), N * (1 - ZCB(TM*3+1) - strike * tao * sum(ZCB(4:3:(TM*3+1)))));












