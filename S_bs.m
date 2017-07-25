function [f] = S_bs(iter,tk)

    global tao;
    
    TM = 7/tao;
    s = 0;
    for j = (tk+1):TM
        s = s + tao * P(iter,tk,j);
    end;
      
    f = (1 - P(iter,tk,TM)) / s;
   
end