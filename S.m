function [f] = S(iter,tk)

    global tao  k;
    
    s = 0;
    for j = (tk+1):k
        s = s + tao * P(iter,tk,j);
    end;
      
    f = (1 - P(iter,tk,k)) / s;
   
end