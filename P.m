function [f] = P(iter, t0, tn)

    global FRM  tao;
    
    m = 1;
    for i = t0:(tn-1)
        m = m / (1 + tao * FRM{iter}(i,t0));
    end;
      
    f = m;
   
end

