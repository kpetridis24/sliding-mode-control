function res = sat(x, e)
    if x > -e && x < e
        res = x / e;
    elseif x <= e
        res = -1;
    else 
        res = 1;  
    end
end
        