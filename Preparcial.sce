function Optimo = minsol(a,b,e,m,f)
    l = b - a
    m = 2               //numero de subdiviciones.
    x = a:(a-b)/m + 1:b        //Intervalo subdividido
    for i = 1:1:size(x,2)
        fval = f(x(i));
    end
    iter = 0;
    while (l >= e)
        iter = iter + 1;
        D(iter,:) = [a,b];
        [fun,pos] = min(fval)
        
        if pos ~= 1
            a = x(pos - 1);
            fval(1,1) = fval(pos - 1);
        else
            a = x(pos);
            fval(1,1) = fval(pos);
        end
        
        if pos ~= size(x,2)
            b = x(pos + 1);
            fval(1,end) = fval(pos +1);
        else
            a = x(pos);
            fval(1,end) = fval(pos);
        end
        
        x = a:(a-b)/(m + 1):b;
        
        for i = 2:1:size(x,2) - 1
            fval(1,i) = f(x(i));
        end 
    end
    Optimo = min(x);
endfunction

