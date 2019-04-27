function y = f(x)
    y = x^2
endfunction

function Optimo = minSol(a,b,e,m,f)
    l = b-a
    while(l>=e)
        del = l/(m+1)
        u = [a:del:b]
        mini = a
        for i = a:del:b
            if(f(mini)>f(i))
                mini = i
            end
        end
        a = mini-del
        b = mini+del
        l = b - a
    end
    Optimo = mini
endfunction
