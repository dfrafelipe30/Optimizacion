function y = f(x)
    y = x^2 + x +1
endfunction
function Optimi = minSol2(f,a,b,e)
    gr = (sqrt(5) + 1)/2
    c = b - (b - a)/gr
    d = a + (b - a)/gr
    l = c - d 
    while (abs(l) >= e)
        if(f(c) < f(d)) 
            b = d
            d = c
            c = b - (b - a)/gr
        else
            a = c
            c = d
            d = a + (b - a)/gr
        end
        l = c - d   
    end
    Optimi = min(a,c,d,b)
endfunction 
