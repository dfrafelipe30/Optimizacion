function y = f(x)
    y = 2* x(1)**4 + (x(1) + x(2) - 1)**2
endfunction

function Optimo = minSol(a,b,e,m,f,x,d)
    dif = b - a
    while(dif >= e)
        h = dif/(m+1)
        u = [a:h:b]
        mini = a + h
        for i = a:h:m 
            t_i = a + i *h
            if f(x + mini*d) > f(x + t_i* d) then
                mini = a + i*h
            end
        end
        a = mini-h
        b = mini+h
        dif = b - a
    end
    Optimo = mini
endfunction

x = [3;4]
d = [-5;-1]
t = minSol(0,1,0.01,2,f,x,d)
disp(t,"t")
