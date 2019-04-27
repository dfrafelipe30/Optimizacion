function y = f(x)
    y = x(1)^2 + x(2)^2
endfunction
function Jf = Jacobiano(f,x)
    Jf = [2*x(1) 0; 0  2*x(2)]
endfunction
function Optimo = MNR(x0,f,e)
    Jf = Jacobiano(f,x0)
    x = x0 - inv(Jf)*f(x0)
    while (norm(x - x0) > e)
         Jf = Jacobiano(f,x)
         x0 = x
         x = x - inv(Jf) * f(x)
    end
    Optimo = xendfunction
 endfunction
