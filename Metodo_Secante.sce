function y = f(x)
    y = x(1)^2 - x(2)^2 + 3*x(3)^3
endfunction
H = [2 0 0;0 -2 0; 0 0 3]
function minimo = Met_Gradiente(f,e,x,H)
    d = -(numderivative(f,x))'
    numerador = (-d' * d)
    denominador = d' * H * d 
    t = numerador/denominador
    xold = x
    x = x + t * d
    while abs(x - xold) > e
        d = -(numderivative(f,x))'
        numerador = (-d' * d)
        denominador = d' * H * d 
        t = numerador/denominador
        xold = x
        x = x + t*d
    end
    minimo = x
endfunction
