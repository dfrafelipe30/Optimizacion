function y = f(x)
    y = -1 * x(1) - 1.4*x(2)
endfunction
function [di,D] = delta(x)
    di = [1/(4 -x(1) -x(2));1/(5.8 - x(1) - 2*x(2));1/(3 - x(1));1/(x(1));1/(x(2))]
    D = diag(di)
endfunction

function tmax = phi(x,d,A,b)
    k = b - A*x
    h = A'*d
    tmax = 1
    for i = 1:size(h)(1)
        residuo = k(i)/h(i)
        if tmax < residuo
            tmax = residuo
        end
    end
endfunction
function tmx = tmax_g(g, x, d, M)
    //
    // tmx = max { t  :  g( x + t d ) <= 0,
    //             0 < t <= M }
    //
    // se supone que  g(x) < 0
    //
    //                0 < M
    //
    // utiliza biseccion (o dicotomia) y secante
    //
    //*************************
    EPS0 = 1.0e-8
    EPS  = 1.0e-6
    //*************************
    
    //escrVectED(x, 'x en tmax'), escrVectED(d, 'd')
    g0 = g(x)
    if max(g0) >= -EPS0
        printf('El punto no es estrictamente factible.\n')
        tmx = -1000
        return
    end
    gM = g(x + M*d)
    //escrVectED(gM, 'g(x+Md)')
    if max(gM) <= EPS0
        tmx = M
        return
    end
    
    a = 0
    b = M
    c0 = M/3
    c1 = 2*M/3
    gc0 = max( g( x + c0*d ) )
    gc1 = max( g( x + c1*d ) )
    
    while (b-a) >= EPS
    //escrVectED([a b], 'a b')
    deno = gc1 - gc0
    if abs(deno) < EPS0
        c2 = c1 + 0.01
    else
        c2 = c1 - gc1*(c1-c0)/deno
    end

    gc2 = max( g( x + c2*d ) )
    if 0 < c2 & c2 <= M & abs(gc2) <= EPS0
        //printf('secante\n')
        tmx = c2
        return
    end
    c0 = c1, gc0 = gc1
    c1 = c2, gc1 = gc2
    
    m = (a+b)/2
    gm = max( g(x + m*d) )
    
    if abs(gm ) <= EPS0
        tmx = m
        return
    end
    
    if gm > 0
        b = m
    else
        a = m
    end
    end
    tmx = (a+b)/2
endfunction // tmax_g

tao = 1
x = [1;1]
A = [1 1;1 2;1 0; -1 0;0 -1]
b= [4;5.8;3;0;0]
c = [-1 ; -1.4]
function y = g(x)
    y1 = 1 * x(1) + x(2) -4
    y2 = x(1) + 2*x(2) - 5.8
    y3 = x(1) - 3
    y4 = -x(1)
    y5 = -x(2)
    y = [y1;y2;y3;y4;y5]
endfunction
for k = 1:20
    disp(k,"k")
    [d,D] = delta(x)
    if(norm(d) < 0.001)
        break
    else
        B0 = A'* d
        B1 = A' * D**2 * A
        f0 = tao* c + B0
        f1 = B1 
        dire = -f1\f0
        disp(dire,"direccion")
        tmax = phi(x,d,A,b)
//        tmax = tmax_g(g,x,dire,1)
        disp(tmax,"tmax")
        tk = min(1,0.99*tmax)
        disp(tk,"tk")
        x = x + dire
        disp(x,"x")
    end
end     
