function y = f(x)
    y = 0.1*x(1)^2 + x(2)**2 + 10*x(3)+100*x(4)**2 - 0.2*x(1) - 2*x(2) - 20*x(3) - 200*x(4)
endfunction

function [tmn,info] = tmejor(f,x,d,tk,eps,h_ini)
    h = h_ini
    phi_k = f(x+tk*d)
    info = 0
    tmn = 1000000
    while(abs(h) > eps)
        if (f(x+(tk+h)*d) < phi_k)
            tmn = tk + h
            info = 1
            break
        else
            if (f(x+(tk-h)*d) < phi_k)
                tmn = tk - h
                info = 1
                break
            end
        end
        h = h/2
    end

endfunction
function [tao1,tao2,tao3,phi_1,phi_2,phi_3] = busq_tres_ptos(f,x,d,tk,t_m,K)
    tao1 = tk
    phi_1 = f(x+tao1*d)
    tao2 = t_m
    phi_2 = f(x+tao2*d)
    tao3 = tao2
    while (norm(tao3) < K)
        tao3 = tao2 + 2*(tao2 - tao1)
        phi_3 = f(x+tao3*d)
        if (phi_3 > phi_2) 
            break
        else
            tao1 = tao2
            tao2 = tao3
            phi_1 = phi_2
            phi_2 = phi_3
        end
            
    end
endfunction

function tp = interpCuadr(t1,t2,t3,f1,f2,f3)
t21 = t2 - t1
t23 = t2 - t3
c1 = t21*(f2 - f3)
c2 = t23*(f2 - f1)
deno = c1 - c2

if abs(deno) < 1e-8
printf('Denominador nulo.\n')
tp = 0
return
end
tp = t2 - 0.5*(t21*c1 - t23*c2)/deno
endfunction

function t_opt = tres_ptos(f,x,d,tk,eps,h_ini,K)
    while (1) //while true
    [t_m,info] = tmejor(f,x,d,tk,eps,h_ini)
    disp(info,"info:")
    if (info == 0)
        t_opt = -1
        disp("Parece que no hay minimizador")
        break
    else
        [tao1,tao2,tao3,phi_1,phi_2,phi_3] = busq_tres_ptos(f,x,d,tk,t_m,K)
        tp = interpCuadr(tao1,tao2,tao3,phi_1,phi_2,phi_3)
        if (f(x+t_m*d) < f(x+tp*d))
            t_opt = t_m
        else
            t_opt = tp
        end
        disp(t_opt,"toptimo")
        tk = t_opt
    end
end
endfunction

function optimo = M_T_P(x,epsilon,max_it)
    for i = 1:max_it
        y1 = x
        gr = numderivative(f,y1)
        if gr <epsilon then
            break
        else
            d1 = -gr
            t1 = tres_ptos(f,y1,d,0,epsilon,10,1e6)
            y2 = y1+ t1 * d1 
            for j=2:
        end
    end
endfunction
