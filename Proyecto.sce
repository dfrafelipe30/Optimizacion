function fx = h(x)
    fx = (x(1)-1)**4 + (x(1)+2*x(2)-3)**2 + (x(1) + 3*x(2) - x(3) + x(4) - 2)**2 + 10
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
function topt = backtracking(f,x,d)
    alpha = 0.25
    bet = 0.5
    topt = 1
    gr = numderivative(f,x)
    gr = gr'
    while (f(x+t*d) > f(x) + alpha*topt*gr'*d)
        topt = bet*t
    end
endfunction

A = [1 2 3 4]
b = [10]
//x_m = [10;0;0;0]
x_m = [0;0;0;2.5]
B = [-2 -3 -4; 1 0 0;0 1 0;0 0 1]
z = [0;0;0]

Como calcular B,xm,z
function x = M_N_Q_I(f,A,b,Max,e) // Metodo de Newton sin desigualdades
    x_m = A\b
    B = kernel(A)
    z = zeros(1,size(B)(2))
    for i = 1:Max
         x = x_m + B*z
        [gr,H] = numderivative(f,x,[],[],"blockmat")
        gr = gr'
        if norm(gr) < e then 
            break
        else
            grz = B' * gr
            Hz = B' * H * B
            d = -Hz\grz
            if f(z  + d) < f(z)then
                z = z + d
            elseif gr * d < 0 then
                tk_t = tres_ptos(f,x,d,0,0.001,10,1e6)
                tk_b = backtracking(f,x,d)
                if tk_m == -1 then
                    z = z + tk_b * d
                else
                    z = z + tk_t * d
                end
            end
        end
    end
endfunction

