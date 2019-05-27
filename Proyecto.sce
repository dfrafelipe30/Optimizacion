/*function y = f(x)
    y = x(1)**4 - 11*x(1)**3 + 41*x(1)**2 - 61*x(1) + 30 + (x(1) - x(2))**2 + (x(2) + x(3) - 1)**2
endfunction*/
/*function y = f(x)
    y = x(1)**2 + x(2)**2 + 5*x(3)**4
endfunction*/
function y = f(x)
    y = x(1)**2 + x(2) + 5*x(3)
endfunction
//A1 =[1 1 0]
//b1 =[5]
//A2 =[1 1 0 ;1 2 1]
//b2 =[5;6]
//A3 =[1 1 1]
//b3 =[10]
//A4 =[1 1 1;1 2 -1]
//b4 =[15;10]
//A5 =[0 1 1]
//b5 =[10]

function sol = fz(z,x_m,B,f)
    x = x_m + B*z
    sol = f(x)
endfunction
function [tmn,info] = tmejor(f,z,d,tk,eps,h_ini,fz,B,x_m)
    h = h_ini
    phi_k = fz(z+tk*d,x_m,B,f)
    info = 0
    tmn = tk
    while(abs(h) > eps)
        if (fz(z+(tk+h)*d,x_m,B,f) < phi_k)
            tmn = tk + h
            info = 1
            break
        else
            if (fz(z+(tk-h)*d,x_m,B,f) < phi_k)
                tmn = tk - h
                info = 1
                break
            end
        end
        h = h/2
    end

endfunction
function [tao1,tao2,tao3,phi_1,phi_2,phi_3] = busq_tres_ptos(f,z,d,tk,t_m,K,fz,B,x_m)
    tao1 = tk
    phi_1 = fz(z+tao1*d,x_m,B,f)
    tao2 = t_m
    phi_2 = fz(z+tao2*d,x_m,B,f)
    tao3 = tao2
    while (norm(tao3) < K)
        tao3 = tao2 + 2*(tao2 - tao1)
        phi_3 = fz(z+tao3*d,x_m,B,f)
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
//    printf('Denominador nulo.\n')
    tp = 0
    return
    end
    tp = t2 - 0.5*(t21*c1 - t23*c2)/deno
endfunction

function t_opt = tres_ptos(f,z,d,tk,eps,h_ini,K,fz,B,x_m)
    while (1) //while true
    [t_m,info] = tmejor(f,z,d,tk,eps,h_ini,fz,B,x_m)
//    disp(info,"info:")
    if (info == 0)
        t_opt = -1
//        disp("Parece que no hay minimizador")
        break
    else
        [tao1,tao2,tao3,phi_1,phi_2,phi_3] = busq_tres_ptos(f,z,d,tk,t_m,K,fz,B,x_m)
        tp = interpCuadr(tao1,tao2,tao3,phi_1,phi_2,phi_3)
        if (fz(z+t_m*d,x_m,B,f) < fz(z+tp*d,x_m,B,f))
            t_opt = t_m
        else
            t_opt = tp
        end
//        disp(t_opt,"toptimo")
        tk = t_opt
    end
    
    end
endfunction
function topt = backtracking(f,z,d,fz,B,x_m)
    alpha = 0.25
    bet = 0.5
    topt = 1
    gr = numderivative(fz,z)
    gr = gr'
    while (fz(z+topt*d,x_m,B,f) > fz(z,x_m,B,f) + alpha*topt*gr'*d)
        topt = bet*topt
    end
endfunction

function x = M_N_Q_I(f,A,b,Max,e) // Metodo de Newton sin desigualdades
    x_m = A\b
    B = kernel(A)
    z = zeros(size(B)(2),1)
    estado = 1
    for i = 1:Max
//        disp(i,"i")
        x = x_m + B*z
        [gr,H] = numderivative(f,x,[],[],"blockmat")
        gr = gr'
        grz = B' * gr
        Hz = B' * H * B
//        disp(norm(grz),"norm_grz")
        if norm(grz) < e then 
            estado = 0
            break
        else
            d = -Hz\grz
//            disp(grz,"grz")
//            disp(Hz,"Hz")
            if f(x_m+B*(z + d)) < f(x_m + B*z)then
                z = z + d
//            disp(z,"z")
        elseif grz'*d < 0 then
                tk_t = tres_ptos(f,z,d,0,0.001,10,1e6,fz,B,x_m)
                tk_b = backtracking(f,z,d)
                if tk_t == -1 then
                    z = z + tk_b * d
                else
                    z = z + tk_t * d
                end
            end
        end
    end
    if(estado) then
        disp "No hay minimizador"
    end
endfunction

xsol = M_N_Q_I(f,A5,b5,10,0.00001)
disp(xsol)
