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
//    disp(info,"info:")
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
//        disp(t_opt,"toptimo")
        tk = t_opt
    end
end
endfunction

function escV(x,tit)
    n = length(x)
    printf('%s :',tit)
    for i = 1:n
        printf('%10.4f',x(i))
    end
    printf('\n')
endfunction

function sigma = Sigma(A)
    temp = []
    B = abs(A)
    for i = 1:size(A)(1)
        sig_i = sum(B(i,:)) - B(i,i)
//        disp(sig_i,"sigma i")
        temp(i) = sig_i
    end
    temp2 = []
    for j = 1:size(temp)(2)
        temp2(j)= temp(1) - A(j,j)
//        disp(temp2(j),"sigma")
    end
    sigma_1 = max(temp2)
//    disp(sigma_1,"sigma")
    sigma = max(sigma_1,0)
endfunction

function x_opt = M_N_M(x,f,epsilon,max_iter)
    for i =1:max_iter
        disp(i,"i")
        escV(x,"xk")
        [gr,H] = numderivative(f,x,[],[],'blockmat')
        escV(f(x),"f(X)")
        escV(gr,"gr")
        disp(H,"H")
        tk = 0
        if norm(gr) < epsilon then
            break
        else
            lamda = 0
            fink = 0
            while(fink == 0)
                //disp("Entro al while")
                M = H + lamda * eye(H)
                d =-M\gr'
                disp(d,"d")
                disp(x +d,"x +d")
                if f(x + d) < f(x) then
                    disp("entro")
                    x = x + d
                    disp(x,"x")
                    fink = 1
                else
                    disp(gr*d,"gr*d")
                    if gr*d < 0 then
                        disp("entro2")
                        escV(gr*d,"direccion")
                        t = tres_ptos(f,x,d,tk,0.001,10,1e6)
                        x = x + t *d
                        disp(x,"x")
                        fink = 1
                     else
                         lamda = Sigma(H)
                         disp(lamda,"lambda")
                    end
                end 
            end
        end
    end
    x_opt = x
endfunction

function y = f(x)
    y = (x(1)+1.5)*(x(1)+0.5)*(x(1)-0.5) + (x(2)-0.5)*(x(2)-1.5)*(x(2)-2.5) + 0.3*x(1)*x(2) + 0.01*(x(1)-3)**4 + 0.01*(x(2) - 4 )**4
endfunction

x = [-2.5;1]

x_optimo = M_N_M(x,f,0.001,20)
disp(x_optimo,"x optimo")


/*[gr,H] = numderivative(f,x,[],[],'blockmat')
M = H + (7.8818182)* eye(H)
d = -M\gr'*/
