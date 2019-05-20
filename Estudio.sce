function y = f(x)
    y = 2 * x(1)**4 + (x(1) + x(2) - 1)**2 
endfunction 

x = [3;4]
d = [-5;- 1]
function escV(x,tit)
    n = length(x)
    printf('%s :',tit)
    for i = 1:n
        printf('%10.4f',x(i))
    end
    printf('\n')
endfunction

function [tmn,info] = tmejor(f,x,d,tk,eps,h_ini)
    h = h_ini
    phi_k = f(x+tk*d)
    escV(phi_k,"fk")
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
        f_x = f(x + (tk + h)*d)
        escV(f_x,"fk")
        h = h/2
    end
endfunction

t_m = tmejor(f,x,d,0,0.0010,10)
disp(t_m,"t_m")
