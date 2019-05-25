function y = f(x)
    //y = x(1) + 1.4*x(2)
    y = -x(1) - 1.4*x(2)
endfunction
function tmx = tmaxu(u,du)
    tmx = 1e10
    n = length(u)
    for i=1:n
        if du(i) < 0 then
            tmx = min(tmx,-u(i)/du(i))
        end
    end
endfunction
//x= [1; 1; 2; 2.8; 2]

x = [0.2;1.45;2.35;2.7;2.8]
//y = [-0.3;-0.6;-0.3]

y = [0;0;0]
//s = [0.2;0.1;0.3;0.6;0.3]

s = [1;1;1;1;1]

c = [-1 ; -1.4; 0 ; 0 ; 0]

A = [1 1 1 0 0;1 2 0 1 0;1 0 0 0 1]

b = [4;5.8;3]

function x_opt = M_P_D_A_F(A,b,c,x,y,s,epsilon,max_iter)
    for i=1:max_iter
        disp(i,"i")
        disp(x,"x")
        disp(y,"y")
        disp(s,"s")
        rd = A'*y + s - c 
        rp = A*x - b
        rc = x.*s
        disp(rc,"rc")
        r = [rd;rp;rc]
        roc = norm(rc)/(1 + norm(x))
        if roc < epsilon then
            break
        else
            n = size(A)(2)  //numero de columnas de A
            m = size(A)(1)
            primero = zeros(n,n)
            segundo = zeros(m,m)
            tercero = zeros(m,n)
            cuarto = zeros(n,m)
            e = eye(n,n)
            matriz = [primero A' e;A segundo tercero; diag(s) cuarto diag(x)]
            sol = -matriz\r
            dx = sol(1:5)
            dy = sol(6:8)
            ds = sol(9:13)
            disp(dx,"dx")
            disp(dy,"dy")
            disp(ds,"ds")
            tma_x = tmaxu(x,dx)
            tma_U = tmaxu(s,ds)
            tmax = min(tma_x,tma_U)
            tk = 0.99 * tmax
            disp(tmax,"tmax")
            disp(tk,"tk")
            x = x + tk*dx
            y = y + tk*dy
            s = s + tk*ds
        end
        x_opt = x
    end
endfunction

x_optimo = M_P_D_A_F(A,b,c,x,y,s,0.001,10)
disp(x_optimo,"X_optimo")
