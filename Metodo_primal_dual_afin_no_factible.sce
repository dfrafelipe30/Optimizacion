function y = f(x)
    //y = x(1) + 1.4*x(2)
    y = -x(1) - 1.4*x(2) //Problema dual
endfunction
function tmx = tmaxX(G,q,x,dx)
    tmx = 1e10
    n = length(x)
    gx = G * x - q
    gd = G *dx
    for i = 1:n
        if gd(i) > 0 then
            tmx = min(tmx,-gx(i)/gd(i))
        end
    end
endfunction

x= [1; 1; 2; 2.8; 2]

y = [-0.3;-0.6;-0.3]

s = [1;1;1;1;1]

c = [-1 ; -1.4; 0 ; 0 ; 0]

A = [1 1 1 0 0;1 2 0 1 0;1 0 0 0 1]

b = [4;5.8;3]

primero = zeros(5,5)
segundo = zeros(3,3)
tercero = zeros(3,5)
cuarto = zeros(5,3)
e = eye(5,5)

rd = A'*y + s - c 
rp = A*x - b
rc = x.*s

r = [rd;rp;rc]

//disp (rd,"rd")
//disp(rp,"rp")
//disp(rc,"rc")

matriz = [primero A' e;A segundo tercero; diag(s) cuarto diag(x)]
//disp(matriz,"matriz")

sol = -matriz\r

//disp(sol)

dx = sol(1:5)
dy = sol(6:8)
ds = sol(9:13)

//disp(dx,"dx")
//disp(dy,"dy")
//disp(ds,"ds")

tma_x = tmaxu(x,dx)
tma_U = tmaxu(s,ds)
tmax = min(tma_x,tma_U)

tk = 0.99 * tmax
disp(tk)

x = x + tk*dx
y = y + tk*dy
s = s + tk*ds

//Realizar las iteraciones para encontrar el valor de x,numero de iteraciones son 8
