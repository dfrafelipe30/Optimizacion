/*function x = M_OL(c,A,b,[x0,y0,s0]ep,ed,ec,maxit)
    [x,y,s] = [x0,y0,s0]
    for i=1:maxit
        rp = A*x - b
        rd = A'*y + s -c 
        rc = x.*s
        r = [rd;rp;rc]
        rop = (norm(rp))/(1 + norm(b))
        rod = (norm(rd))/(1 + norm(c))
        roc = (norm(rc))/(1 + norm(x))
        if rop < ep & rod < ed & roc < ec then
            break
        else
            matriz = [zeros(size(A)(2),size(A)(2))]
        end
endfunction*/
function y = f(x)
    //y = x(1) + 1.4*x(2)
    y = -x(1) - 1.4*x(2) //Problema dual
endfunction
function tmx = tmaxu(u,du)
    tmx = 1e10
    n = length(u)
    for i =1:n
        if du(i) < 0 then
            tmx = min(tmx,-u(i)/du(i))
        end
    end
endfunction


x= [1; 1; 1;1; 1]

y = [0;0;0]

s = [1;1;1;1;1]


c = [-1 ; -1.4; 0 ; 0 ; 0]

A = [1 1 1 0 0;1 2 0 1 0;1 0 0 0 1]

b = [4;5.8;3]


primero = zeros(5,5)
segundo = zeros(3,3)
tercero = zeros(3,5)
cuarto = zeros(5,3)
identidad = eye(5,5)

rd = A'*y + s - c 
rp = A*x - b
rc = x.*s

r = [rd;rp;rc]

matriz = [primero A' identidad;A segundo tercero; diag(s) cuarto diag(x)]

sol = -matriz\r


dxa = sol(1:5)
dya = sol(6:8)
dsa = sol(9:13)

tpa = min(1,tmaxu(x,dxa))
tda = min(1,tmaxu(s,dsa))

x_afin = (x + tpa*dxa)'
s_afin = s + tda*dsa

ua = (x_afin*s_afin)/5  //n: numero de variables



u = (x'*s)/5  //n:nuermo de variables


sigma = (ua/u)**3



e = ones(5,1)
rc = (-sigma*u*e) + (dxa.*dsa)

primero = zeros(1,5)
segundo = zeros(1,3)

r_corre = [primero';segundo';rc]

disp(r_corre)

sol_corr = -matriz\r_corre

dxc = sol_corr(1:5)
dyc = sol_corr(6:8)
dsc = sol_corr(9:13)

dx = dxa + dxc
dy = dya + dyc
ds = dsa + dsc

tp = min(1,0.99*tmaxu(x,dx))
td = min(1,0.99*tmaxu(s,ds))

x = x + tp * dx
y = y + td * dy
s = s + td * ds

//Segunda iteracion
rd = A'*y + s - c 
rp = A*x - b
rc = x.*s

r = [rd;rp;rc]

matriz = [primero A' identidad;A segundo tercero; diag(s) cuarto diag(x)]

sol = -matriz\r


dxa = sol(1:5)
dya = sol(6:8)
dsa = sol(9:13)

tpa = min(1,tmaxu(x,dxa))
tda = min(1,tmaxu(s,dsa))

x_afin = (x + tpa*dxa)'
s_afin = s + tda*dsa

ua = (x_afin*s_afin)/5  //n: numero de variables



u = (x'*s)/5  //n:nuermo de variables


sigma = (ua/u)**3



e = ones(5,1)
rc = (-sigma*u*e) + (dxa.*dsa)

r_corre = [0; 0 ;-rc]

sol_corr = -matriz\r_corre

dxc = sol_corr(1:5)
dyc = sol_corr(6:8)
dsc = sol_corr(8:13)

dx = dxa + dxc
dy = dya + dyc
ds = dsa + dsc

tp = min(1,0.99*tmaxu(x,dx))
td = min(1,0.99*tmaxu(s,ds))

x = x + tp * dx
y = y + td * dy
s = s + td * ds

//Tercera iteracion
rd = A'*y + s - c 
rp = A*x - b
rc = x.*s

r = [rd;rp;rc]

matriz = [primero A' identidad;A segundo tercero; diag(s) cuarto diag(x)]

sol = -matriz\r


dxa = sol(1:5)
dya = sol(6:8)
dsa = sol(9:13)

tpa = min(1,tmaxu(x,dxa))
tda = min(1,tmaxu(s,dsa))

x_afin = (x + tpa*dxa)'
s_afin = s + tda*dsa

ua = (x_afin*s_afin)/5  //n: numero de variables



u = (x'*s)/5  //n:nuermo de variables


sigma = (ua/u)**3



e = ones(5,1)
rc = (-sigma*u*e) + (dxa.*dsa)

r_corre = [0; 0 ;-rc]

sol_corr = -matriz\r_corre

dxc = sol_corr(1:5)
dyc = sol_corr(6:8)
dsc = sol_corr(8:13)

dx = dxa + dxc
dy = dya + dyc
ds = dsa + dsc

tp = min(1,0.99*tmaxu(x,dx))
td = min(1,0.99*tmaxu(s,ds))

x = x + tp * dx
y = y + td * dy
s = s + td * ds
//Cuarta iteracion
rd = A'*y + s - c 
rp = A*x - b
rc = x.*s

r = [rd;rp;rc]

matriz = [primero A' identidad;A segundo tercero; diag(s) cuarto diag(x)]

sol = -matriz\r


dxa = sol(1:5)
dya = sol(6:8)
dsa = sol(9:13)

tpa = min(1,tmaxu(x,dxa))
tda = min(1,tmaxu(s,dsa))

x_afin = (x + tpa*dxa)'
s_afin = s + tda*dsa

ua = (x_afin*s_afin)/5  //n: numero de variables



u = (x'*s)/5  //n:nuermo de variables


sigma = (ua/u)**3



e = ones(5,1)
rc = (-sigma*u*e) + (dxa.*dsa)

r_corre = [0;0;-rc]

sol_corr = -matriz\r_corre

dxc = sol_corr(1:5)
dyc = sol_corr(6:8)
dsc = sol_corr(8:13)

dx = dxa + dxc
dy = dya + dyc
ds = dsa + dsc

tp = min(1,0.99*tmaxu(x,dx))
td = min(1,0.99*tmaxu(s,ds))

x = x + tp * dx
y = y + td * dy
s = s + td * ds
//Quinta iteracion
rd = A'*y + s - c 
rp = A*x - b
rc = x.*s

r = [rd;rp;rc]

matriz = [primero A' identidad;A segundo tercero; diag(s) cuarto diag(x)]

sol = -matriz\r


dxa = sol(1:5)
dya = sol(6:8)
dsa = sol(9:13)

tpa = min(1,tmaxu(x,dxa))
tda = min(1,tmaxu(s,dsa))

x_afin = (x + tpa*dxa)'
s_afin = s + tda*dsa

ua = (x_afin*s_afin)/5  //n: numero de variables



u = (x'*s)/5  //n:nuermo de variables


sigma = (ua/u)**3



e = ones(5,1)
rc = (-sigma*u*e) + (dxa.*dsa)

r_corre = [0;0;-rc]

sol_corr = -matriz\r_corre

dxc = sol_corr(1:5)
dyc = sol_corr(6:8)
dsc = sol_corr(8:13)

dx = dxa + dxc
dy = dya + dyc
ds = dsa + dsc

tp = min(1,0.99*tmaxu(x,dx))
td = min(1,0.99*tmaxu(s,ds))

x = x + tp * dx
y = y + td * dy
s = s + td * ds
//Sexta iteracion
rd = A'*y + s - c 
rp = A*x - b
rc = x.*s

r = [rd;rp;rc]

matriz = [primero A' identidad;A segundo tercero; diag(s) cuarto diag(x)]

sol = -matriz\r


dxa = sol(1:5)
dya = sol(6:8)
dsa = sol(9:13)

tpa = min(1,tmaxu(x,dxa))
tda = min(1,tmaxu(s,dsa))

x_afin = (x + tpa*dxa)'
s_afin = s + tda*dsa

ua = (x_afin*s_afin)/5  //n: numero de variables



u = (x'*s)/5  //n:nuermo de variables


sigma = (ua/u)**3



e = ones(5,1)
rc = (-sigma*u*e) + (dxa.*dsa)

r_corre = [0;0;-rc]

sol_corr = -matriz\r_corre

dxc = sol_corr(1:5)
dyc = sol_corr(6:8)
dsc = sol_corr(8:13)

dx = dxa + dxc
dy = dya + dyc
ds = dsa + dsc

tp = min(1,0.99*tmaxu(x,dx))
td = min(1,0.99*tmaxu(s,ds))

x = x + tp * dx
y = y + td * dy
s = s + td * ds
//Septima iteracion
rd = A'*y + s - c 
rp = A*x - b
rc = x.*s

r = [rd;rp;rc]

matriz = [primero A' identidad;A segundo tercero; diag(s) cuarto diag(x)]

sol = -matriz\r


dxa = sol(1:5)
dya = sol(6:8)
dsa = sol(9:13)

tpa = min(1,tmaxu(x,dxa))
tda = min(1,tmaxu(s,dsa))

x_afin = (x + tpa*dxa)'
s_afin = s + tda*dsa

ua = (x_afin*s_afin)/5  //n: numero de variables



u = (x'*s)/5  //n:nuermo de variables


sigma = (ua/u)**3



e = ones(5,1)
rc = (-sigma*u*e) + (dxa.*dsa)

r_corre = [0;0;-rc]

sol_corr = -matriz\r_corre

dxc = sol_corr(1:5)
dyc = sol_corr(6:8)
dsc = sol_corr(8:13)

dx = dxa + dxc
dy = dya + dyc
ds = dsa + dsc

tp = min(1,0.99*tmaxu(x,dx))
td = min(1,0.99*tmaxu(s,ds))

x = x + tp * dx
y = y + td * dy
s = s + td * ds
//Octava iteracion
rd = A'*y + s - c 
rp = A*x - b
rc = x.*s

r = [rd;rp;rc]

matriz = [primero A' identidad;A segundo tercero; diag(s) cuarto diag(x)]

sol = -matriz\r


dxa = sol(1:5)
dya = sol(6:8)
dsa = sol(9:13)

tpa = min(1,tmaxu(x,dxa))
tda = min(1,tmaxu(s,dsa))

x_afin = (x + tpa*dxa)'
s_afin = s + tda*dsa

ua = (x_afin*s_afin)/5  //n: numero de variables



u = (x'*s)/5  //n:nuermo de variables


sigma = (ua/u)**3



e = ones(5,1)
rc = (-sigma*u*e) + (dxa.*dsa)

r_corre = [0;0;-rc]

sol_corr = -matriz\r_corre

dxc = sol_corr(1:5)
dyc = sol_corr(6:8)
dsc = sol_corr(8:13)

dx = dxa + dxc
dy = dya + dyc
ds = dsa + dsc

tp = min(1,0.99*tmaxu(x,dx))
td = min(1,0.99*tmaxu(s,ds))

x = x + tp * dx
y = y + td * dy
s = s + td * ds

disp(x)

