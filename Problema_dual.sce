function y = f(x)
    y = -x(1) - 1.4*x(2)
endfunction

function G = g(x)
    gd =[-1 0 0 0 0;0 -1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0; 0 0 0 0 -1]
    q = [-0.1; -0.2; 0 ; 0 ; 0]
    G =  gd * x - q
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
function escV(x,tit)
    n = length(x)
    printf('%s :',tit)
    for i = 1:n
        printf('%10.4f',x(i))
    end
    printf('\n')
endfunction

x =[0.2 ; 1.45; 2.35; 2.7;2.8]

G =[-1 0 0 0 0;0 -1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0; 0 0 0 0 -1]

q = [-0.1; -0.2; 0 ; 0 ; 0]

u = [1;1;1;1;1]

v = [0;0;0] 

c  = [-1 ; -1.4 ; 0 ; 0 ;0]

A = [1  1 1 0 0 ;1 2 0 1 0;1 0 0 0 1]

b = [4; 5.8 ; 3]

max_iter = 8
m = 5
function x_opt = M_P_D(x,u,v,miu,epsilon)
    for i = 1:max_iter
        disp(x,"x")
        gx = g(x)
        max_gx = max(gx)
        etam = -u' * gx
        disp(u,"u")
        disp(v,"v")
        disp(gx,"gx")
        disp(max_gx,"max(g(x))")
        disp(etam,"etam")
        tao = (miu * m) / etam
        disp(tao,"tao")
        e = ones(m,1)
        rd = c + G'*u + A'*v
        rc = -u .* (gx) - e/tao
        rp = A*x - b
        r = [rd;rc;rp]
        escV(r,"r_tao")
        norma_r = norm(r)
        escV(norma_r,"||r_tao||")
        if norm(r) < epsilon then
            break
        else
            primer = zeros(5,5)
            segundo = zeros(3,5)
            tercero = zeros(3,3)
            cuarto = zeros(5,3)
            diagU = -diag(u)
            diagG = -diag(gx)
            matriz = [primer G' A' ;diagU.*G diagG cuarto;A segundo tercero]
    
            //direccion solucionando el gradiente
            disp(matriz,"M:")
            sol = -matriz\r
            dx = sol(1:5)
            du = sol(6:10)
            dv = sol(11:13)
            disp(dx,"dx")
            disp(du,"du")
            disp(dv,"dv")
            tmaxU = tmaxu(u,du)
            tmaxk = tmaxX(G,q,x,dx)
            tmax = min(tmaxU,tmaxk)
            disp(tmax,"tmax")
            tk = 0.99 * tmax
            disp(tk,"tk")
            x = x + tk*dx
            u = u + tk*du
            v = v + tk*dv
        end
        x_opt = x
    end
endfunction

x_optimo = M_P_D(x,u,v,10,0.0001)
disp(x_optimo,"x_optimo")
