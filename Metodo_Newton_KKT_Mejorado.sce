//Metodo de Newton con igualdades
H =[6 -1 2 -2 0;
    -1 6 2 -2 0 ;
    2 2 6 1 1 ;
    -2 -2 1 6 2;
    0 0 1 2 7]
    
c= [-5;-5;-12;-5;-10]
//x= [0;4;-2;0;0]
//v= [0;4;-2;0;0]
A= [0 1 1 -1 1;
    0 2 1 2 1 ]
b= [2;6]
function fx=fun_obj(x)
    f1= (x(1)+x(5)-4)**4
    f2 = 0.5* (x'*H*x)
    f3 = c' *x+100
    fx= f1+f2+f3
endfunction

function gx=res_g(x)
    g1 =x(2) +x(3) - x(4) +x(5)
    g2=2 * x(2)+ x(3)+2*x(4)+x(5)
    gx=[g1;g2]
endfunction

function x_opt = cond_kkt(xv,f,A,b,max_iter,epsilon)
    [m,n] = size(A)
    for i=1:max_iter
        disp(i,"i")
        disp(xv,"xv")
        x = xv(1:n)
        v = xv(n + 1:n + m)
        [gr,H] = numderivative(f,x,[],[],"blockmat")
        gr = gr'
        f1 = gr +A'*v
        f2 = A*x - b 
        F1 = [f1;f2]
        disp(F1,"F(xv)")
        disp(norm(F1),"norm(F)")
        F2 = [H   A';A  zeros(m,m)]
        disp(F2,"H")
        if norm(f1) < epsilon then 
            break
        else
            
            d = -F2\F1
            disp(d,"d")
            x = x + d(1:n)
            v = v + d(n+1:n + m)
       end
    end
    x_opt = x
endfunction

xv = [0;4;-2;0;0;0;0]
xsol = cond_kkt(xv,fun_obj,A,b,1,0.001)
