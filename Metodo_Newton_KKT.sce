function y = f(x)
    y = (x(1) - 1)^4 + (x(1) + 2 *x(2) - 3)^2  + (x(1) + 3*x(2) - x(3) + x(4) - 2)^2 + 10
endfunction


A = [1 2 3 4 ]
b = [10 ]
xv  = zeros(5,1)
[m,n] = size(A)
//primera iteracion
x = xv(1:n)

v = xv(n +1:n + m)

[gr,H] = numderivative(f,x,[],[],'blockmat')
gr = gr'

fxv = [gr + A' *v ; A*x - b]

Df = [H A'; A  zeros(m,m)]

d = -Df\fxv

xv = xv + d

for i=1:16
    
    if(norm(fxv) < 0.001)
        break
    else
        x = xv(1:n)
        v = xv(n +1:n + m)
    
        [gr,H] = numderivative(f,x,[],[],'blockmat')
        gr = gr'

        fxv = [gr + A' *v ; A*x - b]

        Df = [H A'; A  zeros(m,m)]

        d = -Df\fxv

        xv = xv + d
        disp (xv)
     end 
end

//Buscanos un punto que cumpla kkt
