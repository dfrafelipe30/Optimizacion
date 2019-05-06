function fx = h(x)
    fx = (x(1)-1)**4 + (x(1)+2*x(2)-3)**2 + (x(1) + 3*x(2) - x(3) + x(4) - 2)**2 + 10
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
            z = z + d
        end
    end
endfunction
/*for k=1:20
    x = x_m + B*z
    [gr,H] = numderivative(h,x,[],[],"blockmat")
    gr = gr'
    Hz = B'*H*B
    grz = B'*gr
    d = -Hz\grz
    z = z + d
end*/
//Método de Newton quitando igualdades
