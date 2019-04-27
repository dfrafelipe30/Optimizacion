function y = f(x)
    y = -1 * x(1) - 1.4*x(2)
endfunction


function [di,D] = delta(x)
    di = [1/(4 -x(1) -x(2));1/(5.8 - x(1) - 2*x(2));1/(3 - x(1));1/(x(1));1/(x(2))]
    D = diag(di)
endfunction

t = 1
x = [1;1]
A = [1 1;1 2;1 0; -1 0;0 -1]
b= [4;5.8;3;0;0]
c = [-1 ; -1.4]

for k = 1:10000
    [d,D] = delta(x)
    if(norm(d) < 0.001)
        break
    else
        B0 = A'* d
        B1 = A' * D**2 * A
        f0 = t* c + B0
        f1 = B1 
        dire = -f1\f0
        disp(dire)
        x = x + dire
        disp(x)
    end
end     
