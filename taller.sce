function y = f(x)
    y = abs(x(1) + x(2))
endfunction

function [tmin,minimo] = min1abFB(f,x,d,a,b,h)
    minimo = f(x + (a*d));
    tmin = a;
    for t = a+h:h:b
        temp = f(x+(t*d));
        if(temp<minimo) then
            disp('Cambio temp')
            minimo=temp;
            tmin = t;
        end
    end
endfunction

a = -2;
b = 2;

x = [-2;2];
d = [-1;0];


[t,val] = min1abFB(f,x,d,a,b,0.0001);

