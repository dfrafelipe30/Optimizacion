function [tmi] = min1abFB(f,x,d,a,b,h)
    minimo = f(x + (a*d));
    tmi = a;
    for t = a+h:h:b
        temp = f(x+(t*d));
        if(temp<minimo) then
            minimo=temp;
            tmi = t;
        end
    end
endfunction
function [tmin] = Reales(f,x,d,h)
    a = -1;
    b = 1;
    error = 0.0001;
    mov = 20;
    it = 1/h;
    count = 0;
    while count <= it
        tmin = min1abFB(f,x,d,a,b,h);
        if (tmin > a & tmin < a + error) then 
            b = a;
            a = a*mov;
        elseif (tmin < b & tmin > b - error) then
            a = b;
            b = b * mov;
        else
            break;
        end
        count = count + 1;
    end
endfunction

function y = MCC(f,x,y,k,j)
    error = 0.01;
    h = 0.0001;
    I = eye(size(x)(1),size(x)(1));
    lanmda = Reales(f,x,I(:,j),h);
    y(:,j+1) = y(:,j)+lanmda*I(:,j);
    if (j < size(x)(1)) then
        j = j + 1;
        disp("Entro a la recu")
        MCC(f,x,y,k,j);
    elseif (j == size(x)(1)) then
        x(:,k+1) = y(:,size(x)(1)+1);
        if (norm(x(:,k+1) - x(:,k)) < error) then
            //minimo = y(:,j+1);
            return 
        else
            y(:,1) = x(:,k+1);
            j = 1;
            k = k + 1;
            MCC(f,x,y,k,j);
        end
    end
endfunction

function y = f(x)
    //y = (1 - x(1))^2 + 5*(x(2) - x(1)^2)^2;
    y = -log(-5 + x(1) + 3*x(2)) - log(-9 + 5 *x(1) - x(2)) - log(15 - 3 *x(1)  - x(2))
endfunction

[minimoF] = MCC(f,[-10;10],[-10;10],1,1);

