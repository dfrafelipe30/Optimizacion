function [tmin] = Reales(f,x,d,h)
    a = -1;
    b = 1;
    error = 0.0001;
    mov = 20;
    it = 1/h;
    count = 0;
    while count <= it
        tmin = minlabFB(f,x,d,a,b,h);
        if  tmin > a & tmin < a + error
            b = a;
            a = a*mov;
        elseif tmin < b & tmin > b - error
            a = b;
            b = b * mov;
        else
            break;
        end
    end
endfunction
