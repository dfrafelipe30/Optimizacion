function y = f(x)
    y = x(1) + 1.1*x(2)
endfunction

//Desigualdades
function gx = g(x)
    y1 = x(1)**2 + x(2)**2 - 2
    y2 = (x(1) - 2)**2 + x(2)**2 - 2
    gx = [y1;y2]
endfunction
//Gradientes
function gr = gr1(x)
    g1 = [2*x(1);2*x(2)]
    g2 = [2*(x(1)-2);2*x(2)]
    gr = [g1 g2] 
endfunction
//Hesianas
H1 = [2 0;0 2]
H2 = [2 0;0 2]

function [gr,HB] = B(x,H)
    disp(H,'H')
    gra = gr1(x)
    disp(gra,'gra')
    gx = g(x)
    n = length(x)
    gr = 0
    HB = zeros(n,n)
    /*
    gr = (grad1/g1) + (grad2/g2)
    gr = -gr
    parte1 = ((grad1*grad1')/g1**2)  + ((grad2*grad2')/g2**2)
    parte2 = (H1/g1) + (H2/g2)
    parte2 = -parte2
    H = parte1 + parte2*/
    for i = 1:size(gr,2)
        gri = gra(:,i)
        gr = gr + gri/gx(i)
        HB = HB + (gri*gri')/(gx(i)**2) 
        HB = HB - H/gx(i)
    end
    gr = -gr
 
        
endfunction

function escV(x,tit)
    n = length(x)
    printf('%s :',tit)
    for i = 1:n
        printf('%10.4f',x(i))
    end
    printf('\n')
endfunction

function tmx = tmax_g(g, x, d, M)
    //
    // tmx = max { t  :  g( x + t d ) <= 0,
    //             0 < t <= M }
    //
    // se supone que  g(x) < 0
    //
    //                0 < M
    //
    // utiliza biseccion (o dicotomia) y secante
    //
    //*************************
    EPS0 = 1.0e-8
    EPS  = 1.0e-6
    //*************************
    
    //escrVectED(x, 'x en tmax'), escrVectED(d, 'd')
    g0 = g(x)
    if max(g0) >= -EPS0
        printf('El punto no es estrictamente factible.\n')
        tmx = -1000
        return
    end
    gM = g(x + M*d)
    //escrVectED(gM, 'g(x+Md)')
    if max(gM) <= EPS0
        tmx = M
        return
    end
    
    a = 0
    b = M
    c0 = M/3
    c1 = 2*M/3
    gc0 = max( g( x + c0*d ) )
    gc1 = max( g( x + c1*d ) )
    
    while (b-a) >= EPS
    //escrVectED([a b], 'a b')
    deno = gc1 - gc0
    if abs(deno) < EPS0
        c2 = c1 + 0.01
    else
        c2 = c1 - gc1*(c1-c0)/deno
    end

    gc2 = max( g( x + c2*d ) )
    if 0 < c2 & c2 <= M & abs(gc2) <= EPS0
        //printf('secante\n')
        tmx = c2
        return
    end
    c0 = c1, gc0 = gc1
    c1 = c2, gc1 = gc2
    
    m = (a+b)/2
    gm = max( g(x + m*d) )
    
    if abs(gm ) <= EPS0
        tmx = m
        return
    end
    
    if gm > 0
        b = m
    else
        a = m
    end
    end
    tmx = (a+b)/2
endfunction // tmax_g

//printf('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n=========================\n')
x = [0.9;0]
//x = [ 0.8905038;-0.4159503]
tao = 10

for k = 1:10
   printf('\n\nK = %d\n',k)
   escV(x,'X')
   gx = g(x)
   escV(gx,'gx')
   fp = numderivative(f,x)
   fp = fp'
   [gradiente,Hessiana] = B(x,H1)
   fpt = tao*fp +  gradiente
   fdpt = Hessiana
   if(norm(fpt) < 0.001) then 
       break
   else
       d = -fdpt\fpt
       //x = x + d
       xm = x + d
       escV(xm,"xm")
       gxm  = g(xm)
       escV(gxm,'gxm')
       tm = tmax_g(g,x,d,1)
       //disp(tm,"el valor de tm")                                    
       if max(gxm) < 0 then 
           x = xm
          ´disp("El t es igual a 1")
       else
           x = x + 0.99*tm * d     //(0.356308)*d
           disp(tm,"El valor de tm")
       end  
       //disp(x)
       //disp(d)
       escV(x,'x')
       escV(d,'d')
   end
end


