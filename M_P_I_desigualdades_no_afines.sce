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

//x = [0.9;0]
x = [ 0.8905038;-0.4159503]
tao = 10

for k = 1:2
    
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
       x = x + d
       xm = x + d
       gxm  = g(xm)
       escV(gxm,'gxm')
                                           
       if max(gxm) < 0 then 
           x = xm
       else
           x = x + (0.356308)*d
           printf("El paso de newton esta muy grande\n")
       end  
       //disp(x)
       //disp(d)
       escV(x,'x')
       escV(d,'d')
   end
end


