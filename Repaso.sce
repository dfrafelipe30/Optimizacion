//function y = f(x)
    //y = exp((x')* A * x)
//endfunction

A = [1 2 ; 2  5];

function y = f(x)
    y = 5 *x(1)^4 + (x(1) + x(2) - 1)^2 + 100 
endfunction

function tp = Opti(t1,t2,t3,x,d,f)
    //La funcion de tao2 y asi sucesivamente vamos a calcularlo es f(x  + t * d)
    div = ((t2 - t1) * (f(x + t2 * d) - f(x + t3 * d)) - (t2 - t3) * (f(x + t2 * d) - f(x + t1 * d)) )
    sup = ((t2 - t1)^2 * (f(x + t2 * d) - f(x + t3 * d)) - (t2 - t3)^2 * (f(x + t2 * d) - f(x + t1 * d)) )
    tp = t2 - 0.5 * (sup/div)
endfunction
x = [3 ; 4];

//d =  numderivative(f,x)

d = [-5 ; -1]


//Tenemos que buscar un |t| > ecilon
t1 = 0.1
aux1 = x + t1 *d   // Punto que cumple que f(x + t *d) < f(x) 
y1 = f(aux1)

t2 = 0.3
aux2 = x + t2 *d
y2 = f(aux2)


t3 = 0.7
aux3 = x + t3 * d
y3 = f(aux3)

t4 = 1.5
aux4 = x + t4 * d
y4 = f(aux4)

//Calculando el nuevo tp para la nueva iteracion

tp = Opti(t2,t3,t4,x,d,f)

yk = f(x + tp * d)   //Valor en tp


//Vamos a volver a hacer la iteracion con t2  ya que t2 es mejor que tp (es el tao de la funcion )


t5 = 0.7
aux5 = x + t5 *d   // Punto que cumple que f(x + t *d) < f(x) 
y5 = f(aux5)


t6 = 0.725        //Este punto es mejor que t2.
aux6 = x + t6 *d
y6 = f(aux6)


t7 = 0.75
aux7 = x + t7 * d
y7 = f(aux7)

//Vamos a calular el nuevo valor de tp 
t8 = Opti(t5,t6,t7,x,d,f)
yk = f(x + t8 * d)

//Para la siguiente iteracion vamos a calcular con tp ya que f(tp) < f(t5)
t9 = 0.7
aux9 = x + t9 *d    
y9 = f(aux9)


t10 = 0.725        
aux10 = x + t10 *d
y10 = f(aux10)


t11 = 0.75
aux11 = x + t11 * d
y11 = f(aux11)

