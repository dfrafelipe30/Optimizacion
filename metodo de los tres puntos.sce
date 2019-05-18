function tp = interpCuadr(t1,t2,t3,f1,f2,f3)
    //
    // interpolacion cuadratica
    //
    t21 = t2 - t1
    t23 = t2 - t3
    c1 = t21*(f2 - f3)
    c2 = t23*(f2 - f1)
    deno = c1 - c2
    //disp(deno, 'deno')
    if abs(deno) < 1e-8
        printf('Denominador nulo.\n')
        tp = 0
        return
    end
    tp = t2 - 0.5*(t21*c1 - t23*c2)/deno        
endfunction

function t = metodos_tres_puntos(x,direc,t1,h_inicio,max_Iter)
    h = h_inicio
    for(i = 0:max_Iter)
        if i == max_Iter then
            disp("No se encontro minimizador")
            break
        end
        
    end
endfunction
