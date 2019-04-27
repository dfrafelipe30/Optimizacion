function y = f(x)
    y = (x(1) - 1)^4 + (x(1) + 2 *x(2) - 3)^2  + (x(1) + 3*x(2) - x(3) + x(4) - 2)^2 + 10
endfunction

A = [1 2 3 4 ]
b = [10 ]

x = [10 ; 0 ; 0 ; 0]

B = [-2 -3 -4; 1 0 0;0 1 0;0 0 1]

z = [0 ; 0 ; 0]


xsol = x + B * z

[gr,H] = numderivative(f,xsol,[],[],'blockmat')

gr = gr'

Hfsol = B' * H * B
grfsol = B'* gr

d = -Hfsol\grfsol

//segunda iteracion
z1 = z + d

xsol1 = x + B * z1

[gr1,H1] = numderivative(f,xsol1,[],[],'blockmat')

gr1 = gr1'

Hfsol1 = B' * H1 * B
grfsol1 = B'* gr1

d1 = -Hfsol1\grfsol1

//tercera iteracion
z2 = z1 + d1

xsol2 = x + B * z2

[gr2,H2] = numderivative(f,xsol2,[],[],'blockmat')

gr2 = gr2'

Hfsol2 = B' * H2 * B
grfsol2 = B'* gr2

d2 = -Hfsol2\grfsol2
//cuarta iteracion
z3 = z2 + d2

xsol3 = x + B * z3

[gr3,H3] = numderivative(f,xsol3,[],[],'blockmat')

gr3 = gr3'

Hfsol3 = B' * H3 * B
grfsol3 = B'* gr3

d3 = -Hfsol3\grfsol3
for i = 1:10
    z3 = z2 + d2

    xsol3 = x + B * z3

    [gr3,H3] = numderivative(f,xsol3,[],[],'blockmat')

    gr3 = gr3'

    Hfsol3 = B' * H3 * B
    grfsol3 = B'* gr3
    z2 = z3
    d3 = -Hfsol3\grfsol3
    d2 = d3 
    disp(xsol,Hfsol3,grfsol3)
end
