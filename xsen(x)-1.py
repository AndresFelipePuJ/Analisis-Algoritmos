import numpy as np
#funcion define la funcion dada
def f(x):
    return x * np.sin(x) - 1
#funcion posfalsa
def posfalsa(x0,x1,e):
    paso = 1
    condicion = True
    while condicion:
        #formula ecuacion de la recta
        x2 = x0 - (x1-x0) * f(x0)/( f(x1) - f(x0) )
        print('Iteracion0-%d, x2 = %0.6f and f(x2) = %0.6f' % (paso, x2, f(x2)))
#valida si son positivos o negativos y se iguala a x2 si es el caso
        if f(x0) * f(x2) < 0:
            x1 = x2
            #si es positivo x0=x2
        else:
            x0 = x2

        paso = paso + 1
        condicion = abs(f(x2)) > e

    print('\nRaiz requerida es: %0.8f' % x2)



x0 = input('Primer intervalo ')
x1 = input('Segundo Intervalo ')
e = input('Error torelable: ')

# casteo a float
x0 = float(x0)
x1 = float(x1)
e = float(e)




if f(x0) * f(x1) > 0.0:
    print('Error.')
    print('Intente denuevo con otros valores.')
else:
    posfalsa(x0,x1,e)