###########################################################

CÓDIGO DE LA TAREA 2 DE MÉTODOS NUMÉRICOS

@Author Daniel Vallejo Aldana (danielvallejo237 on Github)

##########################################################

Especificaciones técnicas:

El programa fue corrido en una máquina con las siguietes características:

OS: Ubuntu 22.04 LTS
Gnome Version: 42.2
g++ version: 11.2.0

Arvhivos generados:

biseccion.cpp -> Archivo con extensión *.cpp que contiene el código de bisección implementado en c++

newton.cpp -> Archivo con extensión *.cpp que contiene el código del método de Newton implementado en c++

Forma de compilación.

Ambos códigos se compilan de la misma forma incluyendo la librería fparser que debe de encontrarse en la misma
carpeta en la que se está ejecutando el código

La forma de compilación de los códigos es la siguiente

g++ (file).cpp  fparser/fparser.cpp

donde file es newton o biseccion

FORMA DE EJECUCIÓN:

La forma de ejecución de los códigos cambia de acuerdo al tipo de código que estemos ejecutando

código de bisección

./a.out 'función' a b

Código de Newton

./a.out 'funcion' 'derivada' xinit

El programa debe de imporimir en consola el valor de la raiz, el error de la condición de paro y el número de iteraciones.


NOTA:

Ambos programas generan un archivo de nombre puntos.txt, dicho archivo contiene las soluciones candidatas generadas a lo
largo de las iteraciones y que son utilizadas posteriormente para graficar en python el error relativo a lo largo de las
iteraciones. El archivo de graficación no se incluye en el archivo comprimido de esta tarea.
