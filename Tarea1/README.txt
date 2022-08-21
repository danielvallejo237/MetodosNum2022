###########################################################
README DE LA TAREA 1
@Author Daniel Vallejo Aldana (danielvallejo237 on Github)
###########################################################

Especificaciones técnicas:

El programa fue corrido en una máquina con las siguietes características:

OS: Ubuntu 22.04 LTS
Gnome Version: 42.2
g++ version: 11.2.0
GTK: 3.0

Archivos:

El código se encuentra en el programa graficador.cpp, dicho código es el que habrá que compilar

Comando de compilación:

El comando para compilar el código es el siguiente

g++ graficador.cpp $(pkg-config --libs gtk+-3.0 --cflags cairo) fparser/fparser.cc
 
Comando de ejecución:

El comando de ejecución del código es

./a.out f a b

Donde f es el string que contiene la función a evaluar y a y b son valores flotantes para delimitar el rango
de la función

NOTA:
El folder que contiene fparser debe de incluirse en la misma carpeta donde se encuentra el código principal


