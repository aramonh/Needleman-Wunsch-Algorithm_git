# Needleman-Wunsch-Algorithm_git
Needleman-Wunsch-Algorithm - Adriano Ramón Hernández

Notes:

Para la paralelizacion del codigo se implementa la tecnica de loop skewing,
para evitar la dependencia de los datos.

# Sequential
Para ejecutar el algoritmo de forma secuencial.

Dentro de la carpeta "Sequential" ejecuta los siguientes comandos:

- Comandos:

1. ``` g++ nw.cpp -o app ```
2. ``` time ./app ```

# Paralell - OPENMP
Para ejecutar el algoritmo de forma paralela con openMP.

Dentro de la carpeta "openMP" ejecuta los siguientes comandos:

- Comandos:
1. ``` g++ nw.cpp -o app -fopenmp```
2. ``` time ./app ```