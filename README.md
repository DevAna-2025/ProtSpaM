# ProtSpaM (Versión Secuencial en C++)

Esta es una adaptación **secuencial** de [ProtSpaM](https://github.com/jschellh/ProtSpaM), una herramienta para estimar distancias filogenéticas entre proteínas basada en **spaced-word matches**.  

La versión original incluye paralelismo con **OpenMP**.  
Esta versión elimina toda la paralelización para ejecutarse de manera **estrictamente secuencial**, manteniendo la misma funcionalidad principal.  

---

## ✨ Cambios respecto a la versión original
- Se creó un archivo `main_sequential.cpp` en lugar del `main.cpp` original.  
- Se eliminaron las dependencias de **OpenMP** (`#include <omp.h>`, `omp_get_wtime`, `omp_set_num_threads`, `#pragma omp parallel for`).  
- El sistema de medición de tiempo ahora utiliza `std::chrono` en lugar de funciones de OpenMP.  
- El parámetro `-t` (número de hilos) se borra ya que no es necesari para un codigo secuencial.  

---

## 🔧 Compilación
- Primero necesitas tener la carpeta /data con todo los archivos FASTA
- Revisa que lo que esta dentro del archivo filelist y que esos archivos existan en /data
- Corre lo que esta en el archivo Makefile

```bash
make
```

### 📂 Uso
El programa acepta los mismos parámetros que la versión original.  

Ejemplo de uso con múltiples archivos de entrada:

```bash
./bin/Debug/protspam -w 6 -d 40 -m 5 -l filelist -p patterns.txt
```

### ⚙️ Opciones principales
- `-w <int>` : Peso del patrón (**default:** 6)  
- `-d <int>` : Número de posiciones "don't-care" (**default:** 40)  
- `-s <int>` : Puntaje mínimo para considerar un spaced-word match como homólogo (**default:** 0)  
- `-m <int>` : Número de patrones (**default:** 5)  
- `-o <file>` : Nombre del archivo de salida con la matriz de distancias (**default:** DMat)  
- `-l <file>` : Lista de archivos de entrada en formato multifasta  
- `-p <file>` : Cargar conjunto de patrones predefinidos  
- `-z` : Guardar patrones generados en `patterns.txt`  
- `-r` : Generar puntajes individuales para cada par de secuencias (spamogramas)  

---

### 📊 Salida
El programa genera los siguientes resultados:  
- **Matriz de distancias** → archivo de salida (`DMat` por defecto).  
- **Spamogramas (scores)** → en el directorio `scores/` si se usa la opción `-r`.  
