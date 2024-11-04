import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from bacteria import bacteria
from chemiotaxis import chemiotaxis

# Parámetros iniciales
ruta_archivo = r"C:\Users\RICKY\OneDrive\Escritorio\BFOA-main\multiFasta.fasta"
numero_de_bacterias = 5
num_bacterias_aleatorias = 1
iteraciones = 30
tumbo = 1
d_atractivo = 0.1
w_atractivo = 0.2
h_repulsivo = d_atractivo
w_repulsivo = 10

# Función para clonar la mejor bacteria
def clonar_mejor(very_best, best):
    very_best.matrix.seqs = np.array(best.matrix.seqs)
    very_best.blosumScore = best.blosumScore
    very_best.fitness = best.fitness
    very_best.interaction = best.interaction

# Función para validar las secuencias
def validar_secuencias(ruta, very_best):
    temp_bacteria = bacteria(ruta)
    temp_bacteria.matrix.seqs = np.array(very_best.matrix.seqs)
    
    for i in range(len(temp_bacteria.matrix.seqs)):
        temp_bacteria.matrix.seqs[i] = temp_bacteria.matrix.seqs[i].replace("-", "")
    
    original = bacteria(ruta)
    for i in range(len(temp_bacteria.matrix.seqs)):
        if temp_bacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("********* Inconsistencia en secuencias *********")
            return False
    return True

# Función para ejecutar la versión mejorada del algoritmo
def ejecutar_algoritmo_mejorado():
    poblacion = [bacteria(ruta_archivo) for _ in range(numero_de_bacterias)]
    chemio = chemiotaxis()
    very_best = bacteria(ruta_archivo)
    global_nfe = 0
    fitness_history = []

    print("Iniciando ejecución de la versión mejorada")

    for iteracion in range(iteraciones):
        print(f"  Iteración {iteracion + 1}/{iteraciones} - Mejorada")
        
        for bacterium in poblacion:
            bacterium.tumboNado(tumbo)
            bacterium.autoEvalua()
        
        chemio.doChemioTaxis(poblacion, d_atractivo, w_atractivo, h_repulsivo, w_repulsivo)
        global_nfe += chemio.parcialNFE
        
        best = max(poblacion, key=lambda x: x.fitness)
        if best.fitness > very_best.fitness:
            clonar_mejor(very_best, best)

        fitness_history.append(very_best.fitness)

        if not validar_secuencias(ruta_archivo, very_best):
            very_best.fitness *= 0.95  # Penalización

        chemio.eliminarClonar(ruta_archivo, poblacion)
        chemio.insertRamdomBacterias(ruta_archivo, num_bacterias_aleatorias, poblacion)
    
    print(f"  Finalizó versión mejorada - Mejor Fitness: {very_best.fitness}, NFE: {global_nfe}")
    return very_best.fitness, global_nfe, fitness_history

# Almacenamiento de resultados de las corridas
resultados_mejorado = []

# Realizar 30 corridas de la versión mejorada
for i in range(30):
    print(f"\nCorrida {i + 1}/30")
    fitness_mej, nfe_mej, fitness_history = ejecutar_algoritmo_mejorado()
    resultados_mejorado.append((fitness_mej, nfe_mej, fitness_history))

# Guardado de resultados en DataFrame para análisis
df_resultados_mejorado = pd.DataFrame({
    'Corrida': range(1, 31),
    'Fitness Mejorado': [f for f, n, _ in resultados_mejorado],
    'NFE Mejorado': [n for f, n, _ in resultados_mejorado],
})

# Guardar los resultados en un archivo CSV
print("\nGuardando resultados en CSV...")
df_resultados_mejorado.to_csv("resultados_mejorados.csv", index=False)
print("Archivo CSV guardado exitosamente.")

# Generación de gráficas para análisis de mejora en Fitness y NFE
print("Generando gráficas...")

# Gráfica comparativa de Fitness
plt.figure(figsize=(12, 6))
plt.boxplot(df_resultados_mejorado['Fitness Mejorado'], labels=['Mejorado'])
plt.title("Comparación de Fitness (Versión Mejorada)")
plt.ylabel("Fitness")
plt.grid(axis='y')
plt.savefig("comparacion_fitness_mejorado.png")
plt.show()

# Gráfica de líneas comparativa para NFE
plt.figure(figsize=(12, 6))
plt.plot(df_resultados_mejorado['Corrida'], df_resultados_mejorado['NFE Mejorado'], label='NFE Mejorado', marker='o', color='orange')
plt.title("Comparación de NFE (Versión Mejorada)")
plt.xlabel("Corrida")
plt.ylabel("NFE")
plt.legend()
plt.grid()
plt.savefig("comparacion_nfe_mejorado.png")
plt.show()

# Generación de una tabla visual con resultados
print("\nResultados de las 30 corridas de la versión mejorada:")
print(df_resultados_mejorado)

# Resumen de resultados
mejor_fitness_mej = df_resultados_mejorado['Fitness Mejorado'].max()
mejor_nfe_mej = df_resultados_mejorado['NFE Mejorado'].min()

print("\nResumen de Resultados (Versión Mejorada):")
print(f"Mejor Fitness Mejorado: {mejor_fitness_mej}")
print(f"Mejor NFE Mejorado: {mejor_nfe_mej}")
