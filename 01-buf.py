import os
import subprocess
    
# Crear la carpeta "archivos_temporales" si no existe
output_folder = "03-proteinortho-graph"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
buff = []
i = 1
block_started = False  # Variable para indicar si se ha comenzado un bloque
    
with open('02-proteinortho/myproject.proteinortho-graph', "r") as f:
    for line in f:
        # Verificar si la línea comienza con '#' seguido de un espacio y un número
        if line.startswith("# ") and line.split()[1][0].isdigit():
            continue  # Saltar esta línea y pasar a la siguiente iteración del bucle
            
        if line.startswith("#") and "CARD.fa" in line:
            block_started = True  # Indicar que se ha comenzado un nuevo bloque
            # Si hay contenido en buff (es decir, si no es el primer bloque), guardarlo en un archivo
            if buff:
                with open(os.path.join(output_folder, f'file{i}.txt'), 'w') as output:
                    # Escribir todas las líneas en buff excepto la última (que contiene '#')
                    output.write(''.join(buff[:-1]))
                i += 1
                buff = []  # Reiniciar el buffer
            buff.append(line)  # Agregar la línea que inicia el bloque
            
        elif block_started:
            # Si ya se ha comenzado un bloque, añadir la línea actual al buffer
            buff.append(line)
            # Si la línea actual comienza con '#' y no es la primera línea del bloque,
            # indicar que el bloque ha terminado y guardar el contenido en un archivo
            if line.startswith("#") and line.strip() != buff[0].strip():
                block_started = False
                with open(os.path.join(output_folder, f'file{i}.txt'), 'w') as output:
                    # Escribir todas las líneas en buff excepto la última (que contiene '#')
                    output.write(''.join(buff[:-1]))
                i += 1
                buff = []  # Reiniciar el buffer
    
# Guardar la sección final si hay contenido en buff
if buff:
    with open(os.path.join(output_folder, f'file{i}.txt'), 'w') as output:
        # Escribir todas las líneas en buff excepto la última (que contiene '#')
        output.write(''.join(buff[:-1]))
   
# Obtener la lista de archivos en la carpeta de salida
file_list = [filename for filename in os.listdir(output_folder)]
    
for filename in file_list:
    with open(os.path.join(output_folder, filename), "r") as f:
        # Leer todas las líneas del archivo
        lines = f.readlines()
        
    with open(os.path.join(output_folder, filename), "w") as f:
        # Escribir todas las líneas excepto la primera con el espacio y el numeral eliminados
        f.writelines(lines[0].lstrip("# "))
        f.writelines(lines[1:])


# Directorio donde se encuentran los archivos
directory = "03-proteinortho-graph"

# Iterar sobre cada archivo que comienza con "file"
for file in os.listdir(directory):
    if file.startswith("file") and file.endswith(".txt"):
        file_path = os.path.join(directory, file)
        output_file = f"{file_path[:-4]}_cut.txt"  # Reemplaza .txt con _cut.txt
        
        # Ejecutar el comando cut para extraer la primera y segunda columna
        with open(output_file, 'w') as out_f:
            subprocess.run(["cut", "-f1,2", file_path], stdout=out_f)

print("Procesamiento completado.")

try:
    subprocess.run(["Rscript", "heatmap-ortos2-SD.R"], check=True)
    print("Script de R ejecutado con éxito.")
except subprocess.CalledProcessError as e:
    print(f"Error al ejecutar el script de R: {e}")
