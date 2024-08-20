# matriz de distâncias 
import numpy as np 
import matplotlib.pyplot as plt

somente_carbono_alfa = True
proteina = "2lzm.pdb"
coords = {}

def distancia_euclidiana(i, j):
	''' Retorna a distância euclidiana entre 2 átomos '''
	return ((i[0]-j[0])**2+(i[1]-j[1])**2+(i[2]-j[2])**2)**(1/2)

cont = 1
with open(proteina) as f:
	linhas = f.readlines()

	for linha in linhas:
		if linha[0:4] == "ATOM":
			x = float(linha[30:38])
			y = float(linha[38:46])
			z = float(linha[47:54])

			atom = int(linha[6:11])
			nome = linha[12:16].strip()
			res = int(linha[22:26])

			if somente_carbono_alfa:
				if nome == "CA":
					coords[cont] = (x,y,z)
					cont+=1
			else:
				coords[atom] = (x,y,z)


total = len(coords) # obtém tamanho
dists = np.zeros((total, total)) # cria uma matriz com zeros

# calcula a distância entre todos os átomos
for i in coords:
	for j in coords:
		dist = distancia_euclidiana(coords[i], coords[j])
		dists[i-1, j-1] = round(dist,2)

print(dists)
# Plota a matriz
plt.imshow(dists, cmap='binary_r') 
# cmap = binary binary_r viridis bwr tab10 inferno Spectral,vmax=30
plt.colorbar()  # Adicionar barra de cores 
plt.show()
