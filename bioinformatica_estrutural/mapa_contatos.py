# matriz de distâncias 
import numpy as np 
import matplotlib.pyplot as plt

somente_carbono_alfa = False
proteina = "2lzm.pdb"
coords = {}
cutoff = 3.9

def distancia_euclidiana(i, j):
	''' Retorna a distância euclidiana entre 2 átomos '''
	return ((i[0]-j[0])**2+(i[1]-j[1])**2+(i[2]-j[2])**2)**(1/2)

cont = 0 # número dos resíduos
with open(proteina) as f:
	linhas = f.readlines()

	for linha in linhas:
		if linha[0:4] == "ATOM":
			x = float(linha[30:38])
			y = float(linha[38:46])
			z = float(linha[47:54])

			atom = int(linha[6:11])
			atom_nome = linha[12:16].strip()
			res = int(linha[22:26]) 
			res_nome = linha[16:20].strip()			

			if atom_nome == "N": # 1º átomo é N
				cont+=1 # corrige o problema de numeração não contínua
			
			coords[atom] = (x,y,z,res_nome,atom_nome,cont)

contatos = np.zeros((cont, cont)) # cria uma matriz com zeros

aceptores = [
	"ALA:O","ARG:O","ASN:O","ASN:OD1","ASP:O","ASP:OD1",
	"ASP:OD2","CYS:O","CYS:SG","GLN:O","GLN:OE1","GLU:O",
	"GLU:OE1","GLU:OE2","GLY:O","HIS:O","ILE:O","LEU:O",
	"LYS:O","MET:O","MET:SD","PHE:O","PRO:O","SER:O",
	"SER:OG","THR:O","THR:OG1","TRP:O","TYR:OH","VAL:O"
]

doadores = [
	"ALA:N","ARG:N","ARG:NE","ARG:NH1","ARG:NH2","ASN:N",
	"ASN:ND2","ASP:N","CYS:N","GLN:N","GLN:NE2","GLU:N",
	"GLY:N","HIS:N","HIS:ND1","HIS:NE2","ILE:N","LEU:N",
	"LYS:N","LYS:NZ","MET:N","PHE:N","SER:N","THR:N",
	"TRP:N","TRP:NE1","TYR:N","VAL:N"
]


# calcula a distância entre todos os átomos
for i in coords:
	for j in coords:
		dist = distancia_euclidiana(coords[i], coords[j])
		atomo1 = coords[i][3]+":"+coords[i][4]
		atomo2 = coords[j][3]+":"+coords[j][4]
		res1 = coords[i][5]
		res2 = coords[j][5]

		if dist < cutoff and (
			(atomo1 in aceptores and atomo2 in doadores) or 
			(atomo2 in aceptores and atomo1 in doadores)
		): 
			# print("Contato:", res1, res2, atomo1, atomo2)
			contatos[res1-1, res2-1] = 1

#print(contatos)
# Plota a matriz
plt.imshow(contatos, cmap='binary') 
# cmap = binary binary_r viridis bwr tab10 inferno Spectral,vmax=30
plt.show()