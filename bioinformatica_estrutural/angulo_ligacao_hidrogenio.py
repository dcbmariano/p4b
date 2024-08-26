import math

def dist(D, H, A):
    '''retorna a distancia euclidiana entre da, dh e ah'''
    DH = [H[i] - D[i] for i in range(3)] 
    DA = [D[i] - A[i] for i in range(3)]
    AH = [H[i] - A[i] for i in range(3)]

    da = math.sqrt(sum(DA[i] ** 2 for i in range(3)))
    dh = math.sqrt(sum(DH[i] ** 2 for i in range(3)))
    ah = math.sqrt(sum(AH[i] ** 2 for i in range(3)))

    return da, dh, ah

def lc(da, dh, ah):
    '''usa a lei dos cossenos para calcular ângulo theta'''
    cos_theta = (dh**2 + ah**2 - da**2) / (2 * dh * ah)
    theta = math.acos(cos_theta)
    theta_grau = math.degrees(theta)

    return round(theta_grau,2)

# PDB-ID: 2LZM
D = [39.395, 13.931, 7.725]  #R145 - NH1
H = [39.286, 14.234, 8.686]  #R145 - HH12
A = [38.869, 13.954, 10.575] #E11 - OE1

da, dh, ah = dist(D, H, A)
angulo = lc(da, dh, ah)

print("Detectando LIGAÇÕES DE HIDROGÊNIO")
print("A distância entre aceptor e doador é", round(da,2),"angstrom")
print("Ângulo da ligação de hidrogênio é", angulo,'graus')

if angulo > 90 and da <= 3.9:
    print("Existe uma ligação de hidrogênio")

print(da,dh,ah)
# validar em https://www.omnicalculator.com/pt/matematica/angulo-triangulos