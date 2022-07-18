"""Programa para resolução do exercício 5.1 do livro."""

# Importações
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import materials as mat
import os
from sympy import *

class plate:
    """Cria a classe da placa."""

    def __init__(self, a, b, t, material):
        """Método de inicialização da placa."""
        self.a = a
        self.b = b
        self.t = t
        self.D = (material.E*t**3)/(1*(1-material.v**2))

    def analysis(self, p, m_max_evaluation, n_max_evaluation,
                 evaluationPoint, converGraph=False, relError=5e-2,
                 printProcess=False, folderName='fig', converGraphColor='b',
                 converGraphName='Convergência', converGraphLegend=None):

                 w = self.functionAnalysis(p, m_max_evalutation,
                                            n_max_evaluation, 100, 'w')

                 error_w = [w(evaluationPoint[0], evaluationPoint[1])]
                 
                 return w, mx, my, mxy, Qx, Qy, Rx, Ry

    def functionAnalysis(self, p, m_max, n_max, num_int, interestFunction):
        """Cria uma função para analisar a deflexção para dada condição."""

        def w(x, y):
            """Função deslocamento transversal."""
            # Resolução principal do problema
            sum = 0
            for m in range(1, m_max+1):
                for n in range(1, n_max+1):
                    sum_Qmn = 0
                    for i in range(num_int+1):
                        for j in range(num_int+1):
                            sum_Qmn += p.subs([(x_s, self.a*i/num_int), \
                            (y_s, self.b*j/num_int)])*sin(m*np.pi*i/num_int)* \
                            sin(n*np.pi*j/num_int)
                    Qmn = 4*sum_Qmn/(self.a*self.b)
                    num = np.sin(np.pi*m*x/self.a) * \
                        np.sin(np.pi*n*y/self.b) * Qmn
                    den = ((m**2)/(self.a)**2 + (n**2)/(self.b**2))**2
                    sum += num/den
            return sum/((np.pi**4)*self.D)

        def m(w_all):

            return 0

        def Q(mx, my, mxy):

            return 0

        def R(x, y):

            return 0

        if interestFunction == 'w':
            func = w
        elif interestFunction == 'm':
            func = m
        elif interestFunction == 'Q':
            func = Q
        elif interestFunction == 'R':
            func = R
        return func


def makeFolder(name):
    """Função para criar a pasta de figuras se não existente."""
    if os.path.isdir(name):
        return True
    else:
        os.mkdir(name)


# ---------------------- CONSTANTES GERAIS DO PROBLEMA ----------------------
a = 1  # [m] Tamanho da placa em x
b = 1  # [m] Tamanho da placa em y
t = 0.01  # [m] Espessura da placa
p0 = 400e3  # [N/m²] Carga transversal distribuída sob a placa
x_s, y_s = symbols('x_s y_s') # Variáveis para descrever p(x, y)
p = p0*sin(x_s*np.pi/a)*sin(y_s*np.pi/b) # Carga distribuída p(x, y)
x1 = a/2 # [m] Posição em x do centro da área de aplicação da carga
y1 = b/2 # [m] Posição em y do centro da área de aplicação da carga
c = a  # [m] Largura (em x) da área de aplicação da carga
d = b  # [m] Comprimento (em y) da área de aplicação da carga
p_area = [x1, y1, c, d] # Área de aplicação da carga distribuída
xnum = 100  # Quantidade de pontos em x
ynum = 100  # Quantidade de pontos em y
M_MAX = 100  # Número máximo de avaliações de m para a convergência
N_MAX = 100  # Número máximo de avaliações de n para a convergência
P_aval = [x1, y1]  # Ponto para avaliar as deflexões (convergência e outros)

ex_num = '51'
figDirectory = 'fig' + ex_num  # Nome da pasta para salvar as imagens

# ---------------------- ANÁLISE PRINCIPAL DA DEFLEXÃO ----------------------
# Cria a placa e a pasta para salvas as figuras
myPlate = plate(a, b, t, mat.mat_exemplos5)
makeFolder(figDirectory)

# Cria os pontos pro gráfico e analisa a placa de entrada
X = np.linspace(0, a, xnum)
Y = np.linspace(0, b, ynum)
X, Y = np.meshgrid(X, Y)
w, m_conv, n_conv = myPlate.analysis(p, M_MAX, N_MAX, relError=5e-2, P_aval
                                     converGraph=True, printProcess=True,
                                     foldername = figDirectory)
plt.savefig(figDirectory + '\\91_convergencia.png')
Z = np.ndarray(shape=(xnum, ynum), dtype=float, order='F')
for i in range(xnum):
    for j in range(ynum):
        Z[i][j] = w(X[i][j], Y[i][j])

# Gráfico 3D
fig = plt.figure('Plot3D')
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, 1e3*Z, cmap=cm.coolwarm, linewidth=0,
                       antialiased=True)
plt.xlabel('Eixo x [m]')
plt.ylabel('Eixo y [m]')
ax.set_zlabel('Deflexão z [mm]')
plt.savefig(figDirectory + '\\91_deflexao3D.png')

# ------------------------------- FINALIZAÇÃO -------------------------------
plt.show()
