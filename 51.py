"""Programa para resolução do exercício 9.1 do livro."""

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

    def
    def analysis(self, p, m_max_evaluation, n_max_evaluation,
                 converGraph=False, relError=5e-2, evaluationPoint,
                 typeOfAnalysis='together', printProcess=False,
                 folderName='fig', converGraphColor='b',
                 converGraphName='Convergencia', converGraphLegend=None):
        """Cria a função para a análise da deflexão da placa."""
        print('Iniciando a análise para a placa.')

        wEval = np.ndarray(shape=m_max_evaluation+1, dtype=float)
        if typeOfAnalysis == 'together':
            w = self.functionAnalysis(p0, N, 1, 1)
            wEval[0] = w(evaluationPoint[0], evaluationPoint[1])

            E_MAX = 0
            k = 1
            mn_values = np.array(range(int((m_max_evaluation+1)/2)))*2 + 1
            print('\tBusca por convergência para Er = {}'.format(relError))
            if printProcess:
                print('\tm = n =  1 \t', end='')
                print('Nov aval: {:7.4}'.format((wEval[0])), end='')
                print(' Aval ant:    -    Erro:     -')
            for e_MAX in mn_values[1:-1]:
                w = self.functionAnalysis(p0, N, e_MAX, e_MAX)
                wEval[k] = w(p[0], p[1])
                new_Error = abs((wEval[k] - wEval[k-1])/wEval[k-1])
                if printProcess:
                    print('\tm = n = {:2} \t'.format(e_MAX), end='')
                    print('Nov aval: {:7.4}'.format((wEval[k])), end='')
                    print(' Aval ant: {:7.4}'.format((wEval[k-1])), end='')
                    print(' Erro: {:7.4}'.format(new_Error))
                if new_Error < relError:
                    E_MAX = e_MAX
                    print('\tConvergência em m = n = {}.\n'.format(E_MAX))
                    break
                if e_MAX == m_max_evaluation:
                    print('Não foi alcaçada convergência')
                    return 0, 0
                k += 1

            if converGraph:
                plt.figure(converGraphName)
                if converGraphLegend is None:
                    plt.plot(mn_values[0:k+1], 1e3*wEval[0:k+1], linewidth=1.5,
                             color=converGraphColor)
                else:
                    plt.plot(mn_values[0:k+1], 1e3*wEval[0:k+1], linewidth=1.5,
                             color=converGraphColor, label=converGraphLegend)

                plt.xlabel(r'$m_{max} = n_{max}$')
                plt.ylabel('Deflexão [mm]')
                plt.grid(b=True, which='minor', color='gray', linestyle='-',
                         linewidth=0.1)
                plt.grid(b=True, which='major', color='gray', linestyle='-',
                         linewidth=0.1)
                plt.minorticks_on()

            return self.functionAnalysis(p0, N, e_MAX, e_MAX), e_MAX, e_MAX

    def functionAnalysis(self, p0, N, m_max, n_max):
        """Cria uma função para analisar a deflexão para dada condição."""
        if not m_max % 2 or not n_max % 2:
            print('Erro: os valores precisam ser números inteiros ímpares.')

        def w(x, y):
            """Função deslocamento transversal."""
            # Resolução principal do problema
            sum = 0
            for m in np.array(range(int((m_max+1)/2)))*2 + 1:
                for n in np.array(range(int((n_max+1)/2)))*2 + 1:
                    num = np.sin(np.pi*m*x/self.a) * \
                        np.sin(np.pi*n*y/self.b)
                    den = m*n*(((m/self.a)**2 + (n/self.b)**2)**2 +
                               (N/self.D)*(m/(np.pi*self.a))**2)
                    sum += num/den

            return sum*16*p0/((np.pi**6)*self.D)

        return w


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
w, m_conv, n_conv = myPlate.analysis(p, M_MAX, N_MAX, relError=5e-2,
                                     converGraph=True, printProcess=True,
                                     evaluationPoint=P_aval)
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
