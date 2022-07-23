"""Programa para resolução do exercício 5.1 do livro."""

# Importações
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import materials as mat
import os
from sympy import *
from scipy.optimize import curve_fit
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors as mcolors

class plate:
    """Cria a classe da placa."""

    def __init__(self, a, b, t, material):
        """Método de inicialização da placa."""
        self.a = a
        self.b = b
        self.t = t
        self.v = material.v
        self.D = (material.E*t**3)/(12*(1-material.v**2))
        print('D = {}'.format(self.D))

    def convAnalysisTogether(self, p0, max_evaluation=100, converGraph=False,
                         relError=5e-2, pointType='center',
                         evaluationPoint=[0, 0], printProcess=False,
                         folderName='fig', converGraphColor='b',
                         converGraphName='Convergencia',
                         converGraphLegend=None):
        """Cria a função para a análise da deflexão da placa."""
        print('Iniciando a análise para a placa.')

        if pointType == 'center':
            p = [self.a/2, self.b/2]
        elif pointType == 'other':
            p = evaluationPoint
        else:
            print('O valor para pointType ({}) não está definido'
                  .format(evaluationPoint))

        wEval = np.ndarray(shape=max_evaluation+1, dtype=float)

        w = self.functionAnalysis(p0, 1, 1, 'w')
        wEval[0] = w(p[0], p[1])
        error = np.ndarray(shape=max_evaluation)

        E_MAX = 0
        k = 1
        mn_values = np.array(range(int((max_evaluation+1)/2)))*2 + 1
        print('\tBusca por convergência para Er = {}'.format(relError))
        if printProcess:
            print('\tm = n =  1 \t', end='')
            print('Nov aval: {:7.4}'.format((wEval[0])), end='')
            print(' Aval ant:    -    Erro:     -')
        for e_MAX in mn_values[1:-1]:
            w = self.functionAnalysis(p0, e_MAX, e_MAX, 'w')
            wEval[k] = w(p[0], p[1])
            error[k-1] = abs((wEval[k] - wEval[k-1])/wEval[k-1])
            if printProcess:
                print('\tm = n = {:2} \t'.format(e_MAX), end='')
                print('Nov aval: {:7.4}'.format((wEval[k])), end='')
                print(' Aval ant: {:7.4}'.format((wEval[k-1])), end='')
                print(' Erro: {:7.4}'.format(error[k-1]))
            if error[k-1] < relError:
                E_MAX = e_MAX
                print('\tConvergência em m = n = {}.\n'.format(E_MAX))
                break
            if e_MAX == max_evaluation:
                print('Não foi alcaçada convergência')
                return 0, 0
            k += 1

        if converGraph:
            def func(x, a, b):
                return a * np.exp(-b * x)

            plt.figure(converGraphName + '1')
            if converGraphLegend is None:
                plt.plot(mn_values[0:k+1], (1/wEval[k])*wEval[0:k+1], 'o-',
                         linewidth=1.5, color=converGraphColor)
            else:
                plt.plot(mn_values[0:k+1], (1/wEval[k])*wEval[0:k+1], 'o-',
                         linewidth=1.5, color=converGraphColor,
                         label=converGraphLegend)

            plt.xlabel('$m_{max} = n_{max}$')
            plt.ylabel('Deflexão $w$')
            plt.grid(b=True, which='minor', color='gray', linestyle='-',
                     linewidth=0.1)
            plt.grid(b=True, which='major', color='gray', linestyle='-',
                     linewidth=0.1)
            plt.minorticks_on()

        return self.functionAnalysis(p0, e_MAX, e_MAX, 'w'), e_MAX, e_MAX

    def functionAnalysis(self, p, m_max, n_max, interestFunction):
        """Cria uma função para analisar a deflexção para dada condição."""

        def w(x, y):
            """Função deslocamento transversal."""
            # Resolução principal do problema
            sum = 0
            for m in np.array(range(int((m_max+1)/2)))*2 + 2:
                for n in np.array(range(int((n_max+1)/2)))*2 + 2:
                    Qmn = (p/np.pi)*(np.sin(np.pi*(m+1))/(m+1) - np.sin(np.pi*(m-1))/(m-1))*(np.sin(np.pi*(n+1))/(n+1) - np.sin(np.pi*(n-1))/(n-1))
                    num = Qmn*np.sin(np.pi*m*x/self.a)*np.sin(np.pi*n*y/self.b)
                    den = ((m/self.a)**2 + (n/self.b)**2)**2
                    sum += num/den
            return sum/((np.pi**4)*self.D)

        def m(x, y):
            wsum = 0
            wsum2 = 0
            for m in np.array(range(int((m_max+1)/2)))*2 + 2:
                for n in np.array(range(int((n_max+1)/2)))*2 + 2:
                    Qmn = (p/np.pi)*(np.sin(np.pi*(m+1))/(m+1) - np.sin(np.pi*(m-1))/(m-1))*(np.sin(np.pi*(n+1))/(n+1) - np.sin(np.pi*(n-1))/(n-1))
                    wnum = Qmn*np.sin(np.pi*m*x/self.a)*np.sin(np.pi*n*y/self.b)
                    wnum2 = Qmn*np.cos(np.pi*m*x/self.a)*np.cos(np.pi*n*y/self.b)
                    wden = (np.pi**4)*self.D*((m/self.a)**2 + (n/self.b)**2)**2
                    wsum2 += wnum2/wden
                    wsum += wnum/wden
            d2wdx2 = -(np.pi*m/self.a)**2*wsum
            d2wdy2 = -(np.pi*n/self.b)**2*wsum
            d2wdxy = np.pi**2*n*m*wsum2/(self.a*self.b)
            mx = -self.D*(d2wdx2 + self.v*d2wdy2)
            my = -self.D*(d2wdy2 + self.v*d2wdx2)
            mxy = -self.D*(1-self.v)*d2wdxy
            return mx, my, mxy

        def Q(x, y, m_f):
            df = 1e-6
            dmxdx = (m_f(x+df, y)[0] - m_f(x-df, y)[0])/(2*df)
            dmxydy = (m_f(x, y+df)[2] - m_f(x, y-df)[2])/(2*df)
            dmxydx = (m_f(x+df, y)[2] - m_f(x-df, y)[2])/(2*df)
            dmydy = (m_f(x, y+df)[1] - m_f(x, y-df)[1])/(2*df)
            Qx = dmxdx + dmxydy
            Qy = dmxydy + dmydy
            return Qx, Qy

        def R(x, y, m_f):
            df = 1e-6
            dmxdx = (m_f(x+df, y)[0] - m_f(x-df, y)[0])/(2*df)
            dmxydy = (m_f(x, y+df)[2] - m_f(x, y-df)[2])/(2*df)
            dmxydx = (m_f(x+df, y)[2] - m_f(x-df, y)[2])/(2*df)
            dmydy = (m_f(x, y+df)[1] - m_f(x, y-df)[1])/(2*df)
            Rx = dmxdx + 2*dmxydy
            Ry = dmydy + 2*dmxydx
            return Rx, Ry

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
p = 1e3 # Carga distribuída p(x, y)
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

ex_num = '52'
figDirectory = 'fig' + ex_num  # Nome da pasta para salvar as imagens

# ---------------------- ANÁLISE PRINCIPAL DA DEFLEXÃO ----------------------
# Cria a placa e a pasta para salvas as figuras
myPlate = plate(a, b, t, mat.mat_exemplos5)
makeFolder(figDirectory)

w, m_conv, n_conv = myPlate.convAnalysisTogether(p, max_evaluation=100, converGraph=True,
                     relError=5e-2, pointType='center',
                     evaluationPoint=P_aval, printProcess=True,
                     folderName=figDirectory, converGraphColor='b',
                     converGraphName='Convergencia',
                     converGraphLegend=None)

plt.savefig(figDirectory + '\\52_convergencia.png')
m = myPlate.functionAnalysis(p, m_conv, n_conv, 'm')
Q = myPlate.functionAnalysis(p, m_conv, n_conv, 'Q')
R = myPlate.functionAnalysis(p, m_conv, n_conv, 'R')

# Cria os pontos pro gráfico e analisa a placa de entrada
X = np.linspace(0, a, xnum)
Y = np.linspace(0, b, ynum)
Rx = np.ndarray(shape=(2, ynum+2), dtype=float, order='F')
Ry = np.ndarray(shape=(2, xnum+2), dtype=float, order='F')
for j in range(ynum):
    Rx[0][j+1] = R(0, Y[j], m)[0]
    Rx[1][j+1] = R(a, Y[j], m)[0]
for i in range(xnum):
    Ry[0][i+1] = R(X[i], 0, m)[1]
    Ry[1][i+1] = R(X[i], b, m)[1]

X, Y = np.meshgrid(X, Y)
P = np.ndarray(shape=(xnum, ynum), dtype=float, order='F')
Z = np.ndarray(shape=(xnum, ynum), dtype=float, order='F')
Mx = np.ndarray(shape=(xnum, ynum), dtype=float, order='F')
My = np.ndarray(shape=(xnum, ynum), dtype=float, order='F')
Mxy = np.ndarray(shape=(xnum, ynum), dtype=float, order='F')
Qx = np.ndarray(shape=(xnum, ynum), dtype=float, order='F')
Qy = np.ndarray(shape=(xnum, ynum), dtype=float, order='F')
print('\tCalculando Deflexoes, Momentos e Esforcos Cortantes.')
for i in range(xnum):
    for j in range(ynum):
        P[i][j] = p*sin(np.pi*X[i][j]/a)*sin(np.pi*Y[i][j]/b)
        Z[i][j] = w(X[i][j], Y[i][j])
        Mx[i][j], My[i][j], Mxy[i][j] = m(X[i][j], Y[i][j])
        Qx[i][j], Qy[i][j] = Q(X[i][j], Y[i][j], m)

print("\n\nValores máximos:")
print('\tMomentos:\tMx = {} Nm\tMy = {} Nm\tMxy = {} Nm'.format(np.max(Mx),
        np.max(My), np.max(Mxy)))
print('\tCortantes:\tQx = {} N/m\tQy = {} N/m'.format(np.max(Qx),
        np.max(Qy)))
print('\tDeflexão:\tw = {} mm'.format(np.max(Z)*1e3))

# Gráfico Deflexão 3D
fig = plt.figure('Plot3D')
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, 1e3*Z, cmap=cm.coolwarm, linewidth=0,
                       antialiased=True)
plt.xlabel('Eixo x [m]')
plt.ylabel('Eixo y [m]')
ax.set_zlabel('Deflexão z [mm]')
plt.savefig(figDirectory + '\\52_deflexao3D.png')

# Gráficos Momentos 3D
# Plot de Mx
fig = plt.figure('PlotMomentos', figsize=plt.figaspect(1/3))
ax = fig.add_subplot(1, 3, 1, projection='3d')
surf = ax.plot_surface(X, Y, Mx, cmap=cm.coolwarm, linewidth=0,
                       antialiased=True)
plt.xlabel('Eixo x [m]')
plt.ylabel('Eixo y [m]')
ax.set_zlabel(r'$M_x$ [Nm]')
ax.set_title(r'$M_x$')
# Plot de My
ax = fig.add_subplot(1, 3, 2, projection='3d')
surf = ax.plot_surface(X, Y, My, cmap=cm.coolwarm, linewidth=0,
                       antialiased=True)
plt.xlabel('Eixo x [m]')
plt.ylabel('Eixo y [m]')
ax.set_zlabel(r'$M_y$ [Nm]')
ax.set_title(r'$M_y$')
# Plot de Mxy
ax = fig.add_subplot(1, 3, 3, projection='3d')
surf = ax.plot_surface(X, Y, Mxy, cmap=cm.coolwarm, linewidth=0,
                       antialiased=True)
plt.xlabel('Eixo x [m]')
plt.ylabel('Eixo y [m]')
ax.set_zlabel(r'$M_{xy}$ [Nm]')
ax.set_title(r'$M_{xy}$')

plt.savefig(figDirectory + '\\52_Momentos.png')

# Gráficos Momentos 3D
# Plot de Qx
fig = plt.figure('PlotCortantes', figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 2, 1, projection='3d')
surf = ax.plot_surface(X, Y, Qx, cmap=cm.coolwarm, linewidth=0,
                       antialiased=True)
plt.xlabel('Eixo x [m]')
plt.ylabel('Eixo y [m]')
ax.set_zlabel(r'$Q_x$ [N/m]')
ax.set_title(r'$Q_x$')
# Plot de Qy
ax = fig.add_subplot(1, 2, 2, projection='3d')
surf = ax.plot_surface(X, Y, Qy, cmap=cm.coolwarm, linewidth=0,
                       antialiased=True)
plt.xlabel('Eixo x [m]')
plt.ylabel('Eixo y [m]')
ax.set_zlabel(r'$Q_y$ [N/m]')
ax.set_title(r'$Q_y$')

plt.savefig(figDirectory + '\\52_Cortante.png')

# Gráfico Forças em 3D
fig = plt.figure('Forcas3D')
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, 0*Z, cmap=cm.coolwarm, linewidth=0,
                       antialiased=True)
surf = ax.plot_surface(X, Y, P, cmap=cm.coolwarm, linewidth=0,
                       antialiased=True, alpha = 0.7)
X = np.linspace(0, a, xnum+2)
X_rx = [np.full((ynum+2,1), 0).reshape(-1), np.full((ynum+2,1), a).reshape(-1)]
Y = np.linspace(0, b, ynum+2)
Y_rx = [np.full((xnum+2,1), 0).reshape(-1), np.full((xnum+2,1), b).reshape(-1)]

pos = [[0, a], [0, b]]
for i in range(2):
    ax.plot(X_rx[i], Y, Rx[i], 'r')
    verts = [list(zip(Y, Rx[i]))]
    poly = PolyCollection(verts, facecolors=[mcolors.to_rgba('r', alpha=0.6)])
    poly.set_alpha(0.7)
    ax.add_collection3d(poly, zs=pos[0][i], zdir='x')

    ax.plot(X, Y_rx[i], Ry[i], 'r')
    verts = [list(zip(X, Ry[i]))]
    poly = PolyCollection(verts, facecolors=[mcolors.to_rgba('r', alpha=0.6)])
    poly.set_alpha(0.7)
    ax.add_collection3d(poly, zs=pos[1][i], zdir='y')

x = [[0, 0], [a, a], [0, 0], [a, a]]
y = [[0, 0], [b, b], [b, b], [0, 0]]
z = [[0, 2*Mxy[0][0]], [0, 2*Mxy[xnum-1][ynum-1]], [0, 2*Mxy[0][ynum-1]], [0, 2*Mxy[xnum-1][0]]]

for i in range(4):
    ax.plot(x[i], y[i], z[i], 'k')

plt.xlabel('Eixo x [m]')
plt.ylabel('Eixo y [m]')
ax.set_zlabel('Forças [N]')
plt.savefig(figDirectory + '\\52_forcas3D.png')

# ------------------------------- FINALIZAÇÃO -------------------------------
plt.show()
