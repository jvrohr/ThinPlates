"""Programa para resolução do exercício 9.1 do livro."""

# Importações
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import materials as mat
import os
from scipy.optimize import curve_fit


class plate:
    """Cria a classe da placa."""

    def __init__(self, a, b, t, material):
        """Método de inicialização da placa."""
        self.a = a
        self.b = b
        self.t = t
        self.D = (material.E*t**3)/(1*(1-material.v**2))

    def analysis(self, p0, N, m_max_evaluation, n_max_evaluation,
                 converGraph=False, relError=5e-2, pointType='center',
                 evaluationPoint=[0, 0], typeOfAnalysis='together',
                 printProcess=False, folderName='fig', converGraphColor='b',
                 converGraphName='Convergencia', converGraphLegend=None):
        """Cria a função para a análise da deflexão da placa."""
        print('Iniciando a análise para a placa.')

        if pointType == 'center':
            p = [self.a/2, self.b/2]
        elif pointType == 'other':
            p = evaluationPoint
        else:
            print('O valor para pointType ({}) não está definido'
                  .format(evaluationPoint))

        wEval = np.ndarray(shape=m_max_evaluation+1, dtype=float)
        if typeOfAnalysis == 'together':
            w = self.functionAnalysis(p0, N, 1, 1)
            wEval[0] = w(p[0], p[1])
            error = np.ndarray(shape=m_max_evaluation)

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
                if e_MAX == m_max_evaluation:
                    print('Não foi alcaçada convergência')
                    return 0, 0
                k += 1

            if converGraph:
                def func(x, a, b):
                    return a * np.exp(-b * x)

                plt.figure(converGraphName + '1')
                if converGraphLegend is None:
                    plt.plot(mn_values[0:k+1], 1e3*wEval[0:k+1],
                             linewidth=1.5, color=converGraphColor)
                else:
                    plt.plot(mn_values[0:k+1], 1e3*wEval[0:k+1], linewidth=1.5,
                             color=converGraphColor, label=converGraphLegend)

                plt.xlabel('$m_{max}$ = $n_{max}$')
                plt.ylabel('Deflexão $w$ [mm]')
                plt.grid(b=True, which='minor', color='gray', linestyle='-',
                         linewidth=0.1)
                plt.grid(b=True, which='major', color='gray', linestyle='-',
                         linewidth=0.1)
                plt.minorticks_on()

                plt.figure(converGraphName + '2')
                if converGraphLegend is None:
                    plt.plot(mn_values[1:k+1], error[0:k], 'o',
                             linewidth=1.5, color=converGraphColor,
                             label='Erros nas iterações')
                    popt, pcov = curve_fit(func, mn_values[1:k+1], error[0:k])
                    plt.plot(np.linspace(mn_values[1], mn_values[k+1], 100),
                             func(np.linspace(mn_values[1], mn_values[k+1],
                                  100), *popt), 'k--',
                             label='Ajuste de curva')
                else:
                    plt.plot(mn_values[1:k+1], error[0:k], 'o', linewidth=1.5,
                             color=converGraphColor, label=converGraphLegend)
                    popt, pcov = curve_fit(func, mn_values[1:k+1], error[0:k])
                    plt.plot(np.linspace(mn_values[1], mn_values[k+1], 100),
                             func(np.linspace(mn_values[1], mn_values[k+1],
                                  100), *popt), converGraphColor + '--')

                plt.xlabel('$m_{max}$ = $n_{max}$')
                plt.ylabel('Deflexão [mm]')
                plt.grid(b=True, which='minor', color='gray', linestyle='-',
                         linewidth=0.1)
                plt.grid(b=True, which='major', color='gray', linestyle='-',
                         linewidth=0.1)
                plt.minorticks_on()

            return self.functionAnalysis(p0, N, e_MAX, e_MAX), e_MAX, e_MAX

    def functionAnalysis(self, p0, N, m_max, n_max):
        """Cria uma função para analisar a deflexção para dada condição."""
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
N0 = 500e3  # [N/m] Carga lateral distribuída sob a placa
xnum = 100  # Quantidade de pontos em x
ynum = 100  # Quantidade de pontos em y
M_MAX = 100  # Número máximo de avaliações de m para a convergência
N_MAX = 100  # Número máximo de avaliações de n para a convergência
P_aval = [a/2, b/2]  # Ponto para avaliar as deflexões (convergência e outros)

pointsN = 100  # Quantidade de pontos para a análise da variação de N
Nstart = -1e6  # Valor inicial do intervalo de N para a análise
Nend = 1e6  # Valor final do intervalo de N para a análise
Npor = [0.8, 0.9, 1, 1.1, 1.2]  # Proporções de N para a análise na conv.
colors = ['r', 'b', 'k', 'y', 'g']  # Cores para diferentes gráficos

pointsAB = 101  # Quantidade de pontos para a análise da variação da razão AB
ABstart = 0.1  # Valor inicial do intervalo de A/B para a análise
ABend = 4  # Valor final do intervalo A/B para a análise

figDirectory = 'fig'  # Nome da pasta para salvar as imagens

# ---------------------- ANÁLISE PRINCIPAL DA DEFLEXÃO ----------------------
# Cria a placa e a pasta para salvas as figuras
myPlate = plate(a, b, t, mat.Al6061)
makeFolder(figDirectory)

# Cria os pontos pro gráfico e analisa a placa de entrada
X = np.linspace(0, a, xnum)
Y = np.linspace(0, b, ynum)
X, Y = np.meshgrid(X, Y)
w, m_conv, n_conv = myPlate.analysis(p0, N0, M_MAX, N_MAX, relError=5e-4,
                                     converGraph=True, printProcess=True,
                                     evaluationPoint=P_aval)
plt.figure('Convergencia1')
plt.savefig(figDirectory + '\\91_convergencia1.png')
plt.figure('Convergencia2')
plt.legend(loc='upper right')
plt.savefig(figDirectory + '\\91_convergencia2.png')

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
ax.set_zlabel('Deflexão $w$ [mm]')
plt.savefig(figDirectory + '\\91_deflexao3D.png')

# ------------- VARIAÇÃO DA CARGA DISTRIBUIDA N NA CONVERGÊNCIA -------------
N_new = [N0*i for i in Npor]
W = np.ndarray(shape=len(Npor))
for i in range(len(Npor)):
    legendName = str(Npor[i]) + ' $N_0$'
    w, m_conv, n_conv = myPlate.analysis(p0, N_new[i], M_MAX, N_MAX,
                                         converGraph=True,
                                         relError=5e-4,
                                         converGraphColor=colors[i],
                                         converGraphName='ConvergenciaN',
                                         converGraphLegend=legendName)
    W[i] = w(P_aval[0], P_aval[1])
plt.figure('ConvergenciaN1')
plt.legend(loc='upper right')
plt.savefig(figDirectory + '\\91_convergenciaN1.png')
plt.figure('ConvergenciaN2')
plt.legend(loc='upper right')
plt.savefig(figDirectory + '\\91_convergenciaN2.png')

# --------------------- VARIAÇÃO DA CARGA DISTRIBUIDA N ---------------------
print('Iniciando análise da variação da carga N.')
N_new = np.linspace(Nstart, Nend, pointsN)
W = np.ndarray(shape=pointsN)
for i in range(pointsN):
    w, m_conv, n_conv = myPlate.analysis(p0, N_new[i], M_MAX, N_MAX)
    W[i] = w(P_aval[0], P_aval[1])

# Gráfico
fig = plt.figure('VariacaoN')
plt.plot(1e-3*N_new, 1e3*W)
plt.xlabel('Força distribuída $N_0$ [kN]')
plt.ylabel('Deflexão $w$ [mm]')
plt.grid(b=True, which='minor', color='gray', linestyle='-', linewidth=0.1)
plt.grid(b=True, which='major', color='gray', linestyle='-', linewidth=0.1)
plt.minorticks_on()
plt.savefig(figDirectory + '\\91_variacaoN.png')

# -------------------------- VARIAÇÃO DA RAZÃO A/B --------------------------
print('Iniciando análise da variação da razão AB.')
ab = np.linspace(ABstart, ABend, pointsAB)
A0 = a*b  # A variação é feita mantendo a área constante
W = np.ndarray(shape=pointsAB)
for i in range(pointsAB):
    newb = np.sqrt(A0/ab[i])
    newa = A0/newb
    myPlate = plate(newa, newb, t, mat.Al6061)
    w, m_conv, n_conv = myPlate.analysis(p0, N0, M_MAX, N_MAX)
    W[i] = w(P_aval[0], P_aval[1])

# Gráfico
plt.figure('VariacaoAB')
plt.plot(ab, 1e3*W)
plt.xlabel('Razão $a/b$')
plt.ylabel('Deflexão $w$ [mm]')
plt.grid(b=True, which='minor', color='gray', linestyle='-', linewidth=0.1)
plt.grid(b=True, which='major', color='gray', linestyle='-', linewidth=0.1)
plt.minorticks_on()
plt.savefig(figDirectory + '\\91_variacaoAB.png')

# ------------------------------- FINALIZAÇÃO -------------------------------
plt.show()
