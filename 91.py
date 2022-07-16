"""Programa para resolução do exercício 9.1 do livro."""

# Importações
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import cm
import materials as mat

# Constantes
xnum = 100
ynum = 100

class plate:
    """Cria a classe da placa."""

    def __init__(self, a, b, t, material):
        """Método de inicialização da placa."""
        self.a = a
        self.b = b
        self.t = t
        self.D = (material.E*t**3)/(1*(1-material.v**2))

    def analysis(self, p0, N, m_max_evaluation, n_max_evaluation,
                 converGraph='False', relError=50e-3, evaluationPoint='center',
                 xEvaluation=0, yEvaluation=0, typeOfAnalysis='together',
                 printProcess='False'):
        """Cria a função para a análise da deflexão da placa."""
        print('Iniciando a análise para a placa.')

        if evaluationPoint == 'center':
            p = [self.a/2, self.b/2]
        else:
            p = [xEvaluation, yEvaluation]

        wEval = np.ndarray(shape=m_max_evaluation+1, dtype=float)
        if typeOfAnalysis == 'together':
            w = self.functionAnalysis(p0, N, 1, 1)
            wEval[0] = w(p[0], p[1])

            E_MAX = 0
            k = 1
            mn_values = np.array(range(int((m_max_evaluation+1)/2)))*2 + 1
            print('Busca por convergência para Er = {}'.format(relError))
            print('\tm = n =  1 \tNov aval: {:7.4}'.format((wEval[0])), end='')
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
                plt.figure('Convergencia')
                plt.plot(mn_values[0:k], wEval[0:k], linewidth=2)
                plt.xlabel('m_max = n_max')
                plt.ylabel('Deflexão [m]')
                plt.grid(b=True, which='minor', color='gray', linestyle='-',
                         linewidth=0.3)
                plt.grid(b=True, which='major', color='k', linestyle='-',
                         linewidth=0.8)
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


# Trocar
a = 1
b = 1

# Cria o material e a placa
myPlate = plate(a, b, 0.01, mat.Al6061)

# Parte principal
X = np.linspace(0, a, xnum)
Y = np.linspace(0, b, ynum)
X, Y = np.meshgrid(X, Y)
w, m_conv, n_conv = myPlate.analysis(4, 500e3, 101, 101)
Z = np.ndarray(shape=(xnum, ynum), dtype=float, order='F')
for i in range(xnum):
    for j in range(ynum):
        Z[i][j] = w(X[i][j], Y[i][j])


# Gráficos
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
plt.show()
