"""Arquivo para salvar as propriedades de diferentes materiais."""


class material:
    """Cria a classe material para guardar os dados do material."""

    def __init__(self, v, E, G=0):
        """Método de inicialização da classe de material."""
        self.v = v
        self.E = E
        if G == 0:
            self.G = E/(2*(1+v))
        else:
            self.G = G


# Lista de materiais disponíveis
Ti6Al4V = material(110e9, 0.31)
Al6061 = material(68e9, 0.33)
