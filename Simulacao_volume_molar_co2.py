import numpy as np

class Questao2():
    
    def __init__(self,a,b,comp,R):
        self.a = a
        self.b = b
        self.comp = comp
        self.R = R
        """Parâmetros para avaliação de desempenho dos métodos numéricos"""
        self.numit_newton = 0
        self.numit_bisseccao = 0
        """Parâmetros para plotagem da função objetivo"""
        self.valorv = []
        self.valorg = []
    
    """Funções auxiliares"""
    
    def funcaoG(self,p,T,V):
        func = p*V + ((self.a)/(V)) - self.b*p - ((self.b*self.a)/(V**2)) - self.R*T
        return func
    def funcaoG_der(self,p,V):
        der_G = p - ((self.a)/(V**2)) + 2*((self.a*self.b)/(V**3))
        return der_G
    
    def error(self,V1,V_1):
        erro = np.absolute((V1 - V_1)/(V1))*100
        return erro
    
    """bloco de funções para cálculo de zeros de função"""
    
    def raiz_newton(self,V_ini,p,T,erro):
        funcG_0 = self.funcaoG(p,T,V_ini)
        derG_0 = self.funcaoG_der(p,V_ini)
        V1 = V_ini - ((funcG_0)/(derG_0))
        V_1 = V_ini
        erro_atual = self.error(V1,V_1)
        if self.funcaoG(p,T,V_1) == 0:
            return V_1
        if self.funcaoG(p,T,V1) == 0:
            return V1
        else:
            while erro_atual >= erro:
                self.numit_newton = self.numit_newton + 1 #desempenho
                func = self.funcaoG(p,T,V1)
                derG = derG_0 = self.funcaoG_der(p,V1)
                V_1 = V1
                V1 = V1 - ((func)/(derG))
                erro_atual = self.error(V1,V_1)
            return V1
    
    #método da bissecção
    def raiz_bisseccao(self,p,erro,V_inf,V_sup,T):
        
        func_inf = self.funcaoG(p,T,V_inf)
        func_sup = self.funcaoG(p,T,V_sup)
        if func_inf == 0:
            return V_inf
        if func_sup == 0:
            return V_sup
        elif func_inf * func_sup < 0:
            erro_e = np.absolute((V_inf - V_sup)/(V_inf))*100
            V1 = V_inf
            V_1 = V_sup
            while erro_e >= erro:
                self.numit_bisseccao = self.numit_bisseccao + 1 #desempenho
                V_novo = (V1+V_1)/2
                if self.funcaoG(p,T,V1) * self.funcaoG(p,T,V_novo) < 0 and self.funcaoG(p,T,V_1) * self.funcaoG(p,T,V_novo) > 0:
                    if self.funcaoG(p,T,V1) > 0 and self.funcaoG(p,T,V_1) < 0:
                        V_1 = V_novo
                        erro_e = np.absolute((V1 - V_1)/(V1))*100
                    elif self.funcaoG(p,T,V1) < 0 and self.funcaoG(p,T,V_1) > 0:
                        V_1 = V_novo
                        erro_e = np.absolute((V1 - V_1)/(V1))*100
                elif self.funcaoG(p,T,V_1) * self.funcaoG(p,T,V_novo) < 0 and self.funcaoG(p,T,V1) * self.funcaoG(p,T,V_novo)> 0:
                    if self.funcaoG(p,T,V_1) > 0 and self.funcaoG(p,T,V1) < 0:
                        V1 = V_novo
                        erro_e = np.absolute((V1 - V_1)/(V1))*100
                    elif self.funcaoG(p,T,V_1) < 0 and self.funcaoG(p,T,V1) > 0:
                        V1 = V_novo
                        erro_e = np.absolute((V1 - V_1)/(V1))*100
            return (V1 + V_1)/2
    
    """Bloco para exibição dos resultados de simulação"""
    
    def relatorio_newton(self,V_ini,p,T,erro):
        V_molar = self.raiz_newton(V_ini,p,T,erro)
        print('Resultados da simulação usando o método de Newton')
        print('Volume molar de '+self.comp+' a P='+str(p)+' atm '+'e T='+str(T)+' K: '+str(V_molar)+' L/mol.')
        print('-------------------------------------')
    
    def relatorio_bisseccao(self,p,erro,V_inf,V_sup,T):
        V_molar = self.raiz_bisseccao(p,erro,V_inf,V_sup,T)
        print('Resultados da simulação usando o método da bissecção')
        print('Volume molar de '+self.comp+' a P='+str(p)+' atm '+'e T='+str(T)+' K: '+str(V_molar)+' L/mol.')
        print('-------------------------------------')
    
    def desempenho(self):
        print('Desempenho dos métodos iterativos')
        print('Newton-Raphson: '+str(self.numit_newton))
        print('Bissecção: '+str(self.numit_bisseccao))
        print('-------------------------------------')

a = 3.592
b = 0.04267
comp = 'CO2'
R = 0.082054 # atmL/molK

V_ini = 5 #L/mol
p  = 10 #atm
T = 300 #K
erro = 0.0001
V_inf = 0.0001
V_sup = 10

obj = Questao2(a,b,comp,R)
obj.relatorio_newton(V_ini,p,T,erro)
obj.relatorio_bisseccao(p,erro,V_inf,V_sup,T)
obj.desempenho()