#-*- coding: utf-8 -*-.
'''
Created on 12/01/2016

@author: Alberto Bosch
'''
import random
import math

import Gnuplot

def generarPesos(N):
    '''
    Todos los pesos son inicializados uniformemente y son vectores unitarios
    '''    
    pesos=[]
    for x in range(1, N+1):
        w1 = random.random()
        w2 = random.random()
        w1 = w1/math.sqrt(w1**2 + w2**2)
        w2 = w2/math.sqrt(w1**2 + w2**2)
        pesos.append([w1,w2])
    return pesos

def generarPoblacion(N, D): 
    individuos = []
    f = open('all_popz.out', 'w')
    for x in range(1, N + 1):
        entradas = []        
        for y in range(1, D + 1):
            entradas.append(random.random())
        individuos.append(entradas)
        resultado = zdt3(entradas)
        f.write(str(resultado[0]) + "\t" + str(resultado[1]) + "\t" + "0.000000e+00\n")
    f.close()
    return individuos    


def generarDistancias(N, m, pesos):
    distancias = []
    for x in range(1, N+1):
        distancia = []
        for y in range(1, N+1):
            suma = 0
            for z in range(0, m):
                suma += (pesos[x-1][z]-pesos[y-1][z])**2
            distancia.append(math.sqrt(suma))
        distancias.append(distancia)
    return distancias

def generarIndicesVecinos(N, T, distancias):
    vecindad = {}
    for x in range(0, N):
        vecinos = sorted(distancias[x])[1:int(T*N)+1]
        indices = []
        for vecino in vecinos:
            indices.append(distancias[x].index(vecino))        
        vecindad[x] = indices
    return vecindad

def redondear(x, xl, xu):
    if x < xl:
        return xl
    elif x > xu:
        return xu
    else:
        return x
    
def zdt3(xreal):
    xreal = [redondear(x, 0.0, 1.0) for x in xreal]
    tmp = sum(xreal[i] for i in range(1,30))
    obj = []
    obj.append(xreal[0])
    g = 1 + ((9*tmp)/29)
    h = 1-math.sqrt(obj[0]/g)-(obj[0]/g)*math.sin(10*math.pi*obj[0])
    obj.append(g*h)
    return obj    
    
# def actualizarDominados(individuos, poblacion, EP, m):    
#     blacklist = []
#     for individuo in poblacion:
#         for objetivo in range (0, m):
#             dominanciaObjetivo = []
#             for cadaindividuo in poblacion:
#                 if individuo[objetivo] > cadaindividuo[objetivo]:
#                     dominanciaObjetivo.append(cadaindividuo[objetivo])
#         if len(dominanciaObjetivo) == m:
#             blacklist.append(poblacion.index(individuo))
#         
#     for individuo in poblacion:
#         indice = poblacion.index(individuo)
#         
#         if indice not in blacklist:
#             if indice not in EP:
#                 EP.append(poblacion[indice])
#     return EP                
                    
def actualizarDominados(poblacion, individuos, EP):
    blacklist = []
    EP = []
    for individuo in poblacion:
        for inviduo2 in poblacion:
            if inviduo2 != individuo:
                if(individuo[0]<=inviduo2[0] and individuo[1]<=inviduo2[1]):
                    indice = poblacion.index(inviduo2)
                    if indice not in blacklist:
                        blacklist.append(indice)
    for i in range(0, len(poblacion)):
        if i not in blacklist:
            EP.append(zdt3(individuos[i]))
        elif individuos[i] in EP:
            EP.delete(zdt3(individuos[i]))
            
    return EP
    
def algoritmo(N=100, G=100, m=2, T=0.2):
    ''' 
    :param N: Población
    :param m: Número de objetivos
    :param T: Tamaño de la vecindad(0-1)
    :param G: Generaciones
    :param pesos: Distribución de vectores peso
    '''
    
    # Inicialización
    pesos = generarPesos(N)
    distancias = generarDistancias(N, m, pesos)
    indices = generarIndicesVecinos(N, T, distancias)
    individuos = generarPoblacion(N, 30)
    valoraciones = [zdt3(individuo) for individuo in individuos]
    mejores_objetivos = [min(valor[i] for valor in valoraciones) for i in range(0, m)]
    EP = actualizarDominados(valoraciones, individuos, [])
    F = 0.5
    
    f = open('all_popz.out', 'a')
    g = Gnuplot.Gnuplot(persist=1)
    g.title('Generacion 1')
    g.plot(valoraciones)
    for gen in range(1, G):
        poblacion_anterior = list(individuos)
        for i in range(0, N):
            # Reproducción
            r1, r2, r3 = random.sample(indices[i], 3)
            x1, x2, x3 = poblacion_anterior[r1], poblacion_anterior[r2], poblacion_anterior[r3]

            y = list(poblacion_anterior[i])
            for componente in range(0, 30):
                if random.random() > 0.5:
                    y[componente] = x1[componente] + F*(x2[componente]-x3[componente])

                if random.random() < 0.05:
                    y[componente] = random.gauss(y[componente], 0.05)


            #Evaluación
            eval = zdt3(y)

            #Actualización de z
            mejores_objetivos = [min(a, b) for a, b in zip(mejores_objetivos, eval)]


            #Actualización de vecinos
            for j in range(0, int(T * N)):
                vector_peso = pesos[indices[i][j]]
                gte_actual = max(u * abs(v-w) for u, v,w in zip(vector_peso, valoraciones[indices[i][j]], mejores_objetivos))
                gte_nuevo = max(u * abs(v-w) for u, v,w in zip(vector_peso, eval, mejores_objetivos))
                if gte_nuevo <= gte_actual:
                    individuos[indices[i][j]] = y
                    valoraciones[indices[i][j]] = eval


        g.title('Generacion '+str(gen+1))
        g.plot(valoraciones)

        #Actualización de EP
        EP = actualizarDominados(valoraciones, individuos, EP)       
        for objetivos in valoraciones:
            f.write(str(objetivos[0]) + "\t" + str(objetivos[1]) + "\t" + "0.000000e+00\n")
    f.close()
    f = open("ga.out", "w")
    for line in EP:
        f.write(str(line[0])+"\t"+str(line[1])+"\n")
    f.close()
    f = open("plot.out", "w")
    for line in valoraciones:
        f.write(str(line[0])+"\t"+str(line[1])+"\n")
    f.close()
    print (mejores_objetivos)
        
        

    
algoritmo()
