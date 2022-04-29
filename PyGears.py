#Author: Juan Gras @CenturySturgeon

import math as mt
import numpy as np
import time
import struct
from io import StringIO
from io import BytesIO

def radToDeg(x):
    grados=x*180/mt.pi
    return grados

def degToRad(x):
    radianes = x * mt.pi/180
    return radianes

def rollAngle(rx,rb):
    teta = mt.sqrt((rx/rb)*(rx/rb) -1)
    return teta

def inv(x):
    invo = mt.tan(x) - x
    return invo

def sigma(rt,rb,dt,t,dp,ap):
    alphat = mt.acos(rb / rt)
    Tt = dt * ((t / dp) + inv(ap) - inv(alphat))
    sigm = Tt / rt
    return sigm


#region Funciones para los helicoidales
def helixturns(ah,height,diameter):
    ph = np.pi * diameter * (np.cos(ah) / np.sin(ah))
    numberofturns = height / ph
    return numberofturns

def lengthofhelix(ah,height,diameter):
    ph = np.pi*diameter*(np.cos(ah)/np.sin(ah))
    numberofturns = height/ph
    print("#ofturns",numberofturns)
    c = np.pi * diameter
    length = numberofturns*np.sqrt((c)**2+(ph)**2)
    return length

def insertodientes22(lista,az,vertijuan,facejuan,vertijuan3,facejuan3,vertieterno,faceterno,anchosep):
    ztotal = sum(lista)
    # lista = concatenada

    if lista[-1] > 0:
        a = 0
    else:
        a = 1
    for k in range(1, len(lista) - a):
        numgears = int(len(vertijuan) / len(vertieterno))
        thetax = az * numgears
        rot = np.array([[np.cos(thetax), -np.sin(thetax), 0], [np.sin(thetax), np.cos(thetax), 0], [0, 0, 1]])
        #vertijuan = np.matmul(vertijuan, rot)
        vertijuan3 = vertijuan3 + np.array([0,0,anchosep*numgears])
        vertijuan3 = np.matmul(vertijuan3, rot)
        for i in range(1, lista[k] + 1):
            theta = az * 2 ** (i - 1)
            rot = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
            f2 = facejuan3 + np.array([len(vertijuan3), len(vertijuan3), len(vertijuan3)])
            facejuan3 = np.append(facejuan3, f2, axis=0)
            v2 = np.matmul(vertijuan3, rot)
            v2 = v2 + np.array([0,0,anchosep*2**(i-1)])
            vertijuan3 = np.append(vertijuan3, v2, axis=0)

        #print(len(vertijuan) / len(vertieterno))
        x = int(len(vertijuan) / len(vertieterno))
        fx = facejuan3 + np.array([len(vertieterno) * x, len(vertieterno) * x, len(vertieterno) * x])
        vertijuan = np.append(vertijuan, vertijuan3, axis=0)
        facejuan = np.append(facejuan, fx, axis=0)
        vertijuan3 = vertieterno
        facejuan3 = faceterno

    if lista[-1] == 0:
        numgears = int(len(vertijuan) / len(vertieterno))
        theta = az * numgears
        rot = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
        vertijuan3 = np.matmul(vertijuan3, rot)
        vertijuan3 = vertijuan3+np.array([0,0,anchosep*numgears])
        vertijuan = np.append(vertijuan, vertijuan3, axis=0)
        x = int(len(facejuan) / len(faceterno))
        fx = facejuan3 + +np.array([len(vertieterno) * x, len(vertieterno) * x, len(vertieterno) * x])
        facejuan = np.append(facejuan, fx, axis=0)



    return vertijuan, facejuan

def insertodientes12(lista,az,vertijuan,facejuan,anchosep):
    #lista = concatenada
    for k in range (0,1):
        for i in range(1, lista[k]+1):
            theta = az*2**(i-1)
            rot = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
            facejuan2 = facejuan + np.array([len(vertijuan), len(vertijuan), len(vertijuan)])
            facejuan = np.append(facejuan, facejuan2, axis=0)
            vertijuan2 = np.matmul(vertijuan, rot)
            vertijuan2 = vertijuan2 + np.array([0,0,anchosep*2**(i-1)])
            vertijuan = np.append(vertijuan, vertijuan2, axis=0)
        pass

    return vertijuan, facejuan

def rotateZ(matrix,angle):
    theta = angle
    rot = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    rotatedmatx = np.matmul(matrix,rot)
    return rotatedmatx

def MirrorX(matrix,angle):
    theta = angle
    rot = np.array([[1, 0, 0],[0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
    rot2 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
    rotatedmatx1 = np.matmul(matrix, rot)
    rotatedmatx = np.matmul(rotatedmatx1, rot2)
    return rotatedmatx

#endregion


#Conversor a binerio es extremadamente importante. Crea una lista con los exponentes para 2^x de manera que controla la cantidad de dientes
def conversoraBinario(x):
    binario = []
    concatenada =[]
    n = []
    a = True
    v=x
    i=0
    while a:
        q = 2**i
        residuo = v % 2
        num = (2**i)*residuo
        if residuo!=0:
            concatenada.append(num)
            n.append(i)
        v=int(v/2)
        binario.append(residuo)
        if v==0:
            a=False
        i = i+1
    reversed_binario = binario
    binario = list(reversed(binario))
    concatenada = list(reversed(concatenada))
    reversed_n = list(reversed(n))
    return reversed_n

def normal_to_surface(v1,v2,v3):
    A = v2 - v1
    B = v3 - v1
    v = np.cross(A, B)
    Normal = v / np.linalg.norm(v)
    return Normal

#Creo el inserto del primer grupo de dientes (recuerda que insertodientes1 afecta directamente la matriz final mientras que insertodientes2 no lo hace hasta el final
def normales1(lista,az,vertices):
    for k in range (0,1):
        for i in range(1, lista[k]+1):
            theta = -az * 2 ** (i - 1)
            rot = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
            vertijuan2 = np.matmul(vertices, rot)
            vertices = np.append(vertices, vertijuan2, axis=0)
        pass

    return vertices

def normales2(lista,az,vertices,vertices2,vertieterno):
    if lista[-1]>0:
        a=0
    else:
        a=1
    for k in range(1,len(lista)-a):
        theta = -az * 2 **lista[k]
        rot = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
        vertices = np.matmul(vertices,rot)
        for i in range(1, lista[k]+1):
            theta = -az * 2 ** (i - 1)
            rot = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
            v2 = np.matmul(vertices2, rot)
            vertices2 = np.append(vertices2, v2, axis=0)

        vertices = np.append(vertices, vertices2, axis=0)
        vertices2 = vertieterno

    if lista[-1]==0:
        theta = az
        rot = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
        vertices2 = np.matmul(vertices2, rot)
        vertices = np.append(vertices, vertices2, axis=0)

    return vertices

#Creo el inserto del primer grupo de dientes (recuerda que insertodientes1 afecta directamente la matriz final mientras que insertodientes2 no lo hace hasta el final
def insertodientes1(lista,az,vertices,faces,vertices2,faces2,vertieterno,faceterno):
    ztotal = sum(lista)
    #lista = concatenada
    for k in range (0,1):
        for i in range(1, lista[k]+1):
            theta = -az * 2 ** (i - 1)
            rot = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
            facejuan2 = faces + np.array([len(vertices), len(vertices), len(vertices)])
            faces = np.append(faces, facejuan2, axis=0)
            vertijuan2 = np.matmul(vertices, rot)
            vertices = np.append(vertices, vertijuan2, axis=0)
        pass
    return vertices, faces

def insertodientes2(lista,az,vertices,faces,vertijuan3,facejuan3,vertieterno,faceterno):
    ztotal = sum(lista)
    #lista = concatenada
    if lista[-1]>0:
        a=0
    else:
        a=1
    for k in range(1,len(lista)-a):
        theta = -az * 2 **lista[k]
        rot = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
        vertices = np.matmul(vertices,rot)
        for i in range(1, lista[k]+1):
            theta = -az * 2 ** (i - 1)
            rot = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
            f2 = facejuan3 + np.array([len(vertijuan3), len(vertijuan3), len(vertijuan3)])
            facejuan3 = np.append(facejuan3, f2, axis=0)
            v2 = np.matmul(vertijuan3, rot)
            vertijuan3 = np.append(vertijuan3, v2, axis=0)

        x = int(len(vertices) / len(vertieterno))
        fx = facejuan3 + np.array([len(vertieterno) * x, len(vertieterno) * x, len(vertieterno) * x])
        vertices = np.append(vertices, vertijuan3, axis=0)
        faces = np.append(faces, fx, axis=0)
        vertijuan3 = vertieterno
        facejuan3 = faceterno

    if lista[-1]==0:
        theta = az
        rot = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
        vertijuan3 = np.matmul(vertijuan3, rot)
        vertices = np.append(vertices, vertijuan3, axis=0)
        x = int(len(faces) / len(faceterno))
        fx = facejuan3 + +np.array([len(vertieterno) * x, len(vertieterno) * x, len(vertieterno) * x])
        faces = np.append(faces, fx, axis=0)


    return vertices, faces


#region Funciones para los agujeros
def cuadrado(vertices, faces, normals, z, radiopoligono, a1, anchoeng, rcircunscrito):
    if radiopoligono != -1 and rcircunscrito <= radiopoligono - .15:
        jx = []
        jy = []
        ix = []
        iy = []
        numcarasorig = len(faces)
        for i in range(0, z):
            jx.append((radiopoligono) * np.cos(i * 2 * np.pi / z))
            jy.append((radiopoligono) * np.sin(i * 2 * np.pi / z))
        for i in range(0, 4):
            ix.append(rcircunscrito * np.cos(i * 2 * np.pi / 4))
            iy.append(rcircunscrito * np.sin(i * 2 * np.pi / 4))
        # vertmax es el número máximo de vértices que puede conectar un punto del hexagono con los del polígono exterior
        vertmax = int((z - 1) / 4) + 2
        # cvgrandes indica qué puntos del hexágono conectan con la cantidad vertmax con vértices del polígono exterior
        cvgrandes = z - int(z / 4) * 4
        # patrones es una lista que indica qué vertices del hexágono conectan con la cantidad vertmax (si cvgrandes=0 entonces todos los vertices del hexágono tienen vertmax)
        patrones = [[0, 1, 2, 3], [0], [0, 2], [0, 1, 2], [0, 1, 2, 3]]

        # Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        # Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        # Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx[-1], jy[-1], a1]], axis=0)
        qj = len(vertices) - 1
        for i in range(0, len(jx)):
            vertices = np.append(vertices, [[jx[i], jy[i], a1]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], a1]], axis=0)
        qi = len(vertices) - 1
        for i in range(0, len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], a1]], axis=0)
        qx = len(vertices) - 1

        #region "Poligon´s inferior faces"
        a = 0
        # Ciclo for que va entre todos los puntos del hexágono
        for i in range(0, 4):
            if i in patrones[cvgrandes]:
                # Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi + 1) - 1 + i, (qi + 1) + i, (qj + 1) - 1 + a]], axis=0)
                for k in range(0, vertmax - 1):
                    # Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qj + 1) - 1 + a + k, (qi + 1) + i, (qj + 1) + a + k]], axis=0)
                # Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                # Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a = a + vertmax - 1
            else:
                # Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi + 1) - 1 + i, (qi + 1) + i, (qj + 1) - 1 + a]], axis=0)
                for m in range(0, vertmax - 2):
                    # Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qj + 1) - 1 + a + m, (qi + 1) + i, (qj + 1) + a + m]], axis=0)
                a = a + vertmax - 2
        #endregion

        # Crea el hexágono en la cara superior
        jx = []
        jy = []
        ix = []
        iy = []
        for i in range(0, z):
            jx.append((radiopoligono) * np.cos(i * 2 * np.pi / z))
            jy.append((radiopoligono) * np.sin(i * 2 * np.pi / z))
        for i in range(0, 4):
            ix.append(rcircunscrito * np.cos(i * 2 * np.pi / 4))
            iy.append(rcircunscrito * np.sin(i * 2 * np.pi / 4))

        # Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        # Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        # Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx[-1], jy[-1], anchoeng]], axis=0)
        qj2 = len(vertices) - 1
        for i in range(0, len(jx)):
            vertices = np.append(vertices, [[jx[i], jy[i], anchoeng]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], anchoeng]], axis=0)
        qi2 = len(vertices) - 1
        for i in range(0, len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], anchoeng]], axis=0)
        qx = len(vertices) - 1
        a = 0
        # Ciclo for que va entre todos los puntos del hexágono
        for i in range(0, 4):
            if i in patrones[cvgrandes]:
                # Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi2 + 1) + i, (qi2 + 1) - 1 + i, (qj2 + 1) - 1 + a]], axis=0)
                for k in range(0, vertmax - 1):
                    # Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qi2 + 1) + i, (qj2 + 1) - 1 + a + k, (qj2 + 1) + a + k]], axis=0)
                # Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                # Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a = a + vertmax - 1
            else:
                # Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi2 + 1) + i, (qi2 + 1) - 1 + i, (qj2 + 1) - 1 + a]], axis=0)
                for m in range(0, vertmax - 2):
                    # Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qi2 + 1) + i, (qj2 + 1) - 1 + a + m, (qj2 + 1) + a + m]], axis=0)
                a = a + vertmax - 2
        # Creo las caras laterales para el hexágono
        pisep = np.linspace(0, 2 * np.pi, 5)
        xi = rcircunscrito * np.cos(pisep)
        yi = rcircunscrito * np.sin(pisep)

        qp = len(vertices) - 1
        for i in range(0, len(pisep)):
            vertices = np.append(vertices, [[xi[i], yi[i], a1]], axis=0)
        qp1 = len(vertices) - 1
        for i in range(0, len(pisep)):
            vertices = np.append(vertices, [[xi[i], yi[i], anchoeng]], axis=0)
        for i in range(0, len(pisep)):
            faces = np.append(faces, [[qp + i + 1, qp + i, qp1 + i]], axis=0)
            faces = np.append(faces, [[qp + i + 1, qp1 + i, qp1 + i + 1]], axis=0)

        # Fixed code for normals
        for i in range(numcarasorig, len(faces)):
            we = faces[i]
            v1 = vertices[we[0]]
            v2 = vertices[we[1]]
            v3 = vertices[we[2]]
            normal = normal_to_surface(v1, v2, v3)
            normals = np.append(normals, [normal], axis=0)
    else:
        print("Radio del polígono circumscrito es muy grande, pruebe con un valor <" + str(radiopoligono))
    return vertices, faces, normals

def hexagono(vertices, faces, normals, z, radiopoligono, a1, anchoeng, rcircunscrito):
    if radiopoligono != -1 and rcircunscrito <= radiopoligono-.15:
        jx=[]
        jy=[]
        ix=[]
        iy=[]
        numcarasorig = len(faces)
        for i in range(0,z):
            jx.append((radiopoligono)*np.cos(i*2*np.pi/z))
            jy.append((radiopoligono)*np.sin(i*2*np.pi/z))
        for i in range(0,6):
            ix.append(rcircunscrito* np.cos(i*2*np.pi/6))
            iy.append(rcircunscrito* np.sin(i*2*np.pi/6))
        #vertmax es el número máximo de vértices que puede conectar un punto del hexagono con los del polígono exterior
        vertmax =int((z-1)/6)+2
        #cvgrandes indica qué puntos del hexágono conectan con la cantidad vertmax con vértices del polígono exterior
        cvgrandes = z-int(z/6)*6
        #patrones es una lista que indica qué vertices del hexágono conectan con la cantidad vertmax (si cvgrandes=0 entonces todos los vertices del hexágono tienen vertmax)
        patrones= [[0,1,2,3,4,5], [0], [0,3], [0,2,4], [0,2,3,5], [0,1,2,3,4]]

        #Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        #Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        #Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx[-1], jy[-1], a1]], axis=0)
        qj = len(vertices)-1
        for i in range(0,len(jx)):
            vertices = np.append(vertices, [[jx[i], jy[i], a1]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], a1]], axis=0)
        qi = len(vertices)-1
        for i in range(0,len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], a1]], axis=0)
        qx= len(vertices)-1

        #region "Inferior faces of the polygon"
        a=0
        #Ciclo for que va entre todos los puntos del hexágono
        for i in range(0,6):
            if i in patrones[cvgrandes]:
                #Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi+1)-1+i, (qi + 1) + i, (qj+1)-1+a]], axis=0)

                for k in range(0,vertmax-1):
                    #Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qj+1)-1+a+k, (qi + 1) + i, (qj+1)+a+k]], axis=0)
                #Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                #Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a=a+vertmax-1
            else:
                #Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi+1)-1+i, (qi + 1) + i, (qj+1)-1+a]], axis=0)

                for m in range(0,vertmax-2):
                    #Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qj+1)-1+a+m, (qi + 1) + i, (qj+1)+a+m]], axis=0)

                a=a+vertmax-2
        #endregion

        #Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        #Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        #Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx[-1], jy[-1], anchoeng]], axis=0)
        qj2 = len(vertices)-1
        for i in range(0,len(jx)):
            vertices = np.append(vertices, [[jx[i], jy[i], anchoeng]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], anchoeng]], axis=0)
        qi2 = len(vertices)-1
        for i in range(0,len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], anchoeng]], axis=0)
        qx= len(vertices)-1
        a=0
        #Ciclo for que va entre todos los puntos del hexágono
        for i in range(0,6):
            if i in patrones[cvgrandes]:
                #Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi2+1)+i,(qi2+1)-1+i,(qj2+1)-1+a]], axis=0)
                for k in range(0,vertmax-1):
                    #Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qi2+1)+i,(qj2+1)-1+a+k,(qj2+1)+a+k]], axis=0)
                #Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                #Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a=a+vertmax-1
            else:
                #Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi2+1)+i,(qi2+1)-1+i,(qj2+1)-1+a]], axis=0)
                for m in range(0,vertmax-2):
                    #Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qi2+1)+i,(qj2+1)-1+a+m,(qj2+1)+a+m]], axis=0)
                a=a+vertmax-2
        #Creo las caras laterales para el hexágono
        pisep=np.linspace(0,2*np.pi,7)
        xi = rcircunscrito*np.cos(pisep)
        yi = rcircunscrito*np.sin(pisep)

        qp=len(vertices)-1
        for i in range(0,len(pisep)):
            vertices = np.append(vertices, [[xi[i], yi[i], a1]], axis=0)
        qp1=len(vertices)-1
        for i in range(0,len(pisep)):
            vertices = np.append(vertices, [[xi[i], yi[i], anchoeng]], axis=0)
        for i in range(0,len(pisep)):
            faces= np.append(faces, [[qp+i+1, qp+i, qp1+i]], axis=0)
            faces = np.append(faces, [[qp+i+1, qp1+i, qp1+i+1]], axis=0)

        #Code for normals
        for i in range(numcarasorig, len(faces)):
            we = faces[i]
            v1 = vertices[we[0]]
            v2 = vertices[we[1]]
            v3 = vertices[we[2]]
            normal = normal_to_surface(v1, v2, v3)
            normals = np.append(normals, [normal], axis=0)
    else:
        print("Radio del polígono circumscrito es muy grande, pruebe con un valor <" + str(radiopoligono))
    return vertices, faces, normals

#para dr y dh
def circulo(vertices,faces,normals,z,radiopoligono,a1,anchoeng,rcircunscrito):
    print("caras niciales",len(faces))
    #el circle_points es porque z es 10, pero esto debe quedar de forma que x/z siempre de un entero
    if radiopoligono != -1 and rcircunscrito <= radiopoligono-.15:


        #ATENCION !!!! TIENES QUE
        # MODFICAR CIRCLE POINTS PARA QUE SEA UNA VARABLE MAS EFCIENTE A LA HORA DE SELECCIONAR LA CANTIDAD DE
        # CARAS PARA EL CIRCULO
        circle_points = int(300/z)*z
        print(circle_points)
        circle_points_per_vertex = circle_points / z
        print(circle_points_per_vertex, mt.floor(circle_points_per_vertex / 2))
        #circle_points_per_vertex=0
        half_angle = (circle_points_per_vertex * 1 * np.pi / circle_points)
        numcarasorig = len(faces)








        #ATENCION !!!!!
        #LO QUE DEJASTE PENDIENTE:
        # El circle_points tiene que cambiar a una variable cuyo valor dependa del numero de dientes porque si no, va a fallar
        # esta nueva variable no la declares dentro de las funciones porque AMBAS FUNCIONES LA USAN
        # Despues, sustituye todos los de circle_points por ya una variable calculada (no tiene que ser la mamada) e intenta cuando z=11,13,16









        # Para esta función estoy usando "mal" a1 porque esa se planeaba para que la cara del agujero no estuviera al mismo nivel que las caras inefirores (como otro agujero o que estuviera sumida)
        # a1 la voy a usar para decirle si son las caras de abajo o la de arriba
        def base_faces(vertices,faces,z,radiopoligono,rcircunscrito,a1,normals_downwards):
            #x=normals_downwards
            jx = []
            jy = []
            ix = []
            iy = []
            for i in range(0, z):
                jx.append((radiopoligono) * np.cos(i * 2 * np.pi / z))
                jy.append((radiopoligono) * np.sin(i * 2 * np.pi / z))
            ## agruego otra vez cuando el angulo=0 para que en el loop al dar la vuelta llegue al ultimo vertice del circulo (que es el inicial porque 2pi=0)
            jx.append((radiopoligono) * np.cos(0 * 2 * np.pi / z))
            jy.append((radiopoligono) * np.sin(0 * 2 * np.pi / z))
            for i in range(0, circle_points):
                ix.append(rcircunscrito * np.cos((i * 2 * np.pi / circle_points) - half_angle))
                iy.append(rcircunscrito * np.sin((i * 2 * np.pi / circle_points) - half_angle))
            ## agruego otra vez cuando el angulo=0 para que en el loop al dar la vuelta llegue al ultimo vertice del circulo (que es el inicial porque 2pi=0)
            ix.append(rcircunscrito * np.cos((0 * 2 * np.pi / circle_points) - half_angle))
            iy.append(rcircunscrito * np.sin((0 * 2 * np.pi / circle_points) - half_angle))

            carasporvt = int(circle_points / z)
            #print("casdñfljks", carasporvt)

            qj = len(vertices) - 0
            for i in range(0, len(jx)):
                vertices = np.append(vertices, [[jx[i], jy[i], a1]], axis=0)
            # ***
            vertices = np.append(vertices, [[ix[-1], iy[-1], a1]], axis=0)
            qi = len(vertices) - 0
            for i in range(0, len(ix)):
                vertices = np.append(vertices, [[ix[i], iy[i], a1]], axis=0)
            qx = len(vertices) - 1

            for i in range(0, z):
                for j in range(0, carasporvt +0):
                    x = i * z + j +2*i
                    x = int(i*circle_points_per_vertex + j)
                    if normals_downwards:
                        faces = np.append(faces, [[qj + i, qi + x, qi + 1 + x]], axis=0)
                    else:
                        faces = np.append(faces, [[qi + x, qj + i, qi + 1 + x]], axis=0)
                if normals_downwards:
                    faces = np.append(faces, [[qj + i + 1, qj + i, qi + 1 + x]], axis=0)
                else:
                    faces = np.append(faces, [[qj + i, qj + 1 + i, qi + 1 + x]], axis=0)
            print(len(vertices))
            return vertices, faces, qi

        vf = base_faces(vertices, faces, z, radiopoligono, rcircunscrito, a1, True)
        vertices = vf[0]
        faces = vf[1]
        qi1 = vf[2]
        vf = base_faces(vertices, faces, z, radiopoligono, rcircunscrito, anchoeng, False)
        vertices = vf[0]
        faces = vf[1]
        qi2 = vf[2]

        def lateral_faces(faces, qi1, qi2, z):
            #   NOTA IMPORTANTE: CAMBIA EL PINCHE circle_points !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            for i in range(0, circle_points):
                faces = np.append(faces, [[qi1 + i, qi2+i, qi1 + i + 1]], axis=0)
                faces = np.append(faces, [[qi1 + i + 1, qi2 + i, qi2 + i + 1]], axis=0)
            return faces
        faces = lateral_faces(faces, qi1, qi2, z)

        # Code for normals
        for i in range(numcarasorig, len(faces)):
            we = faces[i]
            v1 = vertices[we[0]]
            v2 = vertices[we[1]]
            v3 = vertices[we[2]]
            normal = normal_to_surface(v1, v2, v3)
            normals = np.append(normals, [normal], axis=0)

        #print(faces)
        #print("juaaaaaan", qj,qi,qx)
        #print("caras finales", len(faces))
    else:
        print("Radio del polígono circumscrito es muy grande, pruebe con un valor <" + str(radiopoligono))
    return vertices, faces, normals


#ATENCION IMBECIL !!!
#FALTA
#ARREGLAR
#LAS NORMALES !! de las caras de abajo o arriba de keyway
#MEJOR CAMBIA LAS FUNCIONES DEL SPUR GEAR POR LAS DEL DH

def keyway(vertices, faces, normals, z, radiopoligono, anchoeng, r1, alturaC, anchoC):
    # Crea el keyhole
    numcarasorig = len(faces)
    # l1 es el poligono interior y l2 el exterior
    l2 = z
    if z <= 10:
        l1 = 20 * z
    elif 10 < z <= 20:
        l1 = 10 * z
    elif 20 < z <= 50:
        l1 = 5 * z
    elif 50 < z <= 100:
        l1 = 2 * z
    else:
        l1 = z

    # caracteristicas del keyway
    altura = alturaC
    ancho = anchoC

    # radios ocasionados por el keyway
    #r1 = 10
    r2 = radiopoligono

    # puntos del keyway
    px = r1 + altura
    py = ancho / 2
    prx = np.sqrt((r1 ** 2) - ((ancho / 2) ** 2))
    pry = ancho / 2

    pisep = np.linspace(0, 2 * np.pi, l1 + 1)
    pisep2 = np.linspace(0, 2 * np.pi, l2 + 1)

    rx1 = r1 * np.cos(pisep)
    ry1 = r1 * np.sin(pisep)

    rx2 = r2 * np.cos(pisep2)
    ry2 = r2 * np.sin(pisep2)

    # Encuentro el primer vert del pe que está fuera del cuadro
    for i in range(1, len(rx2)):
        a = rx2[i]
        b = ry2[i]
        if a > 0 and b > ancho / 2:
            item = i
            break
    for i in range(1, len(rx2)):
        a = rx2[i]
        b = ry2[i]
        if a > 0 and b < 0 and b >= -ancho / 2:
            item2 = i - 1
            break
    # Encuentro los vertices del pi que están fuera del cuadro
    for i in range(1, len(rx1)):
        a = rx1[i]
        b = ry1[i]
        if a > 0 and b > ancho / 2:
            item3 = i
            break
    for i in range(1, len(rx1)):
        a = rx1[i]
        b = ry1[i]
        if a > 0 and b < -ancho / 2:
            item4 = i

    k = int(l1 / l2)
    a = int(k / 2)
    # qq es para el polígono exterior y qp para el interior
    qq = []
    qp = []
    qk = []
    qk2 = []
    tr = True
    u = a
    for i in range(item, item2 + 1):
        for j in range(u, u + k + 1):
            if j >= item3 and j <= item4:
                qq.append(i)
                qp.append(j)
                # qk es una lista para encontrar el primer vertice del polígono exterior que conecta con uno del interior ya que item2 no necesariamente sera el 1ero
                # con len(qk) puedo saber cuantas lineas tiene: si es 0 o 1 tiene una linea, para los demás numeros tiene la cantidad que indiquen
                if tr:
                    qk.append(i)
                    qk2.append(j)
        tr = False
        u = u + k

    # meto el primer punto del pe a parte porque casi siempre su #de triangulos no es igual al de los demás
    k0 = len(vertices) - 1
    vertices = np.append(vertices, [[rx2[qq[0]], ry2[qq[0]], 0]], axis=0)
    k1 = len(vertices) - 1
    # Para los puntos del pe empiezo en +1 para que no se haga un desmadre ya que el primer grupo suele tener menos separaciones
    for i in range(qq[0] + 1, qq[-1] + 1):
        vertices = np.append(vertices, [[rx2[i], ry2[i], 0]], axis=0)
    k2 = len(vertices) - 1

    # ppi1 es el primer punto del pi en el que se empiezan a generar las separaciones uniformes (o sea todos los triangulos verdes menos el primer y ultimo grupo)
    if len(qk) <= 1:
        ppi1 = 0
    else:
        ppi1 = len(qk) - 1
    # Notese que empieza en ppi y este puede ser 0 cuando el primer grupo es sólo una línea
    for i in range(qp[0] + ppi1, qp[-1] + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k3 = len(vertices) - 1

    # Agrego el primer grupo de triángulos verdes
    for i in range(qp[0], qp[0] + ppi1 + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k4 = len(vertices) - 1

    # Agrego los triángulos sempiternos (laterales magenta con rojo)
    if ppi1 == 0:
        vu = k2
    else:
        vu = k3
    vertices = np.append(vertices, [[px, py, 0]], axis=0)
    vertices = np.append(vertices, [[px, -py, 0]], axis=0)
    k5 = len(vertices) - 1
    # añade los triángulos magenta restantes
    vertices = np.append(vertices, [[prx, pry, 0]], axis=0)
    k6 = len(vertices) - 1

    for i in range(0, qp[0] + 1):
        if ry1[i] > py:
            vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k7 = len(vertices) - 1
    k76 = k7 - k6
    cte = 0
    if k76 != 1:
        cte = 1

    # Añade los triangulos magenta inferiores restantes
    vertices = np.append(vertices, [[prx, -pry, 0]], axis=0)
    k8 = len(vertices) - 1
    for i in range(qp[-1], len(rx1)):
        if ry1[i] < -py:
            vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k9 = len(vertices) - 1
    k97 = k9 - k8

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(0, qq[0] + 1):
        vertices = np.append(vertices, [[rx2[i], ry2[i], 0]], axis=0)
    k10 = len(vertices) - 1

    for i in range(qq[-1], len(rx2)):
        vertices = np.append(vertices, [[rx2[i], ry2[i], 0]], axis=0)
    k11 = len(vertices) - 1

    u = 1
    # Puntos del pi sobrantes
    residuo = (k3 - k2 - 1) % k
    iterac = k3 - k2 - residuo
    # Crea los tríangulos verdes
    for i in range(1, z):
        for j in range(u + 0, u + k + 0):
            if j < (iterac):
                faces = np.append(faces, [[k1 + i, k2 + j + 1, k2 + j]], axis=0)
                p = 0
            if j < (k3 - k2):
                v = u
            else:
                break
        u = u + k
    # En caso de que haya triangulos del pi que sobren ya que en el ultimo punto usado del pe sus triángulos estan incompletos, se ejecuta lo siguiente
    # para conectarlos con el punto del pe anterior
    x = qq[-1] - qq[0]
    if residuo != 0:
        for i in range(0, residuo):
            faces = np.append(faces, [[k1 + x, k2 + i + v+1, k2 + i + v]], axis=0)

    # Creo el primer grupo de triángulos verdes
    for i in range(qq[0], qq[0] + 1):
        for j in range(1, 1 + ppi1):
            faces = np.append(faces, [[k0 + i, k3 + j + 1, k3 + j]], axis=0)
            p = 0
            if ppi1 == 0:
                print("Maldita sea Jimbo, el comunismo no funciona !")

    # Creo los complementos de los triangulos verdes para cerrar
    for i in range(1, qq[-1] - qq[0] + 1):
        faces = np.append(faces, [[k0 + i, k0 + i + 1, k2 + (i - 1) * (k) + 1]], axis=0)

    # Creo los triángulos sempiternos (laterales magenta con rojo)
    faces = np.append(faces, [[k0 + 1, vu + 1, k4 + 1]], axis=0)
    faces = np.append(faces, [[k2, k4 + 2, k3]], axis=0)

    #Agrego los triángulos magenta positivos
    for i in range(1, k7 - k6 + cte):
        faces = np.append(faces, [[k4 + 1, k5 + i + 1, k5 + i]], axis=0)

    # creo el sello para la esquina izquierda del cuadro
    if k7 - k6 == 1:
        faces = np.append(faces, [[k5 + 1, k4 + 1, vu + 1]], axis=0)

    # Añade los triangulos magenta inferiores restantes
    for i in range(1, k97):
        faces = np.append(faces, [[k4 + 2, k7 + i + 2, k7+i+1]], axis=0)
    faces = np.append(faces, [[k4 + 2, k7 + 1, k9]], axis=0)

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(1, qq[0] + 1):
        faces = np.append(faces, [[k4 + 1, k9 + i, k9 + i + 1]], axis=0)
    for i in range(1, len(rx2) - qq[-1]):
        faces = np.append(faces, [[k4 + 2, k10 + i, k10 + i + 1]], axis=0)
    # Agrega la cabeza
    faces = np.append(faces, [[k4 + 1, k4 + 2, k9 + 1]], axis=0)

#Creo el keyway en la parte de arriba
    # meto el primer punto del pe a parte porque casi siempre su #de triangulos no es igual al de los demás
    k02 = len(vertices) - 1
    vertices = np.append(vertices, [[rx2[qq[0]], ry2[qq[0]], anchoeng]], axis=0)
    k12 = len(vertices) - 1
    # Para los puntos del pe empiezo en +1 para que no se haga un desmadre ya que el primer grupo suele tener menos separaciones
    for i in range(qq[0] + 1, qq[-1] + 1):
        vertices = np.append(vertices, [[rx2[i], ry2[i], anchoeng]], axis=0)
    k22 = len(vertices) - 1

    # ppi1 es el primer punto del pi en el que se empiezan a generar las separaciones uniformes (o sea todos los triangulos verdes menos el primer y ultimo grupo)
    if len(qk) <= 1:
        ppi1 = 0
    else:
        ppi1 = len(qk) - 1
    # Notese que empieza en ppi y este puede ser 0 cuando el primer grupo es sólo una línea
    for i in range(qp[0] + ppi1, qp[-1] + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], anchoeng]], axis=0)
    k32 = len(vertices) - 1

    # Agrego el primer grupo de triángulos verdes
    for i in range(qp[0], qp[0] + ppi1 + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], anchoeng]], axis=0)
    k42 = len(vertices) - 1

    # Agrego los triángulos sempiternos (laterales magenta con rojo)
    if ppi1 == 0:
        vu = k22
    else:
        vu = k32
    vertices = np.append(vertices, [[px, py, anchoeng]], axis=0)
    vertices = np.append(vertices, [[px, -py, anchoeng]], axis=0)
    k52 = len(vertices) - 1
    # añade los triángulos magenta restantes
    vertices = np.append(vertices, [[prx, pry, anchoeng]], axis=0)
    k62 = len(vertices) - 1

    for i in range(0, qp[0] + 1):
        if ry1[i] > py:
            vertices = np.append(vertices, [[rx1[i], ry1[i], anchoeng]], axis=0)
    k72 = len(vertices) - 1
    k762 = k72 - k62
    cte = 0
    if k762 != 1:
        cte = 1

    # Añade los triangulos magenta inferiores restantes
    vertices = np.append(vertices, [[prx, -pry, anchoeng]], axis=0)
    k82 = len(vertices) - 1
    for i in range(qp[-1], len(rx1)):
        if ry1[i] < -py:
            vertices = np.append(vertices, [[rx1[i], ry1[i], anchoeng]], axis=0)
    k92 = len(vertices) - 1
    k972 = k92 - k82

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(0, qq[0] + 1):
        vertices = np.append(vertices, [[rx2[i], ry2[i], anchoeng]], axis=0)
    k102 = len(vertices) - 1

    for i in range(qq[-1], len(rx2)):
        vertices = np.append(vertices, [[rx2[i], ry2[i], anchoeng]], axis=0)
    k112 = len(vertices) - 1

    #AQUI ESTA EL ERROR
    u = 1
    #Puntos del pi sobrantes
    residuo = (k32-k22-1)%k
    iterac = k32-k22-residuo
    #Crea los tríangulos verdes
    for i in range(1,  z):
        for j in range(u + 0, u + k + 0):
            if j <(iterac):
                faces = np.append(faces, [[k12 + i, k22 + j + 1, k22 + j]], axis=0)
                p = 0
            if j<(k32-k22):
                v=u
            else:
                break
        u = u + k
    #En caso de que haya triangulos del pi que sobren ya que en el ultimo punto usado del pe sus triángulos estan incompletos, se ejecuta lo siguiente
    # para conectarlos con el punto del pe anterior
    x = qq[-1] - qq[0]
    if residuo!=0:
        for i in range(0, residuo):
            faces = np.append(faces, [[k12 + x, k22 + i + v+1, k22 + i+v]], axis=0)
            p=0

    # Creo el primer grupo de triángulos verdes
    for i in range(qq[0], qq[0] + 1):
        for j in range(1, 1 + ppi1):
            faces = np.append(faces, [[k02 + i, k32 + j + 1, k32 + j]], axis=0)
            p = 0
            if ppi1 == 0:
                print("Maldita sea Jimbo, el comunismo no funciona !")

    # Creo los complementos de los triangulos verdes para cerrar
    for i in range(1, qq[-1] - qq[0] + 1):
        faces = np.append(faces, [[k02 + i, k02 + i + 1, k22 + (i - 1) * (k) + 1]], axis=0)
        p=i


    # Creo los triángulos sempiternos (laterales magenta con rojo)
    faces = np.append(faces, [[k02 + 1, vu + 1, k42 + 1]], axis=0)
    #PISTA 1
    faces = np.append(faces, [[k22, k42 + 2,k32]], axis=0)

    #Creo los triangulos magenta positivos
    for i in range(1, k72 - k62 + cte):
        faces = np.append(faces, [[k42 + 1, k52 + i + 1,k52 + i]], axis=0)
        p=0

    # creo el sello para la esquina izquierda del cuadro
    if k72 - k62 == 1:
        faces = np.append(faces, [[k52 + 1, k42 + 1, vu + 1]], axis=0)
        p=0

    # Añade los triangulos magenta inferiores restantes
    for i in range(1, k972):
        faces = np.append(faces, [[k42 + 2, k72+i+2,k72+i+1]], axis=0)
        p=0
    faces = np.append(faces, [[k42 + 2, k72 + 1, k92]], axis=0)

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(1, qq[0] + 1):
        faces = np.append(faces, [[k42 + 1, k92 + i, k92 + i + 1]], axis=0)
        p=0
    for i in range(1, len(rx2) - qq[-1]):
        faces = np.append(faces, [[k42 + 2, k102 + i, k102 + i + 1]], axis=0)
        p=0
    # Agrega la cabeza
    faces = np.append(faces, [[k42 + 1, k42 + 2, k92 + 1]], axis=0)

#Creo las caras laterales
    # triangulos magenta positivos
    for i in range(0, k7 - k6):
        faces = np.append(faces, [[k6+i, k62+i,k6+i+1]], axis=0)
        faces = np.append(faces, [[k6+i+1, k62+i, k62+i+1]], axis=0)
        p = 0
    kn = k9 - k8
    # triangulos magenta negativos
    if kn > 0:
        for i in range(1, kn):
            faces = np.append(faces, [[k8+i,k82+i,k8+i+1]], axis=0)
            faces = np.append(faces, [[k8+i+1,k82+i,k82+i+1]], axis=0)
            p = 0
    faces = np.append(faces, [[k8,k9,k82]], axis=0)
    faces = np.append(faces, [[k9,k92,k82]], axis=0)
    # casi todos los conectes del pe con casi todos los del pi
    for i in range(1, k3 - k2):
        faces = np.append(faces, [[k2+i, k22+i, k2+i+1]], axis=0)
        faces = np.append(faces, [[k2+i+1, k22+i, k22+i+1]], axis=0)
        p = 0
    # conectes del pi para el primer vertice del pe
    for i in range(1, k4 - k3):
        faces = np.append(faces, [[k3+i, k32+i, k3+i+1]], axis=0)
        faces = np.append(faces, [[k3+i+1, k32+i, k32+i+1]], axis=0)
        p = 0
    #crea las caras laterales para cerrar el keyhole
    #Izquierda
    faces = np.append(faces, [[k4+1,k42+1, k5+1]], axis=0)
    faces = np.append(faces, [[k5+1, k42+1,k52+1]], axis=0)

    #Derecha
    faces = np.append(faces, [[k52, k5,k82]], axis=0)
    faces = np.append(faces, [[ k5,k8,k82]], axis=0)

    #Frontal
    faces = np.append(faces, [[k5,k52, k4+1]], axis=0)
    faces = np.append(faces, [[k52, k42+1, k4+1]], axis=0)

    # Code for normals
    for i in range(numcarasorig, len(faces)):
        we = faces[i]
        v1 = vertices[we[0]]
        v2 = vertices[we[1]]
        v3 = vertices[we[2]]
        normal = normal_to_surface(v1, v2, v3)
        normals = np.append(normals, [normal], axis=0)

    return vertices, faces, normals


def helcuadrado(vertices,faces, normals, z,radiopoligono,rcircunscrito,a1,anchoeng,angroth):
    if radiopoligono!=-1:
        numcarasorig = len(faces)
        # region Esta region es para indicar la rotación de los vértices del pe en la cara superior causada por ser helicoidal
        # listota es para reunir los ángulos de los orígenes de los dientes (el primero es 0 el segundo es 2*pi/z* 1 y así)
        listota = []
        for i in range(0, z):
            listota.append(2 * np.pi / z * i)

        # le resto el absoluto de el angulo de rotacion de la hélice a todos los valores de la lista para simular que rotaron
        listota = [element - abs(angroth) for element in listota]
        # lo que hace es encontrar el valor mínimo en la lista (o sea mas cercano a cero para indicar que ese diente está casi
        # casi sobre el eje x)
        arot = min(listota, key=abs) * -angroth / abs(angroth)
        # endregion

        jx=[]
        jy=[]
        jx2=[]
        jy2=[]
        ix=[]
        iy=[]
        for i in range(0,z):
            jx.append((radiopoligono)*np.cos(i*2*np.pi/z))
            jy.append((radiopoligono)*np.sin(i*2*np.pi/z))
            jx2.append((radiopoligono)*np.cos(i*2*np.pi/z+arot))
            jy2.append((radiopoligono)*np.sin(i*2*np.pi/z+arot))
        for i in range(0,4):
            ix.append(rcircunscrito* np.cos(i*2*np.pi/4))
            iy.append(rcircunscrito* np.sin(i*2*np.pi/4))
        #vertmax es el número máximo de vértices que puede conectar un punto del hexagono con los del polígono exterior
        vertmax =int((z-1)/4)+2
        #cvgrandes indica qué puntos del hexágono conectan con la cantidad vertmax con vértices del polígono exterior
        cvgrandes = z-int(z/4)*4
        #patrones es una lista que indica qué vertices del hexágono conectan con la cantidad vertmax (si cvgrandes=0 entonces todos los vertices del hexágono tienen vertmax)
        patrones= [[0, 1, 2, 3], [0], [0, 2], [0, 1, 2], [0, 1, 2, 3]]

        #Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        #Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        #Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx[-1], jy[-1], a1]], axis=0)
        qj = len(vertices)-1
        for i in range(0,len(jx)):
            vertices = np.append(vertices, [[jx[i], jy[i], a1]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], a1]], axis=0)
        qi = len(vertices)-1
        for i in range(0,len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], a1]], axis=0)
        qx= len(vertices)-1
        a=0
        #Ciclo for que va entre todos los puntos del hexágono
        for i in range(0,4):
            if i in patrones[cvgrandes]:
                #Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi+1)+i, (qj+1)-1+a, (qi+1)-1+i]], axis=0)
                for k in range(0,vertmax-1):
                    #Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qi+1)+i, (qj+1)+a+k, (qj+1)-1+a+k]], axis=0)
                #Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                #Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a=a+vertmax-1
            else:
                #Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi+1)+i, (qj+1)-1+a, (qi+1)-1+i]], axis=0)
                for m in range(0,vertmax-2):
                    #Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qi+1)+i, (qj+1)+a+m, (qj+1)-1+a+m]], axis=0)
                a=a+vertmax-2

        #Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        #Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        #Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx2[-1], jy2[-1], anchoeng]], axis=0)
        qj2 = len(vertices)-1
        for i in range(0,len(jx)):
            vertices = np.append(vertices, [[jx2[i], jy2[i], anchoeng]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], anchoeng]], axis=0)
        qi2 = len(vertices)-1
        for i in range(0,len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], anchoeng]], axis=0)
        qx= len(vertices)-1
        a=0
        #Ciclo for que va entre todos los puntos del hexágono
        for i in range(0,4):
            if i in patrones[cvgrandes]:
                #Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi2+1)+i, (qi2+1)-1+i, (qj2+1)-1+a]], axis=0)
                for k in range(0,vertmax-1):
                    #Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qi2+1)+i, (qj2+1)-1+a+k, (qj2+1)+a+k]], axis=0)
                #Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                #Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a=a+vertmax-1
            else:
                #Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi2+1)+i, (qi2+1)-1+i, (qj2+1)-1+a]], axis=0)
                for m in range(0,vertmax-2):
                    #Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qi2+1)+i, (qj2+1)-1+a+m, (qj2+1)+a+m]], axis=0)
                a=a+vertmax-2
        #Creo las caras laterales para el hexágono
        pisep=np.linspace(0,2*np.pi,5)
        xi = rcircunscrito*np.cos(pisep)
        yi = rcircunscrito*np.sin(pisep)

        qp=len(vertices)-1
        for i in range(0,len(pisep)):
            vertices = np.append(vertices, [[xi[i], yi[i], a1]], axis=0)
        qp1=len(vertices)-1
        for i in range(0,len(pisep)):
            vertices = np.append(vertices, [[xi[i], yi[i], anchoeng]], axis=0)
        for i in range(0,len(pisep)):
            faces= np.append(faces, [[qp+i+1, qp+i, qp1+i]], axis=0)
            faces = np.append(faces, [[qp+i+1, qp1+i, qp1+i+1]], axis=0)

        # Code for normals
        for i in range(numcarasorig, len(faces)):
            we = faces[i]
            v1 = vertices[we[0]]
            v2 = vertices[we[1]]
            v3 = vertices[we[2]]
            normal = normal_to_surface(v1, v2, v3)
            normals = np.append(normals, [normal], axis=0)
    return vertices, faces, normals

def helhexagono(vertices,faces, normals, z,radiopoligono,rcircunscrito,a1,anchoeng,angroth):

    if radiopoligono!=-1:
        numcarasorig = len(faces)
        # region Esta region es para indicar la rotación de los vértices del pe en la cara superior causada por ser helicoidal
        # listota es para reunir los ángulos de los orígenes de los dientes (el primero es 0 el segundo es 2*pi/z* 1 y así)
        listota = []
        for i in range(0, z):
            listota.append(2 * np.pi / z * i)

        # le resto el absoluto de el angulo de rotacion de la hélice a todos los valores de la lista para simular que rotaron
        listota = [element - abs(angroth) for element in listota]
        # lo que hace es encontrar el valor mínimo en la lista (o sea mas cercano a cero para indicar que ese diente está casi
        # casi sobre el eje x)
        arot = min(listota, key=abs) * -angroth / abs(angroth)
        # endregion

        jx=[]
        jy=[]
        jx2=[]
        jy2=[]
        ix=[]
        iy=[]
        for i in range(0,z):
            jx.append((radiopoligono)*np.cos(i*2*np.pi/z))
            jy.append((radiopoligono)*np.sin(i*2*np.pi/z))
            jx2.append((radiopoligono)*np.cos(i*2*np.pi/z+arot))
            jy2.append((radiopoligono)*np.sin(i*2*np.pi/z+arot))
        for i in range(0,6):
            ix.append(rcircunscrito* np.cos(i*2*np.pi/6))
            iy.append(rcircunscrito* np.sin(i*2*np.pi/6))
        #vertmax es el número máximo de vértices que puede conectar un punto del hexagono con los del polígono exterior
        vertmax =int((z-1)/6)+2
        #cvgrandes indica qué puntos del hexágono conectan con la cantidad vertmax con vértices del polígono exterior
        cvgrandes = z-int(z/6)*6
        #patrones es una lista que indica qué vertices del hexágono conectan con la cantidad vertmax (si cvgrandes=0 entonces todos los vertices del hexágono tienen vertmax)
        patrones= [[0,1,2,3,4,5], [0], [0,3], [0,2,4], [0,2,3,5], [0,1,2,3,4]]

        #Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        #Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        #Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx[-1], jy[-1], a1]], axis=0)
        qj = len(vertices)-1
        for i in range(0,len(jx)):
            vertices = np.append(vertices, [[jx[i], jy[i], a1]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], a1]], axis=0)
        qi = len(vertices)-1
        for i in range(0,len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], a1]], axis=0)
        qx= len(vertices)-1
        a=0
        #Ciclo for que va entre todos los puntos del hexágono
        for i in range(0,6):
            if i in patrones[cvgrandes]:
                #Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi+1)+i, (qj+1)-1+a, (qi+1)-1+i]], axis=0)
                for k in range(0,vertmax-1):
                    #Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qi+1)+i, (qj+1)+a+k, (qj+1)-1+a+k]], axis=0)
                #Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                #Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a=a+vertmax-1
            else:
                #Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi+1)+i, (qj+1)-1+a, (qi+1)-1+i]], axis=0)
                for m in range(0,vertmax-2):
                    #Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qi+1)+i, (qj+1)+a+m, (qj+1)-1+a+m]], axis=0)
                a=a+vertmax-2

        #Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        #Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        #Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx2[-1], jy2[-1], anchoeng]], axis=0)
        qj2 = len(vertices)-1
        for i in range(0,len(jx)):
            vertices = np.append(vertices, [[jx2[i], jy2[i], anchoeng]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], anchoeng]], axis=0)
        qi2 = len(vertices)-1
        for i in range(0,len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], anchoeng]], axis=0)
        qx= len(vertices)-1
        a=0
        #Ciclo for que va entre todos los puntos del hexágono
        for i in range(0,6):
            if i in patrones[cvgrandes]:
                #Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi2+1)+i, (qi2+1)-1+i, (qj2+1)-1+a]], axis=0)
                for k in range(0,vertmax-1):
                    #Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qi2+1)+i, (qj2+1)-1+a+k, (qj2+1)+a+k]], axis=0)
                #Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                #Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a=a+vertmax-1
            else:
                #Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi2+1)+i, (qi2+1)-1+i, (qj2+1)-1+a]], axis=0)
                for m in range(0,vertmax-2):
                    #Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qi2+1)+i, (qj2+1)-1+a+m, (qj2+1)+a+m]], axis=0)
                a=a+vertmax-2
        #Creo las caras laterales para el hexágono
        pisep=np.linspace(0,2*np.pi,7)
        xi = rcircunscrito*np.cos(pisep)
        yi = rcircunscrito*np.sin(pisep)

        qp=len(vertices)-1
        for i in range(0,len(pisep)):
            vertices = np.append(vertices, [[xi[i], yi[i], a1]], axis=0)
        qp1=len(vertices)-1
        for i in range(0,len(pisep)):
            vertices = np.append(vertices, [[xi[i], yi[i], anchoeng]], axis=0)
        for i in range(0,len(pisep)):
            faces= np.append(faces, [[qp+i+1, qp+i, qp1+i]], axis=0)
            faces = np.append(faces, [[qp+i+1, qp1+i, qp1+i+1]], axis=0)

        # Code for normals
        for i in range(numcarasorig, len(faces)):
            we = faces[i]
            v1 = vertices[we[0]]
            v2 = vertices[we[1]]
            v3 = vertices[we[2]]
            normal = normal_to_surface(v1, v2, v3)
            normals = np.append(normals, [normal], axis=0)
    return vertices, faces, normals

def helkeyway(vertices,faces, normals, z,radiopoligono,anchoeng,r1,alturaC,anchoC,angroth):
    numcarasorig = len(faces)

    #region Esta region es para indicar la rotación de los vértices del pe en la cara superior causada por ser helicoidal
    #listota es para reunir los ángulos de los orígenes de los dientes (el primero es 0 el segundo es 2*pi/z* 1 y así)
    listota = []
    for i in range(0,z):
        listota.append(2*np.pi/z * i)

    # le resto el absoluto de el angulo de rotacion de la hélice a todos los valores de la lista para simular que rotaron
    listota = [element-abs(angroth) for element in listota]
    #lo que hace es encontrar el valor mínimo en la lista (o sea mas cercano a cero para indicar que ese diente está casi
    # casi sobre el eje x)
    arot = min(listota, key=abs) * -angroth/abs(angroth)
    #endregion

    # Crea el keyhole
    # l1 es el poligono interior y l2 el exterior
    l2 = z
    if z <= 10:
        l1 = 20 * z
    elif 10 < z <= 20:
        l1 = 10 * z
    elif 20 < z <= 50:
        l1 = 5 * z
    elif 50 < z <= 100:
        l1 = 2 * z
    else:
        l1 = z

    # caracteristicas del keyway
    altura = alturaC
    ancho = anchoC

    # radios ocasionados por el keyway
    #r1 = 10
    r2 = radiopoligono

    # puntos del keyway
    px = r1 + altura
    py = ancho / 2
    prx = np.sqrt((r1 ** 2) - ((ancho / 2) ** 2))
    pry = ancho / 2

    pisep = np.linspace(0, 2 * np.pi, l1 + 1)
    pisep2 = np.linspace(0, 2 * np.pi, l2 + 1)

    rx1 = r1 * np.cos(pisep)
    ry1 = r1 * np.sin(pisep)

    rx2 = r2 * np.cos(pisep2)
    ry2 = r2 * np.sin(pisep2)

    rx22 = r2 * np.cos(pisep2+arot)
    ry22 = r2 * np.sin(pisep2+arot)

    # Encuentro el primer vert del pe que está fuera del cuadro
    for i in range(1, len(rx2)):
        a = rx2[i]
        b = ry2[i]
        if a > 0 and b > ancho / 2:
            item = i
            break
    for i in range(1, len(rx2)):
        a = rx2[i]
        b = ry2[i]
        if a > 0 and b < 0 and b >= -ancho / 2:
            item2 = i - 1
            break
    # Encuentro los vertices del pi que están fuera del cuadro
    for i in range(1, len(rx1)):
        a = rx1[i]
        b = ry1[i]
        if a > 0 and b > ancho / 2:
            item3 = i
            break
    for i in range(1, len(rx1)):
        a = rx1[i]
        b = ry1[i]
        if a > 0 and b < -ancho / 2:
            item4 = i

    k = int(l1 / l2)
    a = int(k / 2)
    # qq es para el polígono exterior y qp para el interior
    qq = []
    qp = []
    qk = []
    qk2 = []
    tr = True
    u = a
    for i in range(item, item2 + 1):
        for j in range(u, u + k + 1):
            if j >= item3 and j <= item4:
                qq.append(i)
                qp.append(j)
                # qk es una lista para encontrar el primer vertice del polígono exterior que conecta con uno del interior ya que item2 no necesariamente sera el 1ero
                # con len(qk) puedo saber cuantas lineas tiene: si es 0 o 1 tiene una linea, para los demás numeros tiene la cantidad que indiquen
                if tr:
                    qk.append(i)
                    qk2.append(j)
        tr = False
        u = u + k

    # meto el primer punto del pe a parte porque casi siempre su #de triangulos no es igual al de los demás
    k0 = len(vertices) - 1
    vertices = np.append(vertices, [[rx2[qq[0]], ry2[qq[0]], 0]], axis=0)
    k1 = len(vertices) - 1
    # Para los puntos del pe empiezo en +1 para que no se haga un desmadre ya que el primer grupo suele tener menos separaciones
    for i in range(qq[0] + 1, qq[-1] + 1):
        vertices = np.append(vertices, [[rx2[i], ry2[i], 0]], axis=0)
    k2 = len(vertices) - 1

    # ppi1 es el primer punto del pi en el que se empiezan a generar las separaciones uniformes (o sea todos los triangulos verdes menos el primer y ultimo grupo)
    if len(qk) <= 1:
        ppi1 = 0
    else:
        ppi1 = len(qk) - 1
    # Notese que empieza en ppi y este puede ser 0 cuando el primer grupo es sólo una línea
    for i in range(qp[0] + ppi1, qp[-1] + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k3 = len(vertices) - 1

    # Agrego el primer grupo de triángulos verdes
    for i in range(qp[0], qp[0] + ppi1 + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k4 = len(vertices) - 1

    # Agrego los triángulos sempiternos (laterales magenta con rojo)
    if ppi1 == 0:
        vu = k2
    else:
        vu = k3
    vertices = np.append(vertices, [[px, py, 0]], axis=0)
    vertices = np.append(vertices, [[px, -py, 0]], axis=0)
    k5 = len(vertices) - 1
    # añade los triángulos magenta restantes
    vertices = np.append(vertices, [[prx, pry, 0]], axis=0)
    k6 = len(vertices) - 1

    for i in range(0, qp[0] + 1):
        if ry1[i] > py:
            vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k7 = len(vertices) - 1
    k76 = k7 - k6
    cte = 0
    if k76 != 1:
        cte = 1

    # Añade los triangulos magenta inferiores restantes
    vertices = np.append(vertices, [[prx, -pry, 0]], axis=0)
    k8 = len(vertices) - 1
    for i in range(qp[-1], len(rx1)):
        if ry1[i] < -py:
            vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k9 = len(vertices) - 1
    k97 = k9 - k8

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(0, qq[0] + 1):
        vertices = np.append(vertices, [[rx2[i], ry2[i], 0]], axis=0)
    k10 = len(vertices) - 1

    for i in range(qq[-1], len(rx2)):
        vertices = np.append(vertices, [[rx2[i], ry2[i], 0]], axis=0)
    k11 = len(vertices) - 1

    u = 1
    # Puntos del pi sobrantes
    residuo = (k3 - k2 - 1) % k
    iterac = k3 - k2 - residuo

#Caras inferiores
    # Crea los tríangulos verdes
    for i in range(1, z):
        for j in range(u + 0, u + k + 0):
            if j < (iterac):
                faces = np.append(faces, [[k1 + i, k2 + j, k2 + j +1]], axis=0)
                p = 0
            if j < (k3 - k2):
                v = u
            else:
                break
        u = u + k
    # En caso de que haya triangulos del pi que sobren ya que en el ultimo punto usado del pe sus triángulos estan incompletos, se ejecuta lo siguiente
    # para conectarlos con el punto del pe anterior
    x = qq[-1] - qq[0]
    if residuo != 0:
        for i in range(0, residuo):
            faces = np.append(faces, [[k1 + x, k2 + i + v, k2 + i + v+1]], axis=0)

    # Creo el primer grupo de triángulos verdes
    for i in range(qq[0], qq[0] + 1):
        for j in range(1, 1 + ppi1):
            faces = np.append(faces, [[k0 + i,  k3 + j, k3 + j + 1]], axis=0)
            p = 0
            if ppi1 == 0:
                print("Maldita sea Jimbo, el comunismo no funciona !")

    # Creo los complementos de los triangulos verdes para cerrar
    for i in range(1, qq[-1] - qq[0] + 1):
        faces = np.append(faces, [[k0 + i,  k2 + (i - 1) * (k) + 1, k0 + i + 1]], axis=0)

    # Creo los triángulos sempiternos (laterales magenta con rojo)
    faces = np.append(faces, [[k0 + 1, k4 + 1, vu + 1]], axis=0)
    faces = np.append(faces, [[k2, k3, k4 + 2]], axis=0)

    #Agrego los triángulos magenta positivos
    for i in range(1, k7 - k6 + cte):
        faces = np.append(faces, [[k4 + 1,  k5 + i, k5 + i + 1]], axis=0)

    # creo el sello para la esquina izquierda del cuadro
    if k7 - k6 == 1:
        faces = np.append(faces, [[k5 + 1, vu + 1, k4 + 1]], axis=0)

    # Añade los triangulos magenta inferiores restantes
    for i in range(1, k97):
        faces = np.append(faces, [[k4 + 2, k7+i+1, k7 + i + 2]], axis=0)
    faces = np.append(faces, [[k4 + 2, k9, k7 + 1]], axis=0)

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(1, qq[0] + 1):
        faces = np.append(faces, [[k4 + 1, k9 + i+1, k9 + i]], axis=0)
    for i in range(1, len(rx2) - qq[-1]):
        faces = np.append(faces, [[k4 + 2, k10 + i +1, k10 + i]], axis=0)
    # Agrega la cabeza
    faces = np.append(faces, [[k4 + 1, k9 + 1, k4 + 2]], axis=0)

#Creo el keyway en la parte de arriba
    # meto el primer punto del pe a parte porque casi siempre su #de triangulos no es igual al de los demás
    k02 = len(vertices) - 1
    vertices = np.append(vertices, [[rx22[qq[0]], ry22[qq[0]], anchoeng]], axis=0)
    k12 = len(vertices) - 1
    # Para los puntos del pe empiezo en +1 para que no se haga un desmadre ya que el primer grupo suele tener menos separaciones
    for i in range(qq[0] + 1, qq[-1] + 1):
        vertices = np.append(vertices, [[rx22[i], ry22[i], anchoeng]], axis=0)
    k22 = len(vertices) - 1

    # ppi1 es el primer punto del pi en el que se empiezan a generar las separaciones uniformes (o sea todos los triangulos verdes menos el primer y ultimo grupo)
    if len(qk) <= 1:
        ppi1 = 0
    else:
        ppi1 = len(qk) - 1
    # Notese que empieza en ppi y este puede ser 0 cuando el primer grupo es sólo una línea
    for i in range(qp[0] + ppi1, qp[-1] + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], anchoeng]], axis=0)
    k32 = len(vertices) - 1

    # Agrego el primer grupo de triángulos verdes
    for i in range(qp[0], qp[0] + ppi1 + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], anchoeng]], axis=0)
    k42 = len(vertices) - 1

    # Agrego los triángulos sempiternos (laterales magenta con rojo)
    if ppi1 == 0:
        vu = k22
    else:
        vu = k32
    vertices = np.append(vertices, [[px, py, anchoeng]], axis=0)
    vertices = np.append(vertices, [[px, -py, anchoeng]], axis=0)
    k52 = len(vertices) - 1
    # añade los triángulos magenta restantes
    vertices = np.append(vertices, [[prx, pry, anchoeng]], axis=0)
    k62 = len(vertices) - 1

    for i in range(0, qp[0] + 1):
        if ry1[i] > py:
            vertices = np.append(vertices, [[rx1[i], ry1[i], anchoeng]], axis=0)
    k72 = len(vertices) - 1
    k762 = k72 - k62
    cte = 0
    if k762 != 1:
        cte = 1

    # Añade los triangulos magenta inferiores restantes
    vertices = np.append(vertices, [[prx, -pry, anchoeng]], axis=0)
    k82 = len(vertices) - 1
    for i in range(qp[-1], len(rx1)):
        if ry1[i] < -py:
            vertices = np.append(vertices, [[rx1[i], ry1[i], anchoeng]], axis=0)
    k92 = len(vertices) - 1
    k972 = k92 - k82

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(0, qq[0] + 1):
        vertices = np.append(vertices, [[rx22[i], ry22[i], anchoeng]], axis=0)
    k102 = len(vertices) - 1

    for i in range(qq[-1], len(rx22)):
        vertices = np.append(vertices, [[rx22[i], ry22[i], anchoeng]], axis=0)
    k112 = len(vertices) - 1

    #AQUI ESTA EL ERROR
    u = 1
    #Puntos del pi sobrantes
    residuo = (k32-k22-1)%k
    iterac = k32-k22-residuo
    #Crea los tríangulos verdes
    for i in range(1,  z):
        for j in range(u + 0, u + k + 0):
            if j <(iterac):
                faces = np.append(faces, [[k12 + i, k22 + j + 1, k22 + j]], axis=0)
                p = 0
            if j<(k32-k22):
                v=u
            else:
                break
        u = u + k
    #En caso de que haya triangulos del pi que sobren ya que en el ultimo punto usado del pe sus triángulos estan incompletos, se ejecuta lo siguiente
    # para conectarlos con el punto del pe anterior
    x = qq[-1] - qq[0]
    if residuo!=0:
        for i in range(0, residuo):
            faces = np.append(faces, [[k12 + x, k22 + i + v+1, k22 + i+v]], axis=0)
            p=0

    # Creo el primer grupo de triángulos verdes
    for i in range(qq[0], qq[0] + 1):
        for j in range(1, 1 + ppi1):
            faces = np.append(faces, [[k02 + i, k32 + j + 1, k32 + j]], axis=0)
            p = 0
            if ppi1 == 0:
                print("Maldita sea Jimbo, el comunismo no funciona !")

    # Creo los complementos de los triangulos verdes para cerrar
    for i in range(1, qq[-1] - qq[0] + 1):
        faces = np.append(faces, [[k02 + i, k02 + i + 1, k22 + (i - 1) * (k) + 1]], axis=0)
        p=i

    # Creo los triángulos sempiternos (laterales magenta con rojo)
    faces = np.append(faces, [[k02 + 1, vu + 1, k42 + 1]], axis=0)
    faces = np.append(faces, [[k22, k42 + 2,k32]], axis=0)

    #Creo los triangulos magenta positivos
    for i in range(1, k72 - k62 + cte):
        faces = np.append(faces, [[k42 + 1, k52 + i + 1,k52 + i]], axis=0)
        p=0

    # creo el sello para la esquina izquierda del cuadro
    if k72 - k62 == 1:
        faces = np.append(faces, [[k52 + 1, k42 + 1, vu + 1]], axis=0)
        p=0

    # Añade los triangulos magenta inferiores restantes
    for i in range(1, k972):
        faces = np.append(faces, [[k42 + 2, k72+i+2,k72+i+1]], axis=0)
        p=0
    faces = np.append(faces, [[k42 + 2, k72 + 1, k92]], axis=0)

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(1, qq[0] + 1):
        faces = np.append(faces, [[k42 + 1, k92 + i, k92 + i + 1]], axis=0)
        p=0
    for i in range(1, len(rx2) - qq[-1]):
        faces = np.append(faces, [[k42 + 2, k102 + i, k102 + i + 1]], axis=0)
        p=0
    # Agrega la cabeza
    faces = np.append(faces, [[k42 + 1, k42 + 2, k92 + 1]], axis=0)

#Creo las caras laterales
    # triangulos magenta positivos
    for i in range(0, k7 - k6):
        faces = np.append(faces, [[k6+i, k62+i,k6+i+1]], axis=0)
        faces = np.append(faces, [[k6+i+1, k62+i, k62+i+1]], axis=0)
        p = 0
    kn = k9 - k8
    # triangulos magenta negativos
    if kn > 0:
        for i in range(1, kn):
            faces = np.append(faces, [[k8+i,k82+i,k8+i+1]], axis=0)
            faces = np.append(faces, [[k8+i+1,k82+i,k82+i+1]], axis=0)
            p = 0
    faces = np.append(faces, [[k8,k9,k82]], axis=0)
    faces = np.append(faces, [[k9,k92,k82]], axis=0)
    # casi todos los conectes del pe con casi todos los del pi
    for i in range(1, k3 - k2):
        faces = np.append(faces, [[k2+i, k22+i, k2+i+1]], axis=0)
        faces = np.append(faces, [[k2+i+1, k22+i, k22+i+1]], axis=0)
        p = 0
    # conectes del pi para el primer vertice del pe
    for i in range(1, k4 - k3):
        faces = np.append(faces, [[k3+i, k32+i, k3+i+1]], axis=0)
        faces = np.append(faces, [[k3+i+1, k32+i, k32+i+1]], axis=0)
        p = 0
    #crea las caras laterales para cerrar el keyhole
    #Izquierda
    faces = np.append(faces, [[k4+1,k42+1, k5+1]], axis=0)
    faces = np.append(faces, [[k5+1, k42+1,k52+1]], axis=0)

    #Derecha
    faces = np.append(faces, [[k52, k5,k82]], axis=0)
    faces = np.append(faces, [[ k5,k8,k82]], axis=0)

    #Frontal
    faces = np.append(faces, [[k5,k52, k4+1]], axis=0)
    faces = np.append(faces, [[k52, k42+1, k4+1]], axis=0)

    # Code for normals
    for i in range(numcarasorig, len(faces)):
        we = faces[i]
        v1 = vertices[we[0]]
        v2 = vertices[we[1]]
        v3 = vertices[we[2]]
        normal = normal_to_surface(v1, v2, v3)
        normals = np.append(normals, [normal], axis=0)

    return vertices, faces, normals

def helcirculo(vertices,faces,normals,z,radiopoligono,a1,anchoeng,rcircunscrito,angroth):
    #el circle_points es porque z es 10, pero esto debe quedar de forma que x/z siempre de un entero
    if radiopoligono != -1 and rcircunscrito <= radiopoligono-.15:
        listota = []
        numcarasorig = len(faces)
        for i in range(0, z):
            listota.append(2 * np.pi / z * i)

        # le resto el absoluto de el angulo de rotacion de la hélice a todos los valores de la lista para simular que rotaron
        listota = [element - abs(angroth) for element in listota]
        # lo que hace es encontrar el valor mínimo en la lista (o sea mas cercano a cero para indicar que ese diente está casi
        # casi sobre el eje x)
        arot = min(listota, key=abs) * -angroth / abs(angroth)


        #NOTA: SI LOS PUNTOS DE LA CARA SUPERIOR NO CONECTAN BIEN, PUEDES HACER NEGATIVO EL arot Y ESO DEBERIA ARREGLARLO


        #ATENCION !!!! TIENES QUE
        # MODFICAR CIRCLE POINTS PARA QUE SEA UNA VARABLE MAS EFCIENTE A LA HORA DE SELECCIONAR LA CANTIDAD DE
        # CARAS PARA EL CIRCULO
        circle_points = int(300/z)*z
        circle_points = 120
        circle_points_per_vertex = circle_points / z
        half_angle = (circle_points_per_vertex * 1 * np.pi / circle_points)

        # Para esta función estoy usando "mal" a1 porque esa se planeaba para que la cara del agujero no estuviera al mismo nivel que las caras inefirores (como otro agujero o que estuviera sumida)
        # a1 la voy a usar para decirle si son las caras de abajo o la de arriba
        def base_faces(vertices,faces,z,radiopoligono,rcircunscrito,a1,normals_downwards,arot):
            jx = []
            jy = []
            ix = []
            iy = []
            for i in range(0, z):
                jx.append((radiopoligono) * np.cos(i * 2 * np.pi / z + arot))
                jy.append((radiopoligono) * np.sin(i * 2 * np.pi / z + arot))
            # agrego otra vez cuando el angulo=0 para que en el loop al dar la vuelta llegue al ultimo vertice del circulo (que es el inicial porque 2pi=0)
            jx.append((radiopoligono) * np.cos(0 * 2 * np.pi / z + arot))
            jy.append((radiopoligono) * np.sin(0 * 2 * np.pi / z + arot))
            for i in range(0, circle_points):
                ix.append(rcircunscrito * np.cos((i * 2 * np.pi / circle_points) - half_angle))
                iy.append(rcircunscrito * np.sin((i * 2 * np.pi / circle_points) - half_angle))
            # agrego otra vez cuando el angulo=0 para que en el loop al dar la vuelta llegue al ultimo vertice del circulo (que es el inicial porque 2pi=0)
            ix.append(rcircunscrito * np.cos((0 * 2 * np.pi / circle_points) - half_angle))
            iy.append(rcircunscrito * np.sin((0 * 2 * np.pi / circle_points) - half_angle))

            carasporvt = int(circle_points / z)

            qj = len(vertices) - 0
            for i in range(0, len(jx)):
                vertices = np.append(vertices, [[jx[i], jy[i], a1]], axis=0)
            # ***
            vertices = np.append(vertices, [[ix[-1], iy[-1], a1]], axis=0)
            qi = len(vertices) - 0
            for i in range(0, len(ix)):
                vertices = np.append(vertices, [[ix[i], iy[i], a1]], axis=0)
            qx = len(vertices) - 1

            for i in range(0, z):
                for j in range(0, carasporvt +0):
                    x = int(i*circle_points_per_vertex + j)
                    if normals_downwards:
                        faces = np.append(faces, [[qj + i, qi + x, qi + 1 + x]], axis=0)
                    else:
                        faces = np.append(faces, [[qi + x, qj + i, qi + 1 + x]], axis=0)
                if normals_downwards:
                    faces = np.append(faces, [[qj + i + 1, qj + i, qi + 1 + x]], axis=0)
                else:
                    faces = np.append(faces, [[qj + i, qj + 1 + i, qi + 1 + x]], axis=0)
            return vertices, faces, qi

        def lateral_faces(faces, qi1, qi2, z):
            #   NOTA IMPORTANTE: CAMBIA EL PINCHE circle_points !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            for i in range(0, circle_points):
                faces = np.append(faces, [[qi1 + i, qi2 + i, qi1 + i + 1]], axis=0)
                faces = np.append(faces, [[qi1 + i + 1, qi2 + i, qi2 + i + 1]], axis=0)
            return faces

        vf = base_faces(vertices, faces, z, radiopoligono, rcircunscrito, a1, True, 0)
        vertices = vf[0]
        faces = vf[1]
        qi1 = vf[2]
        vf = base_faces(vertices, faces, z, radiopoligono, rcircunscrito, anchoeng, False, arot)
        vertices = vf[0]
        faces = vf[1]
        qi2 = vf[2]


        faces = lateral_faces(faces, qi1, qi2, z)
        # Code for normals
        for i in range(numcarasorig, len(faces)):
            we = faces[i]
            v1 = vertices[we[0]]
            v2 = vertices[we[1]]
            v3 = vertices[we[2]]
            normal = normal_to_surface(v1, v2, v3)
            normals = np.append(normals, [normal], axis=0)

    else:
        print("Radio del polígono circumscrito es muy grande, pruebe con un valor <" + str(radiopoligono))
    return vertices, faces, normals


def dblhelhexagono(vertices,faces,z,radiopoligono,rcircunscrito,a1,anchoeng,angroth):

    #Este bloque es para indicar la rotación de los vértices del pe en la cara superior causada por ser helicoidal
    if radiopoligono!=-1:
        ahx = 2*np.pi/6

        #Con el siguiente if resolví el bug2 de helhexagono
        if angroth>2*np.pi:
            intangroth = int(angroth/(2*np.pi))
            residuoangroth  =   angroth-intangroth*2*np.pi
            angroth = residuoangroth

        #La siguiente línea resuelve el bug1 de helexagono
        if angroth!=ahx:
            for i in range(0,z):
                a = 2*np.pi/z
                if i*a+angroth>2*np.pi:
                    break
            arot = angroth+(2*np.pi/z)*(i)-2*np.pi
        else:
            arot = angroth
    #Fin del bloque
        jx=[]
        jy=[]
        jx2=[]
        jy2=[]
        ix=[]
        iy=[]
        for i in range(0,z):
            jx.append((radiopoligono)*np.cos(i*2*np.pi/z))
            jy.append((radiopoligono)*np.sin(i*2*np.pi/z))
            jx2.append((radiopoligono)*np.cos(i*2*np.pi/z+arot))
            jy2.append((radiopoligono)*np.sin(i*2*np.pi/z+arot))
        for i in range(0,6):
            ix.append(rcircunscrito* np.cos(i*2*np.pi/6))
            iy.append(rcircunscrito* np.sin(i*2*np.pi/6))
        #vertmax es el número máximo de vértices que puede conectar un punto del hexagono con los del polígono exterior
        vertmax =int((z-1)/6)+2
        #cvgrandes indica qué puntos del hexágono conectan con la cantidad vertmax con vértices del polígono exterior
        cvgrandes = z-int(z/6)*6
        #patrones es una lista que indica qué vertices del hexágono conectan con la cantidad vertmax (si cvgrandes=0 entonces todos los vertices del hexágono tienen vertmax)
        patrones= [[0,1,2,3,4,5], [0], [0,3], [0,2,4], [0,2,3,5], [0,1,2,3,4]]



        '''
        #Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        #Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        #Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx[-1], jy[-1], a1]], axis=0)
        qj = len(vertices)-1
        for i in range(0,len(jx)):
            vertices = np.append(vertices, [[jx[i], jy[i], a1]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], a1]], axis=0)
        qi = len(vertices)-1
        for i in range(0,len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], a1]], axis=0)
        qx= len(vertices)-1
        a=0
        #Ciclo for que va entre todos los puntos del hexágono
        for i in range(0,6):
            if i in patrones[cvgrandes]:
                #Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi+1)+i, (qj+1)-1+a, (qi+1)-1+i]], axis=0)
                for k in range(0,vertmax-1):
                    #Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qi+1)+i, (qj+1)+a+k, (qj+1)-1+a+k]], axis=0)
                #Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                #Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a=a+vertmax-1
            else:
                #Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi+1)+i, (qj+1)-1+a, (qi+1)-1+i]], axis=0)
                for m in range(0,vertmax-2):
                    #Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qi+1)+i, (qj+1)+a+m, (qj+1)-1+a+m]], axis=0)
                a=a+vertmax-2
        '''

        #Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        #Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        #Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx2[-1], jy2[-1], a1]], axis=0)
        qj = len(vertices)-1
        for i in range(0,len(jx)):
            vertices = np.append(vertices, [[jx2[i], jy2[i], a1]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], a1]], axis=0)
        qi = len(vertices)-1
        for i in range(0,len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], a1]], axis=0)
        qx= len(vertices)-1
        a=0
        #Ciclo for que va entre todos los puntos del hexágono
        for i in range(0,6):
            if i in patrones[cvgrandes]:
                #Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi+1)+i, (qj+1)-1+a, (qi+1)-1+i]], axis=0)
                for k in range(0,vertmax-1):
                    #Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qi+1)+i, (qj+1)+a+k, (qj+1)-1+a+k]], axis=0)
                #Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                #Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a=a+vertmax-1
            else:
                #Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi+1)+i, (qj+1)-1+a, (qi+1)-1+i]], axis=0)
                for m in range(0,vertmax-2):
                    #Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qi+1)+i, (qj+1)+a+m, (qj+1)-1+a+m]], axis=0)
                a=a+vertmax-2


        #Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        #Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        #Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx2[-1], jy2[-1], anchoeng]], axis=0)
        qj2 = len(vertices)-1
        for i in range(0,len(jx)):
            vertices = np.append(vertices, [[jx2[i], jy2[i], anchoeng]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], anchoeng]], axis=0)
        qi2 = len(vertices)-1
        for i in range(0,len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], anchoeng]], axis=0)
        qx= len(vertices)-1
        a=0
        #Ciclo for que va entre todos los puntos del hexágono
        for i in range(0,6):
            if i in patrones[cvgrandes]:
                #Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi2+1)+i, (qi2+1)-1+i, (qj2+1)-1+a]], axis=0)
                for k in range(0,vertmax-1):
                    #Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qi2+1)+i, (qj2+1)-1+a+k, (qj2+1)+a+k]], axis=0)
                #Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                #Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a=a+vertmax-1
            else:
                #Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi2+1)+i, (qi2+1)-1+i, (qj2+1)-1+a]], axis=0)
                for m in range(0,vertmax-2):
                    #Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qi2+1)+i, (qj2+1)-1+a+m, (qj2+1)+a+m]], axis=0)
                a=a+vertmax-2


        #Creo las caras laterales para el hexágono
        pisep=np.linspace(0,2*np.pi,7)
        xi = rcircunscrito*np.cos(pisep)
        yi = rcircunscrito*np.sin(pisep)

        qp=len(vertices)-1
        for i in range(0,len(pisep)):
            vertices = np.append(vertices, [[xi[i], yi[i], a1]], axis=0)
        qp1=len(vertices)-1
        for i in range(0,len(pisep)):
            vertices = np.append(vertices, [[xi[i], yi[i], anchoeng]], axis=0)
        for i in range(0,len(pisep)):
            faces= np.append(faces, [[qp+i+1, qp+i, qp1+i]], axis=0)
            faces = np.append(faces, [[qp+i+1, qp1+i, qp1+i+1]], axis=0)

    return vertices, faces

def dblhelcuadrado(vertices,faces,z,radiopoligono,rcircunscrito,a1,anchoeng,angroth):

    #Este bloque es para indicar la rotación de los vértices del pe en la cara superior causada por ser helicoidal
    if radiopoligono!=-1:
        ahx = 2*np.pi/4

        #Con el siguiente if resolví el bug2 de helhexagono
        if angroth>2*np.pi:
            intangroth = int(angroth/(2*np.pi))
            residuoangroth  =   angroth-intangroth*2*np.pi
            angroth = residuoangroth

        #La siguiente línea resuelve el bug1 de helexagono
        if angroth!=ahx:
            for i in range(0,z):
                a = 2*np.pi/z
                if i*a+angroth>2*np.pi:
                    break
            arot = angroth+(2*np.pi/z)*(i)-2*np.pi
        else:
            arot = angroth
    #Fin del bloque
        jx=[]
        jy=[]
        jx2=[]
        jy2=[]
        ix=[]
        iy=[]
        for i in range(0,z):
            jx.append((radiopoligono)*np.cos(i*2*np.pi/z))
            jy.append((radiopoligono)*np.sin(i*2*np.pi/z))
            jx2.append((radiopoligono)*np.cos(i*2*np.pi/z+arot))
            jy2.append((radiopoligono)*np.sin(i*2*np.pi/z+arot))
        for i in range(0,4):
            ix.append(rcircunscrito* np.cos(i*2*np.pi/4))
            iy.append(rcircunscrito* np.sin(i*2*np.pi/4))
        #vertmax es el número máximo de vértices que puede conectar un punto del hexagono con los del polígono exterior
        vertmax =int((z-1)/4)+2
        #cvgrandes indica qué puntos del hexágono conectan con la cantidad vertmax con vértices del polígono exterior
        cvgrandes = z-int(z/4)*4
        #patrones es una lista que indica qué vertices del hexágono conectan con la cantidad vertmax (si cvgrandes=0 entonces todos los vertices del hexágono tienen vertmax)
        patrones= [[0, 1, 2, 3], [0], [0, 2], [0, 1, 2], [0, 1, 2, 3]]



        '''
        #Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        #Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        #Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx[-1], jy[-1], a1]], axis=0)
        qj = len(vertices)-1
        for i in range(0,len(jx)):
            vertices = np.append(vertices, [[jx[i], jy[i], a1]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], a1]], axis=0)
        qi = len(vertices)-1
        for i in range(0,len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], a1]], axis=0)
        qx= len(vertices)-1
        a=0
        #Ciclo for que va entre todos los puntos del hexágono
        for i in range(0,6):
            if i in patrones[cvgrandes]:
                #Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi+1)+i, (qj+1)-1+a, (qi+1)-1+i]], axis=0)
                for k in range(0,vertmax-1):
                    #Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qi+1)+i, (qj+1)+a+k, (qj+1)-1+a+k]], axis=0)
                #Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                #Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a=a+vertmax-1
            else:
                #Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi+1)+i, (qj+1)-1+a, (qi+1)-1+i]], axis=0)
                for m in range(0,vertmax-2):
                    #Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qi+1)+i, (qj+1)+a+m, (qj+1)-1+a+m]], axis=0)
                a=a+vertmax-2
        '''

        #Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        #Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        #Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx2[-1], jy2[-1], a1]], axis=0)
        qj = len(vertices)-1
        for i in range(0,len(jx)):
            vertices = np.append(vertices, [[jx2[i], jy2[i], a1]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], a1]], axis=0)
        qi = len(vertices)-1
        for i in range(0,len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], a1]], axis=0)
        qx= len(vertices)-1
        a=0
        #Ciclo for que va entre todos los puntos del hexágono
        for i in range(0,4):
            if i in patrones[cvgrandes]:
                #Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi+1)+i, (qj+1)-1+a, (qi+1)-1+i]], axis=0)
                for k in range(0,vertmax-1):
                    #Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qi+1)+i, (qj+1)+a+k, (qj+1)-1+a+k]], axis=0)
                #Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                #Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a=a+vertmax-1
            else:
                #Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi+1)+i, (qj+1)-1+a, (qi+1)-1+i]], axis=0)
                for m in range(0,vertmax-2):
                    #Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qi+1)+i, (qj+1)+a+m, (qj+1)-1+a+m]], axis=0)
                a=a+vertmax-2


        #Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        #Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        #Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx2[-1], jy2[-1], anchoeng]], axis=0)
        qj2 = len(vertices)-1
        for i in range(0,len(jx)):
            vertices = np.append(vertices, [[jx2[i], jy2[i], anchoeng]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], anchoeng]], axis=0)
        qi2 = len(vertices)-1
        for i in range(0,len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], anchoeng]], axis=0)
        qx= len(vertices)-1
        a=0
        #Ciclo for que va entre todos los puntos del hexágono
        for i in range(0,4):
            if i in patrones[cvgrandes]:
                #Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi2+1)+i, (qi2+1)-1+i, (qj2+1)-1+a]], axis=0)
                for k in range(0,vertmax-1):
                    #Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qi2+1)+i, (qj2+1)-1+a+k, (qj2+1)+a+k]], axis=0)
                #Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                #Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a=a+vertmax-1
            else:
                #Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi2+1)+i, (qi2+1)-1+i, (qj2+1)-1+a]], axis=0)
                for m in range(0,vertmax-2):
                    #Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qi2+1)+i, (qj2+1)-1+a+m, (qj2+1)+a+m]], axis=0)
                a=a+vertmax-2


        #Creo las caras laterales para el hexágono
        pisep=np.linspace(0,2*np.pi,5)
        xi = rcircunscrito*np.cos(pisep)
        yi = rcircunscrito*np.sin(pisep)

        qp=len(vertices)-1
        for i in range(0,len(pisep)):
            vertices = np.append(vertices, [[xi[i], yi[i], a1]], axis=0)
        qp1=len(vertices)-1
        for i in range(0,len(pisep)):
            vertices = np.append(vertices, [[xi[i], yi[i], anchoeng]], axis=0)
        for i in range(0,len(pisep)):
            faces= np.append(faces, [[qp+i+1, qp+i, qp1+i]], axis=0)
            faces = np.append(faces, [[qp+i+1, qp1+i, qp1+i+1]], axis=0)

    return vertices, faces

def dblhelkeyway(vertices,faces,z,radiopoligono,anchoeng,r1,alturaC,anchoC,angroth):

    #Este bloque es para indicar la rotación de los vértices del pe en la cara superior causada por ser helicoidal
    breaker=0
    for i in range(0,z):
        a = 2*np.pi/z
        if i*a+angroth>2*np.pi:
            breaker=1
            arotq = abs(angroth+(2*np.pi/z)*(i-1)-2*np.pi)
            arot2 = abs(angroth+(2*np.pi/z)*(i)-2*np.pi)
            if arotq>arot2:
                arot = angroth+(2*np.pi/z)*(i)-2*np.pi
            else:
                arot = angroth+(2*np.pi/z)*(i-1)-2*np.pi
            break
    if breaker==0:
        arot = angroth
    #Fin del bloque

    # Crea el keyhole
    # l1 es el poligono interior y l2 el exterior
    l2 = z
    if z <= 10:
        l1 = 20 * z
    elif 10 < z <= 20:
        l1 = 10 * z
    elif 20 < z <= 50:
        l1 = 5 * z
    elif 50 < z <= 100:
        l1 = 2 * z
    else:
        l1 = z

    # caracteristicas del keyway
    altura = alturaC
    ancho = anchoC

    # radios ocasionados por el keyway
    #r1 = 10
    r2 = radiopoligono

    # puntos del keyway
    px = r1 + altura
    py = ancho / 2
    prx = np.sqrt((r1 ** 2) - ((ancho / 2) ** 2))
    pry = ancho / 2

    pisep = np.linspace(0, 2 * np.pi, l1 + 1)
    pisep2 = np.linspace(0, 2 * np.pi, l2 + 1)

    rx1 = r1 * np.cos(pisep)
    ry1 = r1 * np.sin(pisep)

    rx2 = r2 * np.cos(pisep2)
    ry2 = r2 * np.sin(pisep2)

    rx22 = r2 * np.cos(pisep2+arot)
    ry22 = r2 * np.sin(pisep2+arot)

    # Encuentro el primer vert del pe que está fuera del cuadro
    for i in range(1, len(rx2)):
        a = rx2[i]
        b = ry2[i]
        if a > 0 and b > ancho / 2:
            item = i
            break
    for i in range(1, len(rx2)):
        a = rx2[i]
        b = ry2[i]
        if a > 0 and b < 0 and b >= -ancho / 2:
            item2 = i - 1
            break
    # Encuentro los vertices del pi que están fuera del cuadro
    for i in range(1, len(rx1)):
        a = rx1[i]
        b = ry1[i]
        if a > 0 and b > ancho / 2:
            item3 = i
            break
    for i in range(1, len(rx1)):
        a = rx1[i]
        b = ry1[i]
        if a > 0 and b < -ancho / 2:
            item4 = i

    k = int(l1 / l2)
    a = int(k / 2)
    # qq es para el polígono exterior y qp para el interior
    qq = []
    qp = []
    qk = []
    qk2 = []
    tr = True
    u = a
    for i in range(item, item2 + 1):
        for j in range(u, u + k + 1):
            if j >= item3 and j <= item4:
                qq.append(i)
                qp.append(j)
                # qk es una lista para encontrar el primer vertice del polígono exterior que conecta con uno del interior ya que item2 no necesariamente sera el 1ero
                # con len(qk) puedo saber cuantas lineas tiene: si es 0 o 1 tiene una linea, para los demás numeros tiene la cantidad que indiquen
                if tr:
                    qk.append(i)
                    qk2.append(j)
        tr = False
        u = u + k

    # meto el primer punto del pe a parte porque casi siempre su #de triangulos no es igual al de los demás
    k0 = len(vertices) - 1
    vertices = np.append(vertices, [[rx2[qq[0]], ry2[qq[0]], 0]], axis=0)
    k1 = len(vertices) - 1
    # Para los puntos del pe empiezo en +1 para que no se haga un desmadre ya que el primer grupo suele tener menos separaciones
    for i in range(qq[0] + 1, qq[-1] + 1):
        vertices = np.append(vertices, [[rx2[i], ry2[i], 0]], axis=0)
    k2 = len(vertices) - 1

    # ppi1 es el primer punto del pi en el que se empiezan a generar las separaciones uniformes (o sea todos los triangulos verdes menos el primer y ultimo grupo)
    if len(qk) <= 1:
        ppi1 = 0
    else:
        ppi1 = len(qk) - 1
    # Notese que empieza en ppi y este puede ser 0 cuando el primer grupo es sólo una línea
    for i in range(qp[0] + ppi1, qp[-1] + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k3 = len(vertices) - 1

    # Agrego el primer grupo de triángulos verdes
    for i in range(qp[0], qp[0] + ppi1 + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k4 = len(vertices) - 1

    # Agrego los triángulos sempiternos (laterales magenta con rojo)
    if ppi1 == 0:
        vu = k2
    else:
        vu = k3
    vertices = np.append(vertices, [[px, py, 0]], axis=0)
    vertices = np.append(vertices, [[px, -py, 0]], axis=0)
    k5 = len(vertices) - 1
    # añade los triángulos magenta restantes
    vertices = np.append(vertices, [[prx, pry, 0]], axis=0)
    k6 = len(vertices) - 1

    for i in range(0, qp[0] + 1):
        if ry1[i] > py:
            vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k7 = len(vertices) - 1
    k76 = k7 - k6
    cte = 0
    if k76 != 1:
        cte = 1

    # Añade los triangulos magenta inferiores restantes
    vertices = np.append(vertices, [[prx, -pry, 0]], axis=0)
    k8 = len(vertices) - 1
    for i in range(qp[-1], len(rx1)):
        if ry1[i] < -py:
            vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k9 = len(vertices) - 1
    k97 = k9 - k8

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(0, qq[0] + 1):
        vertices = np.append(vertices, [[rx2[i], ry2[i], 0]], axis=0)
    k10 = len(vertices) - 1

    for i in range(qq[-1], len(rx2)):
        vertices = np.append(vertices, [[rx2[i], ry2[i], 0]], axis=0)
    k11 = len(vertices) - 1

    u = 1
    # Puntos del pi sobrantes
    residuo = (k3 - k2 - 1) % k
    iterac = k3 - k2 - residuo


# Creo el keyway en la parte de abajo
    # meto el primer punto del pe a parte porque casi siempre su #de triangulos no es igual al de los demás

    k02 = len(vertices) - 1
    vertices = np.append(vertices, [[rx22[qq[0]], ry22[qq[0]], 0]], axis=0)
    k12 = len(vertices) - 1
    # Para los puntos del pe empiezo en +1 para que no se haga un desmadre ya que el primer grupo suele tener menos separaciones
    for i in range(qq[0] + 1, qq[-1] + 1):
        vertices = np.append(vertices, [[rx22[i], ry22[i], 0]], axis=0)
    k22 = len(vertices) - 1

    # ppi1 es el primer punto del pi en el que se empiezan a generar las separaciones uniformes (o sea todos los triangulos verdes menos el primer y ultimo grupo)
    if len(qk) <= 1:
        ppi1 = 0
    else:
        ppi1 = len(qk) - 1
    # Notese que empieza en ppi y este puede ser 0 cuando el primer grupo es sólo una línea
    for i in range(qp[0] + ppi1, qp[-1] + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k32 = len(vertices) - 1

    # Agrego el primer grupo de triángulos verdes
    for i in range(qp[0], qp[0] + ppi1 + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k42 = len(vertices) - 1

    # Agrego los triángulos sempiternos (laterales magenta con rojo)
    if ppi1 == 0:
        vu = k22
    else:
        vu = k32
    vertices = np.append(vertices, [[px, py, 0]], axis=0)
    vertices = np.append(vertices, [[px, -py, 0]], axis=0)
    k52 = len(vertices) - 1
    # añade los triángulos magenta restantes
    vertices = np.append(vertices, [[prx, pry, 0]], axis=0)
    k62 = len(vertices) - 1

    for i in range(0, qp[0] + 1):
        if ry1[i] > py:
            vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k72 = len(vertices) - 1
    k762 = k72 - k62
    cte = 0
    if k762 != 1:
        cte = 1

    # Añade los triangulos magenta inferiores restantes
    vertices = np.append(vertices, [[prx, -pry, 0]], axis=0)
    k82 = len(vertices) - 1
    for i in range(qp[-1], len(rx1)):
        if ry1[i] < -py:
            vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k92 = len(vertices) - 1
    k972 = k92 - k82

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(0, qq[0] + 1):
        vertices = np.append(vertices, [[rx22[i], ry22[i], 0]], axis=0)
    k102 = len(vertices) - 1

    for i in range(qq[-1], len(rx22)):
        vertices = np.append(vertices, [[rx22[i], ry22[i], 0]], axis=0)
    k112 = len(vertices) - 1

    # AQUI ESTA EL ERROR
    u = 1
    # Puntos del pi sobrantes
    residuo = (k32 - k22 - 1) % k
    iterac = k32 - k22 - residuo
    # Crea los tríangulos verdes
    for i in range(1, z):
        for j in range(u + 0, u + k + 0):
            if j < (iterac):
                faces = np.append(faces, [[k12 + i, k22 + j, k22 + j + 1]], axis=0)
                p = 0
            if j < (k32 - k22):
                v = u
            else:
                break
        u = u + k
    # En caso de que haya triangulos del pi que sobren ya que en el ultimo punto usado del pe sus triángulos estan incompletos, se ejecuta lo siguiente
    # para conectarlos con el punto del pe anterior
    x = qq[-1] - qq[0]
    if residuo != 0:
        for i in range(0, residuo):
            faces = np.append(faces, [[k12 + x, k22 + i + v, k22 + i + v + 1]], axis=0)
            p = 0

    # Creo el primer grupo de triángulos verdes
    for i in range(qq[0], qq[0] + 1):
        for j in range(1, 1 + ppi1):
            faces = np.append(faces, [[k02 + i, k32 + j, k32 + j + 1]], axis=0)
            p = 0
            if ppi1 == 0:
                print("Maldita sea Jimbo, el comunismo no funciona !")

    # Creo los complementos de los triangulos verdes para cerrar
    for i in range(1, qq[-1] - qq[0] + 1):
        faces = np.append(faces, [[k02 + i, k22 + (i - 1) * (k) + 1, k02 + i + 1]], axis=0)
        p = i

    # Creo los triángulos sempiternos (laterales magenta con rojo)
    faces = np.append(faces, [[k02 + 1, k42 + 1, vu + 1]], axis=0)
    faces = np.append(faces, [[k22, k32, k42 + 2]], axis=0)

    # Creo los triangulos magenta positivos
    for i in range(1, k72 - k62 + cte):
        faces = np.append(faces, [[k42 + 1, k52 + i, k52 + i + 1]], axis=0)
        p = 0

    # creo el sello para la esquina izquierda del cuadro
    if k72 - k62 == 1:
        faces = np.append(faces, [[k52 + 1, vu + 1, k42 + 1]], axis=0)
        p = 0

    # Añade los triangulos magenta inferiores restantes
    for i in range(1, k972):
        faces = np.append(faces, [[k42 + 2, k72 + i + 1, k72 + i + 2]], axis=0)
        p = 0
    faces = np.append(faces, [[k42 + 2, k92, k72 + 1]], axis=0)

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(1, qq[0] + 1):
        faces = np.append(faces, [[k42 + 1, k92 + i + 1, k92 + i]], axis=0)
        p = 0
    for i in range(1, len(rx2) - qq[-1]):
        faces = np.append(faces, [[k42 + 2,  k102 + i + 1, k102 + i]], axis=0)
        p = 0
    # Agrega la cabeza
    faces = np.append(faces, [[k42 + 1, k92 + 1, k42 + 2]], axis=0)

#Creo el keyway en la parte de arriba
    # meto el primer punto del pe a parte porque casi siempre su #de triangulos no es igual al de los demás

    k02 = len(vertices) - 1
    vertices = np.append(vertices, [[rx22[qq[0]], ry22[qq[0]], anchoeng]], axis=0)
    k12 = len(vertices) - 1
    # Para los puntos del pe empiezo en +1 para que no se haga un desmadre ya que el primer grupo suele tener menos separaciones
    for i in range(qq[0] + 1, qq[-1] + 1):
        vertices = np.append(vertices, [[rx22[i], ry22[i], anchoeng]], axis=0)
    k22 = len(vertices) - 1

    # ppi1 es el primer punto del pi en el que se empiezan a generar las separaciones uniformes (o sea todos los triangulos verdes menos el primer y ultimo grupo)
    if len(qk) <= 1:
        ppi1 = 0
    else:
        ppi1 = len(qk) - 1
    # Notese que empieza en ppi y este puede ser 0 cuando el primer grupo es sólo una línea
    for i in range(qp[0] + ppi1, qp[-1] + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], anchoeng]], axis=0)
    k32 = len(vertices) - 1

    # Agrego el primer grupo de triángulos verdes
    for i in range(qp[0], qp[0] + ppi1 + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], anchoeng]], axis=0)
    k42 = len(vertices) - 1

    # Agrego los triángulos sempiternos (laterales magenta con rojo)
    if ppi1 == 0:
        vu = k22
    else:
        vu = k32
    vertices = np.append(vertices, [[px, py, anchoeng]], axis=0)
    vertices = np.append(vertices, [[px, -py, anchoeng]], axis=0)
    k52 = len(vertices) - 1
    # añade los triángulos magenta restantes
    vertices = np.append(vertices, [[prx, pry, anchoeng]], axis=0)
    k62 = len(vertices) - 1

    for i in range(0, qp[0] + 1):
        if ry1[i] > py:
            vertices = np.append(vertices, [[rx1[i], ry1[i], anchoeng]], axis=0)
    k72 = len(vertices) - 1
    k762 = k72 - k62
    cte = 0
    if k762 != 1:
        cte = 1

    # Añade los triangulos magenta inferiores restantes
    vertices = np.append(vertices, [[prx, -pry, anchoeng]], axis=0)
    k82 = len(vertices) - 1
    for i in range(qp[-1], len(rx1)):
        if ry1[i] < -py:
            vertices = np.append(vertices, [[rx1[i], ry1[i], anchoeng]], axis=0)
    k92 = len(vertices) - 1
    k972 = k92 - k82

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(0, qq[0] + 1):
        vertices = np.append(vertices, [[rx22[i], ry22[i], anchoeng]], axis=0)
    k102 = len(vertices) - 1

    for i in range(qq[-1], len(rx22)):
        vertices = np.append(vertices, [[rx22[i], ry22[i], anchoeng]], axis=0)
    k112 = len(vertices) - 1

    #AQUI ESTA EL ERROR
    u = 1
    #Puntos del pi sobrantes
    residuo = (k32-k22-1)%k
    iterac = k32-k22-residuo
    #Crea los tríangulos verdes
    for i in range(1,  z):
        for j in range(u + 0, u + k + 0):
            if j <(iterac):
                faces = np.append(faces, [[k12 + i, k22 + j + 1, k22 + j]], axis=0)
                p = 0
            if j<(k32-k22):
                v=u
            else:
                break
        u = u + k
    #En caso de que haya triangulos del pi que sobren ya que en el ultimo punto usado del pe sus triángulos estan incompletos, se ejecuta lo siguiente
    # para conectarlos con el punto del pe anterior
    x = qq[-1] - qq[0]
    if residuo!=0:
        for i in range(0, residuo):
            faces = np.append(faces, [[k12 + x, k22 + i + v+1, k22 + i+v]], axis=0)
            p=0

    # Creo el primer grupo de triángulos verdes
    for i in range(qq[0], qq[0] + 1):
        for j in range(1, 1 + ppi1):
            faces = np.append(faces, [[k02 + i, k32 + j + 1, k32 + j]], axis=0)
            p = 0
            if ppi1 == 0:
                print("Maldita sea Jimbo, el comunismo no funciona !")

    # Creo los complementos de los triangulos verdes para cerrar
    for i in range(1, qq[-1] - qq[0] + 1):
        faces = np.append(faces, [[k02 + i, k02 + i + 1, k22 + (i - 1) * (k) + 1]], axis=0)
        p=i

    # Creo los triángulos sempiternos (laterales magenta con rojo)
    faces = np.append(faces, [[k02 + 1, vu + 1, k42 + 1]], axis=0)
    faces = np.append(faces, [[k22, k42 + 2,k32]], axis=0)

    #Creo los triangulos magenta positivos
    for i in range(1, k72 - k62 + cte):
        faces = np.append(faces, [[k42 + 1, k52 + i + 1,k52 + i]], axis=0)
        p=0

    # creo el sello para la esquina izquierda del cuadro
    if k72 - k62 == 1:
        faces = np.append(faces, [[k52 + 1, k42 + 1, vu + 1]], axis=0)
        p=0

    # Añade los triangulos magenta inferiores restantes
    for i in range(1, k972):
        faces = np.append(faces, [[k42 + 2, k72+i+2,k72+i+1]], axis=0)
        p=0
    faces = np.append(faces, [[k42 + 2, k72 + 1, k92]], axis=0)

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(1, qq[0] + 1):
        faces = np.append(faces, [[k42 + 1, k92 + i, k92 + i + 1]], axis=0)
        p=0
    for i in range(1, len(rx2) - qq[-1]):
        faces = np.append(faces, [[k42 + 2, k102 + i, k102 + i + 1]], axis=0)
        p=0
    # Agrega la cabeza
    faces = np.append(faces, [[k42 + 1, k42 + 2, k92 + 1]], axis=0)

#Creo las caras laterales
    # triangulos magenta positivos
    for i in range(0, k7 - k6):
        faces = np.append(faces, [[k6+i, k62+i,k6+i+1]], axis=0)
        faces = np.append(faces, [[k6+i+1, k62+i, k62+i+1]], axis=0)
        p = 0
    kn = k9 - k8
    # triangulos magenta negativos
    if kn > 0:
        for i in range(1, kn):
            faces = np.append(faces, [[k8+i,k82+i,k8+i+1]], axis=0)
            faces = np.append(faces, [[k8+i+1,k82+i,k82+i+1]], axis=0)
            p = 0
    faces = np.append(faces, [[k8,k9,k82]], axis=0)
    faces = np.append(faces, [[k9,k92,k82]], axis=0)
    # casi todos los conectes del pe con casi todos los del pi
    for i in range(1, k3 - k2):
        faces = np.append(faces, [[k2+i, k22+i, k2+i+1]], axis=0)
        faces = np.append(faces, [[k2+i+1, k22+i, k22+i+1]], axis=0)
        p = 0
    # conectes del pi para el primer vertice del pe
    for i in range(1, k4 - k3):
        faces = np.append(faces, [[k3+i, k32+i, k3+i+1]], axis=0)
        faces = np.append(faces, [[k3+i+1, k32+i, k32+i+1]], axis=0)
        p = 0
    #crea las caras laterales para cerrar el keyhole
    #Izquierda
    faces = np.append(faces, [[k4+1,k42+1, k5+1]], axis=0)
    faces = np.append(faces, [[k5+1, k42+1,k52+1]], axis=0)

    #Derecha
    faces = np.append(faces, [[k52, k5,k82]], axis=0)
    faces = np.append(faces, [[ k5,k8,k82]], axis=0)

    #Frontal
    faces = np.append(faces, [[k5,k52, k4+1]], axis=0)
    faces = np.append(faces, [[k52, k42+1, k4+1]], axis=0)

    return vertices, faces


#Las funciones DO sirven tanto para doble helicoidal como para diente recto
def DOcuadrado(vertices, faces, normals, angroth, z, radiopoligono, a1, anchoeng, rcircunscrito):
    #Es la misma función que la del diente recto
    #La clave está en los angroth, búscalos y encontrarás donde necesitas las modificaciones
    if radiopoligono != -1 and rcircunscrito <= radiopoligono - .15:
        jx = []
        jy = []
        ix = []
        iy = []
        numcarasorig = len(faces)
        for i in range(0, z):
            jx.append((radiopoligono) * np.cos(i * 2 * np.pi / z -angroth))
            jy.append((radiopoligono) * np.sin(i * 2 * np.pi / z -angroth))
        for i in range(0, 4):
            ix.append(rcircunscrito * np.cos(i * 2 * np.pi / 4 -angroth))
            iy.append(rcircunscrito * np.sin(i * 2 * np.pi / 4 -angroth))
        # vertmax es el número máximo de vértices que puede conectar un punto del hexagono con los del polígono exterior
        vertmax = int((z - 1) / 4) + 2
        # cvgrandes indica qué puntos del hexágono conectan con la cantidad vertmax con vértices del polígono exterior
        cvgrandes = z - int(z / 4) * 4
        # patrones es una lista que indica qué vertices del hexágono conectan con la cantidad vertmax (si cvgrandes=0 entonces todos los vertices del hexágono tienen vertmax)
        patrones = [[0, 1, 2, 3], [0], [0, 2], [0, 1, 2], [0, 1, 2, 3]]

        # Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        # Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        # Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx[-1], jy[-1], a1]], axis=0)
        qj = len(vertices) - 1
        for i in range(0, len(jx)):
            vertices = np.append(vertices, [[jx[i], jy[i], a1]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], a1]], axis=0)
        qi = len(vertices) - 1
        for i in range(0, len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], a1]], axis=0)
        qx = len(vertices) - 1

        #region "Poligon´s inferior faces"
        a = 0
        # Ciclo for que va entre todos los puntos del hexágono
        for i in range(0, 4):
            if i in patrones[cvgrandes]:
                # Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi + 1) - 1 + i, (qi + 1) + i, (qj + 1) - 1 + a]], axis=0)
                for k in range(0, vertmax - 1):
                    # Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qj + 1) - 1 + a + k, (qi + 1) + i, (qj + 1) + a + k]], axis=0)
                # Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                # Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a = a + vertmax - 1
            else:
                # Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi + 1) - 1 + i, (qi + 1) + i, (qj + 1) - 1 + a]], axis=0)
                for m in range(0, vertmax - 2):
                    # Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qj + 1) - 1 + a + m, (qi + 1) + i, (qj + 1) + a + m]], axis=0)
                a = a + vertmax - 2
        #endregion

        # Crea el hexágono en la cara superior
        jx = []
        jy = []
        ix = []
        iy = []
        for i in range(0, z):
            jx.append((radiopoligono) * np.cos(i * 2 * np.pi / z + -angroth))
            jy.append((radiopoligono) * np.sin(i * 2 * np.pi / z + -angroth))
        for i in range(0, 4):
            ix.append(rcircunscrito * np.cos(i * 2 * np.pi / 4 + -angroth))
            iy.append(rcircunscrito * np.sin(i * 2 * np.pi / 4 + -angroth))

        # Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        # Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        # Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx[-1], jy[-1], anchoeng]], axis=0)
        qj2 = len(vertices) - 1
        for i in range(0, len(jx)):
            vertices = np.append(vertices, [[jx[i], jy[i], anchoeng]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], anchoeng]], axis=0)
        qi2 = len(vertices) - 1
        for i in range(0, len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], anchoeng]], axis=0)
        qx = len(vertices) - 1
        a = 0
        # Ciclo for que va entre todos los puntos del hexágono
        for i in range(0, 4):
            if i in patrones[cvgrandes]:
                # Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi2 + 1) + i, (qi2 + 1) - 1 + i, (qj2 + 1) - 1 + a]], axis=0)
                for k in range(0, vertmax - 1):
                    # Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qi2 + 1) + i, (qj2 + 1) - 1 + a + k, (qj2 + 1) + a + k]], axis=0)
                # Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                # Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a = a + vertmax - 1
            else:
                # Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi2 + 1) + i, (qi2 + 1) - 1 + i, (qj2 + 1) - 1 + a]], axis=0)
                for m in range(0, vertmax - 2):
                    # Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qi2 + 1) + i, (qj2 + 1) - 1 + a + m, (qj2 + 1) + a + m]], axis=0)
                a = a + vertmax - 2
        # Creo las caras laterales para el hexágono
        pisep = np.linspace(0, 2 * np.pi, 5)
        xi = rcircunscrito * np.cos(pisep-angroth)
        yi = rcircunscrito * np.sin(pisep-angroth)

        qp = len(vertices) - 1
        for i in range(0, len(pisep)):
            vertices = np.append(vertices, [[xi[i], yi[i], a1]], axis=0)
        qp1 = len(vertices) - 1
        for i in range(0, len(pisep)):
            vertices = np.append(vertices, [[xi[i], yi[i], anchoeng]], axis=0)
        for i in range(0, len(pisep)):
            faces = np.append(faces, [[qp + i + 1, qp + i, qp1 + i]], axis=0)
            faces = np.append(faces, [[qp + i + 1, qp1 + i, qp1 + i + 1]], axis=0)

        # Fixed code for normals
        for i in range(numcarasorig, len(faces)):
            we = faces[i]
            v1 = vertices[we[0]]
            v2 = vertices[we[1]]
            v3 = vertices[we[2]]
            normal = normal_to_surface(v1, v2, v3)
            normals = np.append(normals, [normal], axis=0)
    else:
        print("Radio del polígono circumscrito es muy grande, pruebe con un valor <" + str(radiopoligono))
    return vertices, faces, normals

def DOhexagono(vertices, faces, normals, angroth, z, radiopoligono, a1, anchoeng, rcircunscrito):
    if radiopoligono != -1 and rcircunscrito <= radiopoligono-.15:
        jx=[]
        jy=[]
        ix=[]
        iy=[]
        numcarasorig = len(faces)
        for i in range(0,z):
            jx.append((radiopoligono)*np.cos(i*2*np.pi/z - angroth))
            jy.append((radiopoligono)*np.sin(i*2*np.pi/z - angroth))
        for i in range(0,6):
            ix.append(rcircunscrito* np.cos(i*2*np.pi/6 - angroth))
            iy.append(rcircunscrito* np.sin(i*2*np.pi/6 - angroth))
        #vertmax es el número máximo de vértices que puede conectar un punto del hexagono con los del polígono exterior
        vertmax =int((z-1)/6)+2
        #cvgrandes indica qué puntos del hexágono conectan con la cantidad vertmax con vértices del polígono exterior
        cvgrandes = z-int(z/6)*6
        #patrones es una lista que indica qué vertices del hexágono conectan con la cantidad vertmax (si cvgrandes=0 entonces todos los vertices del hexágono tienen vertmax)
        patrones= [[0,1,2,3,4,5], [0], [0,3], [0,2,4], [0,2,3,5], [0,1,2,3,4]]

        #Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        #Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        #Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx[-1], jy[-1], a1]], axis=0)
        qj = len(vertices)-1
        for i in range(0,len(jx)):
            vertices = np.append(vertices, [[jx[i], jy[i], a1]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], a1]], axis=0)
        qi = len(vertices)-1
        for i in range(0,len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], a1]], axis=0)
        qx= len(vertices)-1

        #region "Inferior faces of the polygon"
        a=0
        #Ciclo for que va entre todos los puntos del hexágono
        for i in range(0,6):
            if i in patrones[cvgrandes]:
                #Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi+1)-1+i, (qi + 1) + i, (qj+1)-1+a]], axis=0)

                for k in range(0,vertmax-1):
                    #Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qj+1)-1+a+k, (qi + 1) + i, (qj+1)+a+k]], axis=0)
                #Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                #Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a=a+vertmax-1
            else:
                #Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi+1)-1+i, (qi + 1) + i, (qj+1)-1+a]], axis=0)

                for m in range(0,vertmax-2):
                    #Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qj+1)-1+a+m, (qi + 1) + i, (qj+1)+a+m]], axis=0)

                a=a+vertmax-2
        #endregion

        #Notese que tengo que meter el último punto de jx,jy antes de meter los demás
        #Esto porque en los for de conección para i(x) uso j(x-1) y j(x) por lo que si no lo hago j(x-1) va a fallar entregando otro vertice no deseado
        #Esto también se hizo para ix,iy ***
        vertices = np.append(vertices, [[jx[-1], jy[-1], anchoeng]], axis=0)
        qj2 = len(vertices)-1
        for i in range(0,len(jx)):
            vertices = np.append(vertices, [[jx[i], jy[i], anchoeng]], axis=0)
        # ***
        vertices = np.append(vertices, [[ix[-1], iy[-1], anchoeng]], axis=0)
        qi2 = len(vertices)-1
        for i in range(0,len(ix)):
            vertices = np.append(vertices, [[ix[i], iy[i], anchoeng]], axis=0)
        qx= len(vertices)-1
        a=0
        #Ciclo for que va entre todos los puntos del hexágono
        for i in range(0,6):
            if i in patrones[cvgrandes]:
                #Triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi2+1)+i,(qi2+1)-1+i,(qj2+1)-1+a]], axis=0)
                for k in range(0,vertmax-1):
                    #Crea los triangulos que tienen como base caras del polígono exterior
                    faces = np.append(faces, [[(qi2+1)+i,(qj2+1)-1+a+k,(qj2+1)+a+k]], axis=0)
                #Notese que 'a' me da el último vértice del polígono exterior usado, pero depende de si el vertice del héxagono usa vertmax o no
                #Por eso tiene otra parte en el 'else' para que avance de manera no lineal
                a=a+vertmax-1
            else:
                #Crea el resto de los triángulos que tienen como base al hexágono
                faces = np.append(faces, [[(qi2+1)+i,(qi2+1)-1+i,(qj2+1)-1+a]], axis=0)
                for m in range(0,vertmax-2):
                    #Crea el resto de los triángulos que tienen como base al polígono exterior
                    faces = np.append(faces, [[(qi2+1)+i,(qj2+1)-1+a+m,(qj2+1)+a+m]], axis=0)
                a=a+vertmax-2
        #Creo las caras laterales para el hexágono
        pisep=np.linspace(0,2*np.pi,7)
        xi = rcircunscrito*np.cos(pisep - angroth)
        yi = rcircunscrito*np.sin(pisep - angroth)

        qp=len(vertices)-1
        for i in range(0,len(pisep)):
            vertices = np.append(vertices, [[xi[i], yi[i], a1]], axis=0)
        qp1=len(vertices)-1
        for i in range(0,len(pisep)):
            vertices = np.append(vertices, [[xi[i], yi[i], anchoeng]], axis=0)
        for i in range(0,len(pisep)):
            faces= np.append(faces, [[qp+i+1, qp+i, qp1+i]], axis=0)
            faces = np.append(faces, [[qp+i+1, qp1+i, qp1+i+1]], axis=0)

        #Code for normals
        for i in range(numcarasorig, len(faces)):
            we = faces[i]
            v1 = vertices[we[0]]
            v2 = vertices[we[1]]
            v3 = vertices[we[2]]
            normal = normal_to_surface(v1, v2, v3)
            normals = np.append(normals, [normal], axis=0)
    else:
        print("Radio del polígono circumscrito es muy grande, pruebe con un valor <" + str(radiopoligono))
    return vertices, faces, normals

def DOkeyway(vertices, faces, normals, angroth, z, radiopoligono, anchoeng, r1, alturaC, anchoC):
    # Creates the keyhole on the bottom
    # Sets faces for normals calculations
    numcarasorig = len(faces)
    # l1 es el poligono interior y l2 el exterior
    l2 = z
    if z <= 10:
        l1 = 20 * z
    elif 10 < z <= 20:
        l1 = 10 * z
    elif 20 < z <= 50:
        l1 = 5 * z
    elif 50 < z <= 100:
        l1 = 2 * z
    else:
        l1 = z

    # caracteristicas del keyway
    altura = alturaC
    ancho = anchoC

    # radios ocasionados por el keyway
    # r1 = 10
    r2 = radiopoligono

    # puntos del keyway
    px = r1 + altura
    py = ancho / 2
    prx = np.sqrt((r1 ** 2) - ((ancho / 2) ** 2))
    pry = ancho / 2

    pisep = np.linspace(0, 2 * np.pi, l1 + 1)
    pisep2 = np.linspace(0, 2 * np.pi, l2 + 1)

    rx1 = r1 * np.cos(pisep)
    ry1 = r1 * np.sin(pisep)

    rx2 = r2 * np.cos(pisep2)
    ry2 = r2 * np.sin(pisep2)
    # Encuentro el primer vert del pe que está fuera del cuadro
    for i in range(1, len(rx2)):
        a = rx2[i]
        b = ry2[i]
        if a > 0 and b > ancho / 2:
            item = i
            break
    for i in range(1, len(rx2)):
        a = rx2[i]
        b = ry2[i]
        if a > 0 and b < 0 and b >= -ancho / 2:
            item2 = i - 1
            break
    # Encuentro los vertices del pi que están fuera del cuadro
    for i in range(1, len(rx1)):
        a = rx1[i]
        b = ry1[i]
        if a > 0 and b > ancho / 2:
            item3 = i
            break
    for i in range(1, len(rx1)):
        a = rx1[i]
        b = ry1[i]
        if a > 0 and b < -ancho / 2:
            item4 = i

    k = int(l1 / l2)
    a = int(k / 2)
    # qq es para el polígono exterior y qp para el interior
    qq = []
    qp = []
    qk = []
    qk2 = []
    tr = True
    u = a
    for i in range(item, item2 + 1):
        for j in range(u, u + k + 1):
            if j >= item3 and j <= item4:
                qq.append(i)
                qp.append(j)
                # qk es una lista para encontrar el primer vertice del polígono exterior que conecta con uno del interior ya que item2 no necesariamente sera el 1ero
                # con len(qk) puedo saber cuantas lineas tiene: si es 0 o 1 tiene una linea, para los demás numeros tiene la cantidad que indiquen
                if tr:
                    qk.append(i)
                    qk2.append(j)
        tr = False
        u = u + k

    # region "Bloque para rotar el keyway"

    """
    La lógica que usé en algunos de mis loops, si bien no está mal, no es correcta tampoco. Esto lo explico
    mejor a continuación, pero en teoría todo este bloque de validación debería estar arriba (o sea sin
    sustituir las variables originales por unas nuevas que están mejor escritas).

    Algunas características no se trasladan bien al sistema radial (ve el diagrama)
    Con esto me refiero que pej usé puntos precisos en vez de coordenadas radiales pej el npy era originalmente
    -py, eso porque el keyway como tal es simétrico respecto al eje, pero a la hora de rotar lo ni de pedo 
    puede ser -py, ahí pasa a ser npy (negative py).

    Este tipo de cosas muestran que no es que esté mal mi lógica, sólo que hay que cambiar las expresiones
    matemáticas por unas mejores que si se trasladen al rotar el agujero.
    """

    anguloi = np.arctan((ancho / 2) / (r1 + altura))
    tetai = np.arcsin((ancho / 2) / r1)
    radioAlPy = np.sqrt((ancho / 2) ** 2 + (r1 + altura) ** 2)

    """
    No todos los loops y condicionales están arriba de este bloque, tons en algunos tengo que usar el valor
    de las variables cuando el agujero es simétrico respecto al eje X (o sea el valor original)
     """
    ry10 = ry1
    py0 = py

    # angroth = 0
    px = radioAlPy * np.cos(anguloi - angroth)
    py = radioAlPy * np.sin(anguloi - angroth)
    prx = r1 * np.cos(tetai - angroth)
    pry = r1 * np.sin(tetai - angroth)

    npx = radioAlPy * np.cos(-anguloi - angroth)
    npy = radioAlPy * np.sin(-anguloi - angroth)

    nprx = r1 * np.cos(-tetai - angroth)
    npry = r1 * np.sin(-tetai - angroth)

    rx1 = r1 * np.cos(pisep - angroth)
    ry1 = r1 * np.sin(pisep - angroth)

    rx2 = r2 * np.cos(pisep2 - angroth)
    ry2 = r2 * np.sin(pisep2 - angroth)

    # endregion

    # meto el primer punto del pe a parte porque casi siempre su #de triangulos no es igual al de los demás
    k0 = len(vertices) - 1
    vertices = np.append(vertices, [[rx2[qq[0]], ry2[qq[0]], 0]], axis=0)
    k1 = len(vertices) - 1
    # Para los puntos del pe empiezo en +1 para que no se haga un desmadre ya que el primer grupo suele tener menos separaciones
    for i in range(qq[0] + 1, qq[-1] + 1):
        vertices = np.append(vertices, [[rx2[i], ry2[i], 0]], axis=0)
    k2 = len(vertices) - 1

    # ppi1 es el primer punto del pi en el que se empiezan a generar las separaciones uniformes (o sea todos los triangulos verdes menos el primer y ultimo grupo)
    if len(qk) <= 1:
        ppi1 = 0
    else:
        ppi1 = len(qk) - 1
    # Notese que empieza en ppi y este puede ser 0 cuando el primer grupo es sólo una línea
    for i in range(qp[0] + ppi1, qp[-1] + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k3 = len(vertices) - 1

    # Agrego el primer grupo de triángulos verdes
    for i in range(qp[0], qp[0] + ppi1 + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k4 = len(vertices) - 1

    # Agrego los triángulos sempiternos (laterales magenta con rojo)
    if ppi1 == 0:
        vu = k2
    else:
        vu = k3
    vertices = np.append(vertices, [[px, py, 0]], axis=0)
    vertices = np.append(vertices, [[npx, npy, 0]], axis=0)
    k5 = len(vertices) - 1
    # añade los triángulos magenta restantes
    vertices = np.append(vertices, [[prx, pry, 0]], axis=0)
    k6 = len(vertices) - 1

    for i in range(0, qp[0] + 1):
        if ry10[i] > py0:
            vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k7 = len(vertices) - 1
    k76 = k7 - k6
    cte = 0
    if k76 != 1:
        cte = 1

    # Añade los triangulos magenta inferiores restantes
    vertices = np.append(vertices, [[nprx, npry, 0]], axis=0)
    k8 = len(vertices) - 1
    print("k8", k8)
    for i in range(qp[-1], len(rx1)):
        if ry10[i] < -py0:
            vertices = np.append(vertices, [[rx1[i], ry1[i], 0]], axis=0)
    k9 = len(vertices) - 1
    k97 = k9 - k8

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(0, qq[0] + 1):
        vertices = np.append(vertices, [[rx2[i], ry2[i], 0]], axis=0)
    k10 = len(vertices) - 1

    for i in range(qq[-1], len(rx2)):
        vertices = np.append(vertices, [[rx2[i], ry2[i], 0]], axis=0)
    k11 = len(vertices) - 1

    u = 1
    # Puntos del pi sobrantes
    residuo = (k3 - k2 - 1) % k
    iterac = k3 - k2 - residuo
    # Crea los tríangulos verdes
    for i in range(1, z):
        for j in range(u + 0, u + k + 0):
            if j < (iterac):
                faces = np.append(faces, [[k2 + j + 1, k1 + i, k2 + j]], axis=0)
                p = 0
            if j < (k3 - k2):
                v = u
            else:
                break
        u = u + k
    # En caso de que haya triangulos del pi que sobren ya que en el ultimo punto usado del pe sus triángulos estan incompletos, se ejecuta lo siguiente
    # para conectarlos con el punto del pe anterior
    x = qq[-1] - qq[0]
    if residuo != 0:
        for i in range(0, residuo):
            faces = np.append(faces, [[k2 + i + v + 1, k1 + x, k2 + i + v]], axis=0)

    # Creo el primer grupo de triángulos verdes
    for i in range(qq[0], qq[0] + 1):
        for j in range(1, 1 + ppi1):
            faces = np.append(faces, [[k3 + j + 1, k0 + i, k3 + j]], axis=0)
            p = 0
            if ppi1 == 0:
                print("Maldita sea Jimbo, el comunismo no funciona !")

    # Creo los complementos de los triangulos verdes para cerrar
    for i in range(1, qq[-1] - qq[0] + 1):
        faces = np.append(faces, [[k0 + i + 1, k0 + i, k2 + (i - 1) * (k) + 1]], axis=0)

    # Creo los triángulos sempiternos (laterales magenta con rojo)
    faces = np.append(faces, [[vu + 1, k0 + 1, k4 + 1]], axis=0)
    faces = np.append(faces, [[k4 + 2, k2, k3]], axis=0)

    # Agrego los triángulos magenta positivos
    for i in range(1, k7 - k6 + cte):
        faces = np.append(faces, [[k5 + i + 1, k4 + 1, k5 + i]], axis=0)

    # creo el sello para la esquina izquierda del cuadro
    if k7 - k6 == 1:
        faces = np.append(faces, [[k4 + 1, k5 + 1, vu + 1]], axis=0)

    # Añade los triangulos magenta inferiores restantes
    for i in range(1, k97):
        faces = np.append(faces, [[k7 + i + 2, k4 + 2, k7 + i + 1]], axis=0)
        print("magenta", i)
    faces = np.append(faces, [[k7 + 1, k4 + 2, k9]], axis=0)

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(1, qq[0] + 1):
        faces = np.append(faces, [[k9 + i, k4 + 1, k9 + i + 1]], axis=0)
    for i in range(1, len(rx2) - qq[-1]):
        faces = np.append(faces, [[k10 + i, k4 + 2, k10 + i + 1]], axis=0)
    # Agrega la cabeza
    faces = np.append(faces, [[k4 + 2, k4 + 1, k9 + 1]], axis=0)

    # Creo el keyway en la parte de arriba
    # meto el primer punto del pe a parte porque casi siempre su #de triangulos no es igual al de los demás
    k02 = len(vertices) - 1
    vertices = np.append(vertices, [[rx2[qq[0]], ry2[qq[0]], anchoeng]], axis=0)
    k12 = len(vertices) - 1
    # Para los puntos del pe empiezo en +1 para que no se haga un desmadre ya que el primer grupo suele tener menos separaciones
    for i in range(qq[0] + 1, qq[-1] + 1):
        vertices = np.append(vertices, [[rx2[i], ry2[i], anchoeng]], axis=0)
    k22 = len(vertices) - 1

    # ppi1 es el primer punto del pi en el que se empiezan a generar las separaciones uniformes (o sea todos los triangulos verdes menos el primer y ultimo grupo)
    if len(qk) <= 1:
        ppi1 = 0
    else:
        ppi1 = len(qk) - 1
    # Notese que empieza en ppi y este puede ser 0 cuando el primer grupo es sólo una línea
    for i in range(qp[0] + ppi1, qp[-1] + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], anchoeng]], axis=0)
    k32 = len(vertices) - 1

    # Agrego el primer grupo de triángulos verdes
    for i in range(qp[0], qp[0] + ppi1 + 1):
        vertices = np.append(vertices, [[rx1[i], ry1[i], anchoeng]], axis=0)
    k42 = len(vertices) - 1

    # Agrego los triángulos sempiternos (laterales magenta con rojo)
    if ppi1 == 0:
        vu = k22
    else:
        vu = k32
    vertices = np.append(vertices, [[px, py, anchoeng]], axis=0)
    vertices = np.append(vertices, [[npx, npy, anchoeng]], axis=0)
    k52 = len(vertices) - 1
    # añade los triángulos magenta restantes
    vertices = np.append(vertices, [[prx, pry, anchoeng]], axis=0)
    k62 = len(vertices) - 1

    for i in range(0, qp[0] + 1):
        if ry10[i] > py0:
            vertices = np.append(vertices, [[rx1[i], ry1[i], anchoeng]], axis=0)
    k72 = len(vertices) - 1
    k762 = k72 - k62
    cte = 0
    if k762 != 1:
        cte = 1

    # Añade los triangulos magenta inferiores restantes
    vertices = np.append(vertices, [[nprx, npry, anchoeng]], axis=0)
    k82 = len(vertices) - 1
    for i in range(qp[-1], len(rx1)):
        if ry10[i] < -py0:
            vertices = np.append(vertices, [[rx1[i], ry1[i], anchoeng]], axis=0)
    k92 = len(vertices) - 1
    k972 = k92 - k82

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(0, qq[0] + 1):
        vertices = np.append(vertices, [[rx2[i], ry2[i], anchoeng]], axis=0)
    k102 = len(vertices) - 1

    for i in range(qq[-1], len(rx2)):
        vertices = np.append(vertices, [[rx2[i], ry2[i], anchoeng]], axis=0)
    k112 = len(vertices) - 1

    # AQUI ESTA EL ERROR
    u = 1
    # Puntos del pi sobrantes
    residuo = (k32 - k22 - 1) % k
    iterac = k32 - k22 - residuo
    # Crea los tríangulos verdes
    for i in range(1, z):
        for j in range(u + 0, u + k + 0):
            if j < (iterac):
                faces = np.append(faces, [[k12 + i, k22 + j + 1, k22 + j]], axis=0)
                p = 0
            if j < (k32 - k22):
                v = u
            else:
                break
        u = u + k
    # En caso de que haya triangulos del pi que sobren ya que en el ultimo punto usado del pe sus triángulos estan incompletos, se ejecuta lo siguiente
    # para conectarlos con el punto del pe anterior
    x = qq[-1] - qq[0]
    if residuo != 0:
        for i in range(0, residuo):
            faces = np.append(faces, [[k12 + x, k22 + i + v + 1, k22 + i + v]], axis=0)
            p = 0

    # Creo el primer grupo de triángulos verdes
    for i in range(qq[0], qq[0] + 1):
        for j in range(1, 1 + ppi1):
            faces = np.append(faces, [[k02 + i, k32 + j + 1, k32 + j]], axis=0)
            p = 0
            if ppi1 == 0:
                print("Maldita sea Jimbo, el comunismo no funciona !")

    # Creo los complementos de los triangulos verdes para cerrar
    for i in range(1, qq[-1] - qq[0] + 1):
        faces = np.append(faces, [[k02 + i, k02 + i + 1, k22 + (i - 1) * (k) + 1]], axis=0)
        p = i

    # Creo los triángulos sempiternos (laterales magenta con rojo)
    faces = np.append(faces, [[k02 + 1, vu + 1, k42 + 1]], axis=0)
    # PISTA 1
    faces = np.append(faces, [[k22, k42 + 2, k32]], axis=0)

    # Creo los triangulos magenta positivos
    for i in range(1, k72 - k62 + cte):
        faces = np.append(faces, [[k42 + 1, k52 + i + 1, k52 + i]], axis=0)
        p = 0

    # creo el sello para la esquina izquierda del cuadro
    if k72 - k62 == 1:
        faces = np.append(faces, [[k52 + 1, k42 + 1, vu + 1]], axis=0)
        p = 0

    # Añade los triangulos magenta inferiores restantes
    for i in range(1, k972):
        faces = np.append(faces, [[k42 + 2, k72 + i + 2, k72 + i + 1]], axis=0)
        p = 0
    faces = np.append(faces, [[k42 + 2, k72 + 1, k92]], axis=0)

    # Agrega la cabeza y el resto de triangulos del pe que no están
    for i in range(1, qq[0] + 1):
        faces = np.append(faces, [[k42 + 1, k92 + i, k92 + i + 1]], axis=0)
        p = 0
    for i in range(1, len(rx2) - qq[-1]):
        faces = np.append(faces, [[k42 + 2, k102 + i, k102 + i + 1]], axis=0)
        p = 0
    # Agrega la cabeza
    faces = np.append(faces, [[k42 + 1, k42 + 2, k92 + 1]], axis=0)

    # Creo las caras laterales
    # triangulos magenta positivos
    for i in range(0, k7 - k6):
        faces = np.append(faces, [[k6 + i, k62 + i, k6 + i + 1]], axis=0)
        faces = np.append(faces, [[k6 + i + 1, k62 + i, k62 + i + 1]], axis=0)
        p = 0
    kn = k9 - k8
    # triangulos magenta negativos
    if kn > 0:
        for i in range(1, kn):
            faces = np.append(faces, [[k8 + i, k82 + i, k8 + i + 1]], axis=0)
            faces = np.append(faces, [[k8 + i + 1, k82 + i, k82 + i + 1]], axis=0)
            p = 0
    faces = np.append(faces, [[k8, k9, k82]], axis=0)
    faces = np.append(faces, [[k9, k92, k82]], axis=0)
    # casi todos los conectes del pe con casi todos los del pi
    for i in range(1, k3 - k2):
        faces = np.append(faces, [[k2 + i, k22 + i, k2 + i + 1]], axis=0)
        faces = np.append(faces, [[k2 + i + 1, k22 + i, k22 + i + 1]], axis=0)
        p = 0
    # conectes del pi para el primer vertice del pe
    for i in range(1, k4 - k3):
        faces = np.append(faces, [[k3 + i, k32 + i, k3 + i + 1]], axis=0)
        faces = np.append(faces, [[k3 + i + 1, k32 + i, k32 + i + 1]], axis=0)
        p = 0
    # crea las caras laterales para cerrar el keyhole
    # Izquierda
    faces = np.append(faces, [[k4 + 1, k42 + 1, k5 + 1]], axis=0)
    faces = np.append(faces, [[k5 + 1, k42 + 1, k52 + 1]], axis=0)

    # Derecha
    faces = np.append(faces, [[k52, k5, k82]], axis=0)
    faces = np.append(faces, [[k5, k8, k82]], axis=0)

    # Frontal
    faces = np.append(faces, [[k5, k52, k4 + 1]], axis=0)
    faces = np.append(faces, [[k52, k42 + 1, k4 + 1]], axis=0)

    # Code for normals
    for i in range(numcarasorig, len(faces)):
        we = faces[i]
        v1 = vertices[we[0]]
        v2 = vertices[we[1]]
        v3 = vertices[we[2]]
        normal = normal_to_surface(v1, v2, v3)
        normals = np.append(normals, [normal], axis=0)

    return vertices, faces, normals

def DOcirculo(vertices,faces,normals, angroth, z,radiopoligono,a1,anchoeng,rcircunscrito):

    #el circle_points es porque z es 10, pero esto debe quedar de forma que x/z siempre de un entero
    if radiopoligono != -1 and rcircunscrito <= radiopoligono-.15:
        listota = []
        numcarasorig = len(faces)


        #NOTA: SI LOS PUNTOS DE LA CARA SUPERIOR NO CONECTAN BIEN, PUEDES HACER NEGATIVO EL arot Y ESO DEBERIA ARREGLARLO


        #ATENCION !!!! TIENES QUE
        # MODFICAR CIRCLE POINTS PARA QUE SEA UNA VARABLE MAS EFCIENTE A LA HORA DE SELECCIONAR LA CANTIDAD DE
        # CARAS PARA EL CIRCULO
        circle_points = int(300/z)*z
        circle_points_per_vertex = circle_points / z
        half_angle = (circle_points_per_vertex * 1 * np.pi / circle_points)


        # Para esta función estoy usando "mal" a1 porque esa se planeaba para que la cara del agujero no estuviera al mismo nivel que las caras inefirores (como otro agujero o que estuviera sumida)
        # a1 la voy a usar para decirle si son las caras de abajo o la de arriba
        def base_faces(vertices,faces,z,radiopoligono,rcircunscrito,a1,normals_downwards,arot):
            #x=normals_downwards
            jx = []
            jy = []
            ix = []
            iy = []
            for i in range(0, z):
                jx.append((radiopoligono) * np.cos(i * 2 * np.pi / z + angroth))
                jy.append((radiopoligono) * np.sin(i * 2 * np.pi / z + angroth))
            # agrego otra vez cuando el angulo=0 para que en el loop al dar la vuelta llegue al ultimo vertice del circulo (que es el inicial porque 2pi=0)
            jx.append((radiopoligono) * np.cos(0 * 2 * np.pi / z + angroth))
            jy.append((radiopoligono) * np.sin(0 * 2 * np.pi / z + angroth))
            for i in range(0, circle_points):
                ix.append(rcircunscrito * np.cos((i * 2 * np.pi / circle_points) - half_angle + angroth))
                iy.append(rcircunscrito * np.sin((i * 2 * np.pi / circle_points) - half_angle + angroth))
            # agrego otra vez cuando el angulo=0 para que en el loop al dar la vuelta llegue al ultimo vertice del circulo (que es el inicial porque 2pi=0)
            ix.append(rcircunscrito * np.cos((0 * 2 * np.pi / circle_points) - half_angle + angroth))
            iy.append(rcircunscrito * np.sin((0 * 2 * np.pi / circle_points) - half_angle + angroth))

            carasporvt = int(circle_points / z)

            qj = len(vertices) - 0
            for i in range(0, len(jx)):
                vertices = np.append(vertices, [[jx[i], jy[i], a1]], axis=0)
            # ***
            vertices = np.append(vertices, [[ix[-1], iy[-1], a1]], axis=0)
            qi = len(vertices) - 0
            for i in range(0, len(ix)):
                vertices = np.append(vertices, [[ix[i], iy[i], a1]], axis=0)
            qx = len(vertices) - 1

            for i in range(0, z):
                for j in range(0, carasporvt +0):
                    x = i * z + j +2*i
                    x = int(i*circle_points_per_vertex + j)
                    if normals_downwards:
                        faces = np.append(faces, [[qj + i, qi + x, qi + 1 + x]], axis=0)
                    else:
                        faces = np.append(faces, [[qi + x, qj + i, qi + 1 + x]], axis=0)
                if normals_downwards:
                    faces = np.append(faces, [[qj + i + 1, qj + i, qi + 1 + x]], axis=0)
                else:
                    faces = np.append(faces, [[qj + i, qj + 1 + i, qi + 1 + x]], axis=0)
            return vertices, faces, qi

        def lateral_faces(faces, qi1, qi2, z):
            for i in range(0, circle_points):
                faces = np.append(faces, [[qi1 + i, qi2 + i, qi1 + i + 1]], axis=0)
                faces = np.append(faces, [[qi1 + i + 1, qi2 + i, qi2 + i + 1]], axis=0)
            return faces

        vf = base_faces(vertices, faces, z, radiopoligono, rcircunscrito, a1, True, angroth)
        vertices = vf[0]
        faces = vf[1]
        qi1 = vf[2]
        vf = base_faces(vertices, faces, z, radiopoligono, rcircunscrito, anchoeng, False, angroth)
        vertices = vf[0]
        faces = vf[1]
        qi2 = vf[2]


        faces = lateral_faces(faces, qi1, qi2, z)
        # Code for normals
        for i in range(numcarasorig, len(faces)):
            we = faces[i]
            v1 = vertices[we[0]]
            v2 = vertices[we[1]]
            v3 = vertices[we[2]]
            normal = normal_to_surface(v1, v2, v3)
            normals = np.append(normals, [normal], axis=0)
    else:
        print("Radio del polígono circumscrito es muy grande, pruebe con un valor <" + str(radiopoligono))
    return vertices, faces, normals



#endregion

#Master class
class HoleData:

    def __init__(self, isHoleNeeded):
        self.Exists = isHoleNeeded

#region Hole sub classes
class PolygonalHole(HoleData):
    def __init__(self, isHoleNeeded, Type, Circumradius):
        super().__init__(isHoleNeeded)
        self.type = Type
        self.circumradius = Circumradius


class KeywayHole(HoleData):
    def __init__(self, isHoleNeeded, Type, BoreDiameter, KeywayWidth, BorePlusKwDepth):
        super().__init__(isHoleNeeded)
        self.type = Type
        self.boreDiameter = BoreDiameter
        self.keywayWidth = KeywayWidth
        self.BPKwD = BorePlusKwDepth
        # AlturaCuñero ver la imagen KeyWayParameters
        self._alturaC = BorePlusKwDepth - BoreDiameter
        # Desde el centro del bore hasta una de las esquinas del cuñero
        self._circumradius = np.sqrt((BoreDiameter + self._alturaC)**2 + (KeywayWidth/2)**2)




#endregion

#Master class
class Gear:

    def __init__(self, modulo, apdeg, z, anchoeng, HoleData):
        self.bytesIOstl = None
        self.modulo = modulo
        self.apdeg = apdeg
        self.z = z
        self.anchoeng = anchoeng
        self.__HoleData = HoleData

        #region Variables que puedo cambiar
            #a1 se tiene que eliminar del código
        self.__a1 = 0

            # separaciones define la calidad de la involuta
        self.__separaciones = 10

        self.__sep2 = 4
        self.__sep3 = 10
        X = 0

        #endregion

        #region Variables inmodificables
        self.__ap = degToRad(apdeg)
        self.__T = mt.pi * modulo / 2

        #endregion

        # region Diámetros y radios
        self.__dp = modulo * z
        self.__db = self.__dp * mt.cos(self.__ap)
        self.__da = self.__dp + 2 * modulo
        self.__df = self.__dp - 2.5 * modulo

        self.__dv = self.__dp + 2 * X * modulo
        self.__dva = self.__dv + 2 * modulo
        self.__dvf = self.__dv - 2.5 * modulo

            # Radios
        self.__rb = self.__db / 2
        self.__rf = self.__df / 2
        self.__rva = self.__dva / 2
        self.__rvf = self.__dvf / 2

            # Propiedades
        self.dp = self.__dv
        self.da = self.__dva
        self.db = self.__db
        self.df = self.__dvf

        # endregion

        #region Ángulos
        self.__alpha = (mt.sqrt((self.__dp * self.__dp) - (self.__db * self.__db))) / (self.__db) - self.__ap
        self.__beta = (mt.pi / (2 * self.z))
        self.__tetarva = rollAngle(self.__rva, self.__rb)
        self.__angrotp = 2 * mt.pi - self.__alpha - self.__beta

            # region Ángulos para la separación equidistante de la involuta
            # self.__L es la longitud de arco de la involuta.
            # Teniendo la longitud total, es sencillo encontrar self.__separaciones equidistantes y por lo tanto, los angulos equidistantes (lo que vienen siendo self.__tetarvasep).
        self.__L = (self.__rb / 2) * self.__tetarva ** 2
        self.__Lsep = np.linspace(0, self.__L, self.__separaciones)
        listan = []
        for i in range(0, self.__separaciones):
            angulo = mt.sqrt(self.__Lsep[i] * 2 / self.__rb)
            listan.append(angulo)
        self.__tetarvasep = np.array(listan)
            # endregion

        self.__az = 2 * np.pi / self.z
        self.__interdentsep = np.linspace((-self.__alpha - self.__beta), (self.__alpha + self.__beta), self.__sep2)
        self.__sigm = sigma(self.__rva, self.__rb, self.__dva, self.__T, self.__dv, self.__ap)
        self.__sigsep = np.linspace(-self.__sigm / 2, self.__sigm / 2, self.__sep3)
        if self.__rvf < self.__rb:
            self.__tsep = np.linspace(self.__alpha + self.__beta, self.__az - (self.__alpha + self.__beta), self.__sep3)
        else:
            self.__ttt = sigma(self.__rvf, self.__rb, self.__dvf, self.__T, self.__dv, self.__ap)
            self.__tsep = np.linspace(self.__ttt / 2, self.__az - self.__ttt / 2, self.__sep3)
            self.__tetarvf = rollAngle(self.__rvf, self.__rb)

            # region Ángulos para la separación equidistante de la involuta
            # self.__L es la longitud de arco de la involuta.
            # Teniendo la longitud total, es sencillo encontrar self.__separaciones equidistantes y por lo tanto, los angulos equidistantes (lo que vienen siendo self.__tetarvasep).
            self.__L = (self.__rb / 2) * self.__tetarva ** 2
            # L2 es la longitud del arco al rf (que en este caso está por encima del de base, entonces causaría interferencias) y como la involuta empieza ahí y no en el self.__rb se lo resto.
            L2 = (self.__rb / 2) * self.__tetarvf ** 2
            self.__Lsep = np.linspace(L2, self.__L, self.__separaciones)
            listan = []
            for i in range(0, self.__separaciones):
                angulo = mt.sqrt(self.__Lsep[i] * 2 / self.__rb)
                listan.append(angulo)
            self.__tetarvasep = np.array(listan)
            # endregion

        self.__puntox = self.__rb * (np.cos(self.__tetarvasep + self.__angrotp) + self.__tetarvasep * np.sin(
            self.__tetarvasep + self.__angrotp))
        self.__puntor = self.__puntox[-1]
        # endregion

    def GetFacets(self):
        pass

    def SetHole(self):
        pass


    def savetoBinary(self, route, name):
        self.GetFacets()
        a = "SpurGear made by GFSI [JMGO] all rights reserved."
        header = "{:<80}".format(a)
        lenfaces = len(self._Gear__faces)
        print("caras",lenfaces)
        print("normales", len(self._Gear__normals))
        bytelf = lenfaces.to_bytes(4, "little")
        ruta = route +"\\"+name + ".stl"
        #ruta = r"C:\Users\juang\Documents\Programación\Web Projects\mythreejsproject\TEXTFILE2.stl"

        # 'w' es por write only
        with open(ruta, "w") as writer:
            writer.write(header)

        # 'ab' es por append binary https://stackoverflow.com/questions/16208206/confused-by-python-file-mode-w
        with open(ruta, "ab") as writer:
            writer.write(bytelf)
            for i in range(0, len(self._Gear__faces)):
                w2 = self._Gear__normals[i]
                n1 = w2[0]
                n2 = w2[1]
                n3 = w2[2]
                # writer.write(struct.pack('%sf' % len(w2), *w2))
                writer.write(struct.pack('<f', n1))
                writer.write(struct.pack('<f', n2))
                writer.write(struct.pack('<f', n3))
                v1 = self._Gear__vertices[self._Gear__faces[i][0]]
                v1x = v1[0]
                v1y = v1[1]
                v1z = v1[2]
                writer.write(struct.pack('<f', v1x))
                writer.write(struct.pack('<f', v1y))
                writer.write(struct.pack('<f', v1z))

                v2 = self._Gear__vertices[self._Gear__faces[i][1]]
                v2x = v2[0]
                v2y = v2[1]
                v2z = v2[2]
                writer.write(struct.pack('<f', v2x))
                writer.write(struct.pack('<f', v2y))
                writer.write(struct.pack('<f', v2z))

                v3 = self._Gear__vertices[self._Gear__faces[i][2]]
                v3x = v3[0]
                v3y = v3[1]
                v3z = v3[2]
                writer.write(struct.pack('<f', v3x))
                writer.write(struct.pack('<f', v3y))
                writer.write(struct.pack('<f', v3z))
                writer.write(b'  ')

    def savetoBinaryIO(self, route, name):
        writer = BytesIO()
        self.GetFacets()
        a = "SpurGear made by GFSI [JMGO] all rights reserved."
        # El header es un requisito de 80 caracteres del stl, es preferible que sea legible (aunque puedes cifrarlo).
        header = "{:<80}".format(a)
        lenfaces = len(self._Gear__faces)
        faces_amount = lenfaces.to_bytes(4, "little")
        ruta = route + "\\" + name + ".stl"

        # nota como "conviertes" el header a bytes aunque realmente sigue siendo el texto como lo conoces
        writer.write(bytes(header, 'utf-8'))

        writer.write(faces_amount)
        for i in range(0, len(self._Gear__faces)):
            w2 = self._Gear__normals[i]
            n1 = w2[0]
            n2 = w2[1]
            n3 = w2[2]
            # writer.write(struct.pack('%sf' % len(w2), *w2))
            writer.write(struct.pack('<f', n1))
            writer.write(struct.pack('<f', n2))
            writer.write(struct.pack('<f', n3))
            v1 = self._Gear__vertices[self._Gear__faces[i][0]]
            v1x = v1[0]
            v1y = v1[1]
            v1z = v1[2]
            writer.write(struct.pack('<f', v1x))
            writer.write(struct.pack('<f', v1y))
            writer.write(struct.pack('<f', v1z))

            v2 = self._Gear__vertices[self._Gear__faces[i][1]]
            v2x = v2[0]
            v2y = v2[1]
            v2z = v2[2]
            writer.write(struct.pack('<f', v2x))
            writer.write(struct.pack('<f', v2y))
            writer.write(struct.pack('<f', v2z))

            v3 = self._Gear__vertices[self._Gear__faces[i][2]]
            v3x = v3[0]
            v3y = v3[1]
            v3z = v3[2]
            writer.write(struct.pack('<f', v3x))
            writer.write(struct.pack('<f', v3y))
            writer.write(struct.pack('<f', v3z))
            writer.write(b'  ')

        self.bytesIOstl = writer

        # Para guardar directo en unarchivo puedes usar el writer directamente o el objeto
        with open(ruta, "wb") as f:
            #f.write(writer.getbuffer())
            f.write(self.bytesIOstl.getbuffer())
        # NOTA IMPORTANTE:
        #   Para el Gear Maker lo que vas a regresar el self.bytesIOstl y luego le vas a dar self.bytesIOstl.getvalue()

    def savetoAscii(self, route, name):
        self.GetFacets()
        a = "SpurGear made by GFSI [JMGO] all rights reserved."
        header = "{:<80}".format(a)
        lenfaces = len(self._Gear__faces)
        print("caras",lenfaces)
        print("normales", len(self._Gear__normals))
        bytelf = lenfaces.to_bytes(4, "little")
        ruta = route +"\\"+name + ".stl"
        #ruta = r"C:\Users\juang\Documents\Programación\Web Projects\mythreejsproject\TEXTFILE2.stl"

        # 'w' es por write only
        with open(ruta, "w") as writer:
            writer.write("solid Juan \n")

        # 'ab' es por append binary https://stackoverflow.com/questions/16208206/confused-by-python-file-mode-w
        with open(ruta, "a") as writer:
            #writer.write(bytelf)
            for i in range(0, len(self._Gear__faces)):
                w2 = self._Gear__normals[i]

                n1 = round(w2[0],6)
                n2 = round(w2[1],6)
                n3 = round(w2[2],6)
                n1 = str("{:.6f}".format(n1))
                n2 = str("{:.6f}".format(n2))
                n3 = str("{:.6f}".format(n3))
                # writer.write(struct.pack('%sf' % len(w2), *w2))
                writer.write("facet normal " + n1 + " " + n2 + " " + n3)
                writer.write("\n")
                writer.write("outer loop\n")
                #writer.write(struct.pack('<f', n2))
                #writer.write(struct.pack('<f', n3))
                v1 = self._Gear__vertices[self._Gear__faces[i][0]]
                v1x = round(v1[0],6)
                v1y = round(v1[1],6)
                v1z = round(v1[2],6)
                v1x = str("{:.6f}".format(v1x))
                v1y = str("{:.6f}".format(v1y))
                v1z = str("{:.6f}".format(v1z))
                writer.write("vertex " + v1x + " " + v1y + " " + v1z)
                writer.write("\n")
                #writer.write(struct.pack('<f', v1y))
                #writer.write(struct.pack('<f', v1z))

                v2 = self._Gear__vertices[self._Gear__faces[i][1]]
                v2x = round(v2[0],6)
                v2y = round(v2[1],6)
                v2z = round(v2[2],6)
                v2x = str("{:.6f}".format(v2x))
                v2y = str("{:.6f}".format(v2y))
                v2z = str("{:.6f}".format(v2z))
                writer.write("vertex " + v2x + " " + v2y + " " + v2z)
                writer.write("\n")
                #writer.write(struct.pack('<f', v2x))
                #writer.write(struct.pack('<f', v2y))
                #writer.write(struct.pack('<f', v2z))

                v3 = self._Gear__vertices[self._Gear__faces[i][2]]
                v3x = round(v3[0],6)
                v3y = round(v3[1],6)
                v3z = round(v3[2],6)
                v3x = str("{:.6f}".format(v3x))
                v3y = str("{:.6f}".format(v3y))
                v3z = str("{:.6f}".format(v3z))
                #writer.write(struct.pack('<f', v3x))
                #writer.write(struct.pack('<f', v3y))
                #writer.write(struct.pack('<f', v3z))
                #writer.write(b'  ')
                writer.write("vertex " + v3x + " " + v3y + " " + v3z)
                writer.write("\n")
                writer.write("endloop")
                writer.write("\n")
                writer.write("endfacet\n")

            writer.write("endsolid Juan")

#region Gear Sub Classes
class SpurGear(Gear):

    def GetFacets(self):
        start_time2 = time.time()
        # region ESTO AQUI VA
        # debido a la inserción de las cuatro involutas en un mismo for, pasado de cierto punto (que lo marca quality) se requiere insertarlas en un orden de magnitud menor (ver el cuaderno)
        # o sea que cuando i llegaba la mitad de las separaciones empezaba a comerse caras que no eran de las cuatro involutas
        quality = (np.ceil(self._Gear__separaciones / 2) - 1)
        a1 = 0
        if self._Gear__HoleData.Exists:
            self.radiopoligono = self._Gear__rvf - .5
        else:
            self.radiopoligono = 0
        # endregion

        # region Código enorme
        # region Declaración de vertices y caras
        # Comienzo a crear el engrane
        self._Gear__vertices = np.array([[0, 0, a1]])
        self._Gear__normals = np.array([[0, 0, 0]])
        self._Gear__faces = np.array([[0, 0, 0]])

        az = 2 * np.pi / self.z
        k = 1

        # Involuta parte hacia arriba (que en realidad va abajo)
        xinvar = self._Gear__rb * (np.cos(self._Gear__tetarvasep + self._Gear__angrotp + k * az) + self._Gear__tetarvasep * np.sin(
            self._Gear__tetarvasep + self._Gear__angrotp + k * az))
        yinvar = self._Gear__rb * (np.sin(self._Gear__tetarvasep + self._Gear__angrotp + k * az) - self._Gear__tetarvasep * np.cos(
            self._Gear__tetarvasep + self._Gear__angrotp + k * az))
        # Involuta parte hacia abajo (que en realidad va arriba)
        xinvab = self._Gear__rb * (np.cos(self._Gear__tetarvasep + self._Gear__angrotp - k * az) + self._Gear__tetarvasep * np.sin(
            self._Gear__tetarvasep + self._Gear__angrotp - k * az))
        yinvab = -self._Gear__rb * (np.sin(self._Gear__tetarvasep + self._Gear__angrotp - k * az) - self._Gear__tetarvasep * np.cos(
            self._Gear__tetarvasep + self._Gear__angrotp - k * az))
        # Arco de adendo para el diente
        xarc = self._Gear__rva * np.cos(self._Gear__sigsep + k * az)
        yarc = self._Gear__rva * np.sin(self._Gear__sigsep + k * az)
        # Arco de dedendo interdental
        xarcf = self._Gear__rvf * np.cos(self._Gear__interdentsep + k * az)
        yarcf = self._Gear__rvf * np.sin(self._Gear__interdentsep + k * az)
        # Creo los vertices para el arco contiguo
        xarco = self._Gear__rvf * np.cos(self._Gear__tsep + k * az)
        yarco = self._Gear__rvf * np.sin(self._Gear__tsep + k * az)

        # SECCION self._Gear__vertices
        lastijuan = len((self._Gear__vertices))
        for i in range(0, len(xinvar)):
            self._Gear__vertices = np.append(self._Gear__vertices, [[xinvar[i], yinvar[i], a1]], axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices, [[xinvab[i], yinvab[i], a1]], axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices, [[xinvar[i], yinvar[i], self.anchoeng]], axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices, [[xinvab[i], yinvab[i], self.anchoeng]], axis=0)
            # j=i*4+1
            beta = 8 * i + lastijuan - 1
            if i < quality:
                # Caras de abajo
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 2, beta + 5]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 2, beta + 6, beta + 5]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 5, beta + 6, beta + 9]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 6, beta + 10, beta + 9]], axis=0)
                # Caras de arriba
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 3, beta + 7, beta + 4]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 4, beta + 7, beta + 8]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 7, beta + 11, beta + 8]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 8, beta + 11, beta + 12]], axis=0)
                # Laterales de izquierda (viendo de frente)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 5, beta + 3]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 3, beta + 5, beta + 7]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 5, beta + 9, beta + 7]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 7, beta + 9, beta + 11]], axis=0)
                # Laterales de derecha (viendo de frente)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 2, beta + 4, beta + 6]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 4, beta + 8, beta + 6]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 6, beta + 8, beta + 10]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 8, beta + 12, beta + 10]], axis=0)
            elif i == quality:
                # Caras de abajo
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 2, beta + 5]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 2, beta + 6, beta + 5]], axis=0)
                # Caras de arriba
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 3, beta + 7, beta + 4]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 4, beta + 7, beta + 8]], axis=0)
                # Laterales de izquierda (viendo de frente)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 5, beta + 3]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 3, beta + 5, beta + 7]], axis=0)
                # Laterales de derecha (viendo de frente)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 2, beta + 4, beta + 6]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 4, beta + 8, beta + 6]], axis=0)
                pass

        # Arco de adendo
        self._Gear__vertices = np.append(self._Gear__vertices,
                                    [[self._Gear__puntor * np.cos(k * az), self._Gear__puntor * np.sin(k * az), a1]], axis=0)
        self._Gear__vertices = np.append(self._Gear__vertices,
                                    [[self._Gear__puntor * np.cos(k * az), self._Gear__puntor * np.sin(k * az), self.anchoeng]],
                                    axis=0)
        lastijuan = len(self._Gear__vertices)
        for i in range(0, len(xarc)):
            self._Gear__vertices = np.append(self._Gear__vertices, [[xarc[i], yarc[i], a1]], axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices, [[xarc[i], yarc[i], self.anchoeng]], axis=0)
            if i < len(xarc) - 1:
                beta = i * 2 + lastijuan
                # Caras de abajo
                self._Gear__faces = np.append(self._Gear__faces, [[lastijuan - 2, beta + 2, beta]], axis=0)
                # Caras de arriba
                self._Gear__faces = np.append(self._Gear__faces, [[lastijuan - 1, beta + 1, beta + 3]], axis=0)
                # Caras frontales con la base arriba
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 2, beta + 3]], axis=0)
                # Caras frontales con base abajo
                self._Gear__faces = np.append(self._Gear__faces, [[beta, beta + 2, beta + 1]], axis=0)

        lastijuan = len(self._Gear__vertices)
        for i in range(0, len(xarco)):
            self._Gear__vertices = np.append(self._Gear__vertices, [[xarco[i], yarco[i], a1]], axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices, [[xarco[i], yarco[i], self.anchoeng]], axis=0)
            if i < len(xarc) - 1:
                beta = i * 2 + lastijuan
                # Caras de abajo
                self._Gear__faces = np.append(self._Gear__faces, [[lastijuan + 2 * len(xarco), beta + 2, beta]], axis=0)
                # Caras de arriba
                self._Gear__faces = np.append(self._Gear__faces, [[lastijuan + 2 * len(xarco) + 1, beta + 1, beta + 3]],
                                         axis=0)
                # Caras frontales con la base arriba
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 2, beta + 3]], axis=0)
                # Caras frontales con base abajo
                self._Gear__faces = np.append(self._Gear__faces, [[beta, beta + 2, beta + 1]], axis=0)

        # Meto los centros para crear el polígono exterior (debe haber un if por si el radio es 0)
        self._Gear__vertices = np.append(self._Gear__vertices,
                                    [[(self.radiopoligono) * np.cos(k * az), (self.radiopoligono) * np.sin(k * az), a1]], axis=0)
        self._Gear__vertices = np.append(self._Gear__vertices, [
            [(self.radiopoligono) * np.cos(k * az), (self.radiopoligono) * np.sin(k * az), self.anchoeng]], axis=0)

        # Creo el triangulo del origen rpol a las involutas
        self._Gear__faces = np.append(self._Gear__faces, [
            [lastijuan + 2 * len(xarco), lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 1,
             lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 2]], axis=0)
        self._Gear__faces = np.append(self._Gear__faces, [
            [lastijuan + 2 * len(xarco) + 1, lastijuan - 2 * len(xarc) - 4 * len(xinvar),
             lastijuan - 2 * len(xarc) - 4 * len(xinvar) + 1]], axis=0)

        lastijuan2 = len(self._Gear__vertices)
        # Selladores de los arcos contiguos con los origenes del self.radiopoligono
        # Puntos del  siguiente arco contiguo
        self._Gear__vertices = np.append(self._Gear__vertices, [[xarcf[0], yarcf[0], a1]], axis=0)
        self._Gear__vertices = np.append(self._Gear__vertices, [[xarcf[0], yarcf[0], self.anchoeng]], axis=0)

        if self._Gear__rvf < self._Gear__rb:
            # Caras del radio de raiz a las involutas

            # Caras laterales
            # Derecha base abajo
            self._Gear__faces = np.append(self._Gear__faces, [[lastijuan, lastijuan - 2 * len(xarc) - 4 * len(xinvar) + 1,
                                                     lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 1]], axis=0)
            # Derecha base arriba
            self._Gear__faces = np.append(self._Gear__faces,
                                     [[lastijuan + 1, lastijuan - 2 * len(xarc) - 4 * len(xinvar) + 1, lastijuan]],
                                     axis=0)
            # Izquierda
            # Base abajo
            self._Gear__faces = np.append(self._Gear__faces, [[lastijuan2, lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 2,
                                                     lastijuan - 2 * len(xarc) - 4 * len(xinvar)]], axis=0)
            # Base arriba
            self._Gear__faces = np.append(self._Gear__faces,
                                     [[lastijuan2, lastijuan - 2 * len(xarc) - 4 * len(xinvar), lastijuan2 + 1]],
                                     axis=0)

        if self.radiopoligono != 0 and self.radiopoligono != -1:
            # Selladores del radio de raiz a la primera involuta
            # Caras de abajo
            self._Gear__faces = np.append(self._Gear__faces, [
                [lastijuan + 2 * len(xarco), lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 2, lastijuan2]], axis=0)
            # Derecha
            self._Gear__faces = np.append(self._Gear__faces, [
                [lastijuan + 2 * len(xarco), lastijuan, lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 1]], axis=0)

            # Caras de arriba
            # Izquierda
            self._Gear__faces = np.append(self._Gear__faces, [
                [lastijuan + 2 * len(xarco) + 1, lastijuan2 + 1, lastijuan - 2 * len(xarc) - 4 * len(xinvar)]], axis=0)
            # Derecha
            self._Gear__faces = np.append(self._Gear__faces, [
                [lastijuan + 2 * len(xarco) + 1, lastijuan - 2 * len(xarc) - 4 * len(xinvar) + 1, lastijuan + 1]],
                                     axis=0)

            # Creo los complementos triangulares para los arcos contiguos que siguen
            self._Gear__vertices = np.append(self._Gear__vertices, [
                [(self.radiopoligono) * np.cos(k * az + az), (self.radiopoligono) * np.sin(k * az + az), a1]], axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices, [
                [(self.radiopoligono) * np.cos(k * az + az), (self.radiopoligono) * np.sin(k * az + az), self.anchoeng]], axis=0)

            # Triangulotes entre centros rpol
            # El de abajo
            self._Gear__faces = np.append(self._Gear__faces, [[lastijuan + 2 * len(xarco), lastijuan2 + 2, lastijuan2 - 4]],
                                     axis=0)
            # El de arriba
            self._Gear__faces = np.append(self._Gear__faces,
                                     [[lastijuan + 2 * len(xarco) + 1, lastijuan2 - 3, lastijuan2 + 3]], axis=0)

        # Elimino la primera cara porque es inutil (unía el origen con el origen con el origen)
        self._Gear__faces = np.delete(self._Gear__faces, 0, 0)

        # endregion

        for i in range(0, len(self._Gear__faces)):
            we = self._Gear__faces[i]
            v1 = self._Gear__vertices[we[0]]
            v2 = self._Gear__vertices[we[1]]
            v3 = self._Gear__vertices[we[2]]
            normal = normal_to_surface(v1, v2, v3)
            self._Gear__normals = np.append(self._Gear__normals, [normal], axis=0)

        self._Gear__normals = np.delete(self._Gear__normals, 0, 0)
        norfinito = self._Gear__normals
        # Conversor a binerio es extremadamente importante. Crea una lista con los exponentes para 2**x de manera que controla la cantidad de dientes
        lista = conversoraBinario(self.z)

        self._Gear__normals = normales1(lista, az, self._Gear__normals)

        # region Inserción de dientes
        # Esta parte es crucial para los insertos de dientes para que se introduzcan en una matriz temporal y luego en la de verdad
        facejuan3 = self._Gear__faces
        vertijuan3 = self._Gear__vertices
        vertieterno = self._Gear__vertices
        faceterno = self._Gear__faces

        # Creo el inserto del primer grupo de dientes (recuerda que insertodientes1 afectamente la matriz final mientras que insertodientes2 no lo hace hasta el final
        dientes = insertodientes1(lista, az, self._Gear__vertices, self._Gear__faces, vertijuan3, facejuan3, vertieterno,
                                  faceterno)
        self._Gear__vertices = dientes[0]
        self._Gear__faces = dientes[1]
        # print("sas", len(self._Gear__normals), len(self._Gear__faces))
        self._Gear__normals = normales2(lista, az, self._Gear__normals, norfinito, norfinito)

        # Creo los insertos para el resto de los dientes
        dientes = insertodientes2(lista, az, self._Gear__vertices, self._Gear__faces, vertijuan3, facejuan3, vertieterno,faceterno)
        self._Gear__vertices = dientes[0]
        self._Gear__faces = dientes[1]
        # endregion

        # Imprime en la pantalla la cantidad de dientes que tiene el engrane realmente
        # print("zzz",len(self._Gear__vertices)/len(vertieterno),len(self._Gear__faces))

        # print("Optimizado2 en %s seconds" % (time.time() - start_time2))
        # Bloque para hacer la malla
        # endregion

        #AGUJERO
        #listaaguj = cuadrado(self._Gear__vertices, self._Gear__faces, self._Gear__normals,self.z,self.radiopoligono,0,self.anchoeng,5)
        #listaaguj = DOcirculo(self._Gear__vertices, self._Gear__faces, self._Gear__normals, 0, self.z,self.radiopoligono, 0, self.anchoeng, 5)
        #r1, alturaC, anchoC

        listaaguj = [self._Gear__vertices, self._Gear__faces, self._Gear__normals]
        # IMPORTANTE !! ESTO JALARIA COMO FUNCION O METODO porque DH Y DR comparten algoritmos DO solo HL es diferente
        # NOTA: HAZ QUE MATCHEEN CON LO QUE ARROJA TU UI -> HEX, SQR, CIR, KW, None
        if self._Gear__HoleData.Exists and self._Gear__HoleData.type == "HEX":
            listaaguj = DOhexagono(self._Gear__vertices, self._Gear__faces, self._Gear__normals, 0, self.z, self.radiopoligono, 0, self.anchoeng, self._Gear__HoleData.circumradius)
            pass
        elif self._Gear__HoleData.Exists and self._Gear__HoleData.type == "SQR":
            listaaguj = DOcuadrado(self._Gear__vertices, self._Gear__faces, self._Gear__normals, 0, self.z, self.radiopoligono, 0, self.anchoeng, self._Gear__HoleData.circumradius)
            pass
        elif self._Gear__HoleData.Exists and self._Gear__HoleData.type == "CIR":
            listaaguj = DOcirculo(self._Gear__vertices, self._Gear__faces, self._Gear__normals, 0, self.z, self.radiopoligono, 0, self.anchoeng, self._Gear__HoleData.circumradius)
            pass
        elif self._Gear__HoleData.Exists and self._Gear__HoleData.type == "KW":
            #print("ES Keyway")
            listaaguj = DOkeyway(self._Gear__vertices, self._Gear__faces, self._Gear__normals, 0, self.z, self.radiopoligono, self.anchoeng, self._Gear__HoleData.boreDiameter, self._Gear__HoleData._alturaC, self._Gear__HoleData.keywayWidth)
            pass
        else:
            # el else solo aplica cuando es el agujero default o no lleva.
            pass

        self._Gear__vertices = listaaguj[0]
        self._Gear__faces = listaaguj[1]
        self._Gear__normals = listaaguj[2]


        print("Spur Gear calculations in %s seconds" % (time.time() - start_time2))

class HelicalGear(Gear):
    def __init__(self, modulo, apdeg, ahdeg, z, anchoeng, HoleData):
        super().__init__(modulo, apdeg, z, anchoeng, HoleData)
        self.ahdeg = ahdeg

        self.__ah = degToRad(ahdeg)

        # EN ESTE BLOQUE ESTÁ LA FUNCIÓN DE DIVISIÓN PARA LAS HÉLICES ORIGINAL
        # self.__seps es la variable que me indica la cantidad de self.__separaciones que tiene el engrane helicoidal (numero de engranes huecos)
        # divisor = 2
        # if 35<ahdeg>=40:
        #     divisor = 1.5
        # elif ahdeg>40:
        #     divisor = 1
        # esta es la separación para la hélice original
        # seps1 = int(lengthofhelix(self.__ah, self.anchoeng, self.__dp) / divisor)

        """
        NOTA MUY IMPORTANTE:
                Para arreglar el bug de que no agarra el agujero en los hel y dblhel cambie la precision de la helice. Si en algun punto vuelve a fallar,
                descomenta la solucion que ahuevo debe de jalar pero es menos eficiente.
        """

        """Si el angulo de rotación helicoidal es menor a tres el engrane se va a hacer mal (recuerda que se reservan self.__separaciones para las dos tapas)
        El número multiplicador para que siempre se vea en HD es 125
        Nótese que tuviste que poner un abs por la cuestión del mirror, esto porque seps1 a veces salía negativo por lo del anghrot
        (checa el comment de la parte del agujero para que entiendas mejor)"""
        seps1 = abs(int(helixturns(self.__ah, self.anchoeng, self._Gear__dp) * 100))
        # bug fixer seps1 = abs(int(helixturns(self.__ah, self.anchoeng, self._Gear__dp) * 500))
        #print("Juanasdasd", seps1)

        if seps1 >= 6:
            self.__seps = seps1
        else:
            self.__seps = 6

        self.__interdentsep = np.linspace((-self._Gear__alpha - self._Gear__beta), (self._Gear__alpha + self._Gear__beta), self._Gear__sep2)
        self.__sigm = sigma(self._Gear__rva, self._Gear__rb, self._Gear__dva, self._Gear__T, self._Gear__dv, self._Gear__ap)
        self.__sigsep = np.linspace(-self.__sigm / 2, self.__sigm / 2, self._Gear__sep3)

        self.__ph = np.pi * self._Gear__dv * (np.cos(self.__ah) / np.sin(self.__ah))
        self.__anchosep = self.anchoeng / self.__seps
        self.__anghrot = (self.anchoeng * 2 * np.pi) / self.__ph
        print("HELICOIDAL", self.__anghrot)
        self.__ahsep = self.__anghrot / self.__seps

    def GetFacets(self):
        #debido a la inserción de las cuatro involutas en un mismo for, pasado de cierto punto (que lo marca quality) se requiere insertarlas en un orden de magnitud menor (ver el cuaderno)
        #o sea que cuando i llegaba la mitad de las self.__separaciones empezaba a comerse caras que no eran de las cuatro involutas
        quality = (np.ceil(self._Gear__separaciones / 2) - 1)
        # region Declaración de medida para radiopolígono
        if self._Gear__HoleData.Exists:
            self.radiopoligono = self._Gear__rvf - .5
        else:
            self.radiopoligono = 0
        # endregion

        # region Este es el GETFACETS del helicoidal
        start_time2 = time.time()
        # region Declaración de vertices, caras, tapas e inserción de la tapa inferior en el engrane
        # Comienzo a crear el engrane
        self._Gear__vertices = np.array([[0, 0, self._Gear__a1]])
        self._Gear__normals = np.array([[0, 0, 0]])
        self._Gear__faces = np.array([[0, 0, 0]])

        # k es el residuo de haber utilizado un for para crear los dientes en la primer versión
        # k es innecesario, necesita refactorizarse y ser eliminado
        k = 1
        # Involuta parte hacia arriba (que en realidad va abajo)
        xinvar = self._Gear__rb * (
                np.cos(
                    self._Gear__tetarvasep + self._Gear__angrotp + k * self._Gear__az) + self._Gear__tetarvasep * np.sin(
            self._Gear__tetarvasep + self._Gear__angrotp + k * self._Gear__az))
        yinvar = self._Gear__rb * (
                np.sin(
                    self._Gear__tetarvasep + self._Gear__angrotp + k * self._Gear__az) - self._Gear__tetarvasep * np.cos(
            self._Gear__tetarvasep + self._Gear__angrotp + k * self._Gear__az))
        # Involuta parte hacia abajo (que en realidad va arriba)
        xinvab = self._Gear__rb * (
                np.cos(
                    self._Gear__tetarvasep + self._Gear__angrotp - k * self._Gear__az) + self._Gear__tetarvasep * np.sin(
            self._Gear__tetarvasep + self._Gear__angrotp - k * self._Gear__az))
        yinvab = -self._Gear__rb * (
                np.sin(
                    self._Gear__tetarvasep + self._Gear__angrotp - k * self._Gear__az) - self._Gear__tetarvasep * np.cos(
            self._Gear__tetarvasep + self._Gear__angrotp - k * self._Gear__az))
        # Arco de adendo para el diente
        xarc = self._Gear__rva * np.cos(self._Gear__sigsep + k * self._Gear__az)
        yarc = self._Gear__rva * np.sin(self._Gear__sigsep + k * self._Gear__az)
        # Arco de dedendo interdental
        xarcf = self._Gear__rvf * np.cos(self.__interdentsep + k * self._Gear__az)
        yarcf = self._Gear__rvf * np.sin(self.__interdentsep + k * self._Gear__az)
        # Creo los vertices para el arco contiguo
        xarco = self._Gear__rvf * np.cos(self._Gear__tsep + k * self._Gear__az)
        yarco = self._Gear__rvf * np.sin(self._Gear__tsep + k * self._Gear__az)

        rotavangulo = self.__ahsep
        rotavert = np.array(
            [[np.cos(rotavangulo), -np.sin(rotavangulo), 0], [np.sin(rotavangulo), np.cos(rotavangulo), 0],
             [0, 0, 1]])

        # SECCION self._Gear__vertices
        lastijuan = len((self._Gear__vertices))

        # Diente con listas para eliminar caras
        sincarasarriba = []
        SCabajo = []
        for i in range(0, len(xinvar)):
            self._Gear__vertices = np.append(self._Gear__vertices, [[xinvar[i], yinvar[i], self._Gear__a1]], axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices, [[xinvab[i], yinvab[i], self._Gear__a1]], axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices,
                                             np.matmul([[xinvar[i], yinvar[i], self.__anchosep]], rotavert), axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices,
                                             np.matmul([[xinvab[i], yinvab[i], self.__anchosep]], rotavert), axis=0)
            # j=i*4+1
            beta = 8 * i + lastijuan - 1
            # LA REGLA PARA i parece que es (separaciones/2) -1 donde la división se redondea hacia arriba entonces para separaciones=10 -> i = 4
            # y para separaciones=11 i = 5
            if i < quality:
                # Caras de abajo
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 2, beta + 5]], axis=0)
                SCabajo.append(len(self._Gear__faces) - 1)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 2, beta + 6, beta + 5]], axis=0)
                SCabajo.append(len(self._Gear__faces) - 1)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 5, beta + 6, beta + 9]], axis=0)
                SCabajo.append(len(self._Gear__faces) - 1)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 6, beta + 10, beta + 9]], axis=0)
                SCabajo.append(len(self._Gear__faces) - 1)
                # Caras de arriba
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 3, beta + 7, beta + 4]], axis=0)
                sincarasarriba.append(len(self._Gear__faces) - 1)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 4, beta + 7, beta + 8]], axis=0)
                sincarasarriba.append(len(self._Gear__faces) - 1)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 7, beta + 11, beta + 8]], axis=0)
                sincarasarriba.append(len(self._Gear__faces) - 1)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 8, beta + 11, beta + 12]], axis=0)
                sincarasarriba.append(len(self._Gear__faces) - 1)
                # Laterales de izquierda (viendo de frente)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 5, beta + 3]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 3, beta + 5, beta + 7]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 5, beta + 9, beta + 7]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 7, beta + 9, beta + 11]], axis=0)
                # Laterales de derecha (viendo de frente)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 2, beta + 4, beta + 6]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 4, beta + 8, beta + 6]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 6, beta + 8, beta + 10]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 8, beta + 12, beta + 10]], axis=0)
            elif i == quality:
                # Caras de abajo
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 2, beta + 5]], axis=0)
                SCabajo.append(len(self._Gear__faces) - 1)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 2, beta + 6, beta + 5]], axis=0)
                SCabajo.append(len(self._Gear__faces) - 1)
                # Caras de arriba
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 3, beta + 7, beta + 4]], axis=0)
                sincarasarriba.append(len(self._Gear__faces) - 1)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 4, beta + 7, beta + 8]], axis=0)
                sincarasarriba.append(len(self._Gear__faces) - 1)
                # Laterales de izquierda (viendo de frente)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 5, beta + 3]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 3, beta + 5, beta + 7]], axis=0)
                # Laterales de derecha (viendo de frente)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 2, beta + 4, beta + 6]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 4, beta + 8, beta + 6]], axis=0)
                pass

        # Arco de adendo
        self._Gear__vertices = np.append(self._Gear__vertices, [
            [self._Gear__puntor * np.cos(k * self._Gear__az), self._Gear__puntor * np.sin(k * self._Gear__az),
             self._Gear__a1]], axis=0)
        self._Gear__vertices = np.append(self._Gear__vertices,
                                         np.matmul([[self._Gear__puntor * np.cos(k * self._Gear__az),
                                                     self._Gear__puntor * np.sin(k * self._Gear__az),
                                                     self.__anchosep]],
                                                   rotavert),
                                         axis=0)
        lastijuan = len(self._Gear__vertices)
        for i in range(0, len(xarc)):
            self._Gear__vertices = np.append(self._Gear__vertices, [[xarc[i], yarc[i], self._Gear__a1]], axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices,
                                             np.matmul([[xarc[i], yarc[i], self.__anchosep]], rotavert), axis=0)
            if i < len(xarc) - 1:
                beta = i * 2 + lastijuan
                # Caras de abajo
                self._Gear__faces = np.append(self._Gear__faces, [[lastijuan - 2, beta + 2, beta]], axis=0)
                SCabajo.append(len(self._Gear__faces) - 1)
                # Caras de arriba
                self._Gear__faces = np.append(self._Gear__faces, [[lastijuan - 1, beta + 1, beta + 3]], axis=0)
                sincarasarriba.append(len(self._Gear__faces) - 1)
                # Caras frontales con la base arriba
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 2, beta + 3]], axis=0)
                # Caras frontales con base abajo
                self._Gear__faces = np.append(self._Gear__faces, [[beta, beta + 2, beta + 1]], axis=0)

        lastijuan = len(self._Gear__vertices)
        for i in range(0, len(xarco)):
            self._Gear__vertices = np.append(self._Gear__vertices, [[xarco[i], yarco[i], self._Gear__a1]], axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices,
                                             np.matmul([[xarco[i], yarco[i], self.__anchosep]], rotavert), axis=0)
            if i < len(xarc) - 1:
                beta = i * 2 + lastijuan
                # Caras de abajo
                self._Gear__faces = np.append(self._Gear__faces, [[lastijuan + 2 * len(xarco), beta + 2, beta]],
                                              axis=0)
                SCabajo.append(len(self._Gear__faces) - 1)
                # Caras de arriba
                self._Gear__faces = np.append(self._Gear__faces,
                                              [[lastijuan + 2 * len(xarco) + 1, beta + 1, beta + 3]], axis=0)
                sincarasarriba.append(len(self._Gear__faces) - 1)
                # Caras frontales con la base arriba
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 2, beta + 3]], axis=0)
                # Caras frontales con base abajo
                self._Gear__faces = np.append(self._Gear__faces, [[beta, beta + 2, beta + 1]], axis=0)

        # Meto los centros para crear el polígono exterior (debe haber un if por si el radio es 0)
        self._Gear__vertices = np.append(self._Gear__vertices, [
            [(self.radiopoligono) * np.cos(k * self._Gear__az), (self.radiopoligono) * np.sin(k * self._Gear__az),
             self._Gear__a1]], axis=0)
        self._Gear__vertices = np.append(self._Gear__vertices,
                                         np.matmul([[(self.radiopoligono) * np.cos(k * self._Gear__az),
                                                     (self.radiopoligono) * np.sin(k * self._Gear__az),
                                                     self.__anchosep]],
                                                   rotavert), axis=0)

        # Creo el triangulo del origen rpol a las involutas
        self._Gear__faces = np.append(self._Gear__faces,
                                      [[lastijuan + 2 * len(xarco), lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 1,
                                        lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 2]], axis=0)
        SCabajo.append(len(self._Gear__faces) - 1)
        self._Gear__faces = np.append(self._Gear__faces,
                                      [[lastijuan + 2 * len(xarco) + 1, lastijuan - 2 * len(xarc) - 4 * len(xinvar),
                                        lastijuan - 2 * len(xarc) - 4 * len(xinvar) + 1]], axis=0)
        sincarasarriba.append(len(self._Gear__faces) - 1)

        lastijuan2 = len(self._Gear__vertices)
        # Selladores de los arcos contiguos con los origenes del self.radiopoligono
        # Puntos del  siguiente arco contiguo
        self._Gear__vertices = np.append(self._Gear__vertices, [[xarcf[0], yarcf[0], self._Gear__a1]], axis=0)
        self._Gear__vertices = np.append(self._Gear__vertices,
                                         np.matmul([[xarcf[0], yarcf[0], self.__anchosep]], rotavert), axis=0)

        if self._Gear__rvf < self._Gear__rb:
            # Selladores del radio de raiz a la primera involuta
            # Caras laterales
            # Derecha base abajo
            self._Gear__faces = np.append(self._Gear__faces,
                                          [[lastijuan, lastijuan - 2 * len(xarc) - 4 * len(xinvar) + 1,
                                            lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 1]], axis=0)
            # Derecha base arriba
            self._Gear__faces = np.append(self._Gear__faces,
                                          [[lastijuan + 1, lastijuan - 2 * len(xarc) - 4 * len(xinvar) + 1,
                                            lastijuan]],
                                          axis=0)
            # Izquierda
            # Base abajo
            self._Gear__faces = np.append(self._Gear__faces, [
                [lastijuan2, lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 2,
                 lastijuan - 2 * len(xarc) - 4 * len(xinvar)]],
                                          axis=0)
            # Base arriba
            self._Gear__faces = np.append(self._Gear__faces,
                                          [[lastijuan2, lastijuan - 2 * len(xarc) - 4 * len(xinvar),
                                            lastijuan2 + 1]],
                                          axis=0)

        if self.radiopoligono != 0 and self.radiopoligono != -1:
            # Selladores del radio de raiz a la primera involuta
            # Caras de abajo
            self._Gear__faces = np.append(self._Gear__faces, [
                [lastijuan + 2 * len(xarco), lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 2, lastijuan2]], axis=0)
            SCabajo.append(len(self._Gear__faces) - 1)
            # Derecha
            self._Gear__faces = np.append(self._Gear__faces,
                                          [[lastijuan + 2 * len(xarco), lastijuan,
                                            lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 1]],
                                          axis=0)
            SCabajo.append(len(self._Gear__faces) - 1)

            # Caras de arriba
            # Izquierda
            self._Gear__faces = np.append(self._Gear__faces, [
                [lastijuan + 2 * len(xarco) + 1, lastijuan2 + 1, lastijuan - 2 * len(xarc) - 4 * len(xinvar)]],
                                          axis=0)
            sincarasarriba.append(len(self._Gear__faces) - 1)
            # Derecha
            self._Gear__faces = np.append(self._Gear__faces, [
                [lastijuan + 2 * len(xarco) + 1, lastijuan - 2 * len(xarc) - 4 * len(xinvar) + 1, lastijuan + 1]],
                                          axis=0)
            sincarasarriba.append(len(self._Gear__faces) - 1)

            # Creo los complementos triangulares para los arcos contiguos que siguen
            self._Gear__vertices = np.append(self._Gear__vertices,
                                             [[(self.radiopoligono) * np.cos(k * self._Gear__az + self._Gear__az),
                                               (self.radiopoligono) * np.sin(k * self._Gear__az + self._Gear__az),
                                               self._Gear__a1]],
                                             axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices, np.matmul(
                [[(self.radiopoligono) * np.cos(k * self._Gear__az + self._Gear__az),
                  (self.radiopoligono) * np.sin(k * self._Gear__az + self._Gear__az), self.__anchosep]], rotavert),
                                             axis=0)

            # Triangulotes entre centros rpol
            # El de abajo
            self._Gear__faces = np.append(self._Gear__faces,
                                          [[lastijuan + 2 * len(xarco), lastijuan2 + 2, lastijuan2 - 4]], axis=0)
            SCabajo.append(len(self._Gear__faces) - 1)
            # El de arriba
            self._Gear__faces = np.append(self._Gear__faces,
                                          [[lastijuan + 2 * len(xarco) + 1, lastijuan2 - 3, lastijuan2 + 3]],
                                          axis=0)
            sincarasarriba.append(len(self._Gear__faces) - 1)

        # Appendo la primera cara(0) porque es inutil (unía el origen con el origen con el origen)
        sincarasarriba.append(0)

        # Creo SCplanas antes de que SCabajo appende el 0 para evitar que se elimine una cara extra que si es relevante
        SCplanas = []
        SCplanas.extend(sincarasarriba)
        SCplanas.extend(SCabajo)

        SCabajo.append(0)

        faceSCarriba = np.delete(self._Gear__faces, sincarasarriba, 0)
        faceSCabajo = np.delete(self._Gear__faces, SCabajo, 0)
        faceSCplanas = np.delete(self._Gear__faces, SCplanas, 0)
        self._Gear__faces = faceSCarriba
        facesceterno = faceSCplanas
        # Esta parte es crucial para los insertos de dientes para que se introduzcan en una matriz temporal y luego en la de verdad
        # (es meramente demostrativo porque lo que meta en la función es lo que importa)
        facejuan3 = self._Gear__faces
        vertijuan3 = self._Gear__vertices
        vertieterno = self._Gear__vertices
        faceterno = self._Gear__faces

        # endregion

        lista = conversoraBinario(self.z)
        liston = conversoraBinario(self.__seps - 2)

        # Creo las partes del diente sin tapas
        dientes = insertodientes12(liston, self.__ahsep, vertieterno, faceSCplanas, self.__anchosep)
        vertsinplanas = dientes[0]
        facesinplanas = dientes[1]

        dientes = insertodientes22(liston, self.__ahsep, vertsinplanas, facesinplanas, vertieterno, faceSCplanas,
                                   vertieterno,
                                   faceSCplanas, self.__anchosep)
        vertsinplanas = dientes[0]
        facesinplanas = dientes[1]

        # Inserto la sección media a la tapa de abajo (recuerda que la tapa de abajo ya está puesta por default en self._Gear__vertices)
        vertsinplanas = vertsinplanas + np.array([0, 0, self.__anchosep])
        vertsinplanas = rotateZ(vertsinplanas, self.__ahsep)
        zum = len(vertieterno)
        facesinplanas = facesinplanas + np.array([zum, zum, zum])
        self._Gear__vertices = np.append(self._Gear__vertices, vertsinplanas, axis=0)
        self._Gear__faces = np.append(self._Gear__faces, facesinplanas, axis=0)

        # Inserto la tapa de arriba
        vertarriba = vertieterno + np.array([0, 0, self.__anchosep * (self.__seps - 1)])
        vertarriba = rotateZ(vertarriba, self.__ahsep * (self.__seps - 1))
        zum = len(self._Gear__vertices)
        facearriba = faceSCabajo + np.array([zum, zum, zum])

        self._Gear__vertices = np.append(self._Gear__vertices, vertarriba, axis=0)
        self._Gear__faces = np.append(self._Gear__faces, facearriba, axis=0)
        # Aquí ya tengo un diente completo

        for i in range(0, len(self._Gear__faces)):
            we = self._Gear__faces[i]
            v1 = self._Gear__vertices[we[0]]
            v2 = self._Gear__vertices[we[1]]
            v3 = self._Gear__vertices[we[2]]
            normal = normal_to_surface(v1, v2, v3)
            self._Gear__normals = np.append(self._Gear__normals, [normal], axis=0)

        self._Gear__normals = np.delete(self._Gear__normals, 0, 0)
        norfinito = self._Gear__normals
        # Conversor a binerio es extremadamente importante. Crea una lista con los exponentes para 2**x de manera que controla la cantidad de dientes
        lista = conversoraBinario(self.z)

        self._Gear__normals = normales1(lista, self._Gear__az, self._Gear__normals)

        # Inserto el resto de los dientes
        vertieterno = self._Gear__vertices
        faceterno = self._Gear__faces
        facejuan3 = self._Gear__faces
        vertijuan3 = self._Gear__vertices
        lista = conversoraBinario(self.z)

        # Creo el inserto del primer grupo de dientes (recuerda que insertodientes1 afectamente la matriz final mientras que insertodientes2 no lo hace hasta el final
        dientes = insertodientes1(lista, self._Gear__az, self._Gear__vertices, self._Gear__faces, vertijuan3,
                                  facejuan3, vertieterno,
                                  faceterno)
        self._Gear__vertices = dientes[0]
        self._Gear__faces = dientes[1]
        self._Gear__normals = normales2(lista, self._Gear__az, self._Gear__normals, norfinito, norfinito)

        # Creo los insertos para el resto de los dientes
        dientes = insertodientes2(lista, self._Gear__az, self._Gear__vertices, self._Gear__faces, vertijuan3,
                                  facejuan3, vertieterno,
                                  faceterno)
        self._Gear__vertices = dientes[0]
        self._Gear__faces = dientes[1]

        # endregion

        #For some reason I have the helix mirrored, which means that I have to use the negative of the anghrot. This doesn't mean it's wrong, it's just rotating in an oposite direction
        #Hole = helhexagono(self._Gear__vertices, self._Gear__faces, self._Gear__normals, self.z, self.radiopoligono, 4, 0, self.anchoeng, -self.__anghrot)
        #Hole = helkeyway(self._Gear__vertices, self._Gear__faces, self._Gear__normals,self.z,self.radiopoligono,self.anchoeng,3,2,2,-self.__anghrot)
        listaaguj = [self._Gear__vertices, self._Gear__faces, self._Gear__normals]
        # IMPORTANTE !! ESTO JALARIA COMO FUNCION O METODO porque DH Y DR comparten algoritmos DO solo HL es diferente
        # NOTA: HAZ QUE MATCHEEN CON LO QUE ARROJA TU UI -> HEX, SQR, CIR, KW, None
        if self._Gear__HoleData.Exists and self._Gear__HoleData.type == "HEX":
            listaaguj = helhexagono(self._Gear__vertices, self._Gear__faces, self._Gear__normals, self.z,
                                    self.radiopoligono, self._Gear__HoleData.circumradius, 0, self.anchoeng,
                                    -self.__anghrot)
            pass
        elif self._Gear__HoleData.Exists and self._Gear__HoleData.type == "SQR":
            listaaguj = helcuadrado(self._Gear__vertices, self._Gear__faces, self._Gear__normals, self.z,
                                    self.radiopoligono, self._Gear__HoleData.circumradius, 0, self.anchoeng,
                                    -self.__anghrot)
            pass
        elif self._Gear__HoleData.Exists and self._Gear__HoleData.type == "CIR":
            listaaguj = helcirculo(self._Gear__vertices, self._Gear__faces, self._Gear__normals, self.z,
                                   self.radiopoligono, 0, self.anchoeng, self._Gear__HoleData.circumradius,
                                   -self.__anghrot)
            pass
        elif self._Gear__HoleData.Exists and self._Gear__HoleData.type == "KW":
            # print("ES Keyway")
            listaaguj = helkeyway(self._Gear__vertices, self._Gear__faces, self._Gear__normals, self.z,
                                  self.radiopoligono, self.anchoeng, self._Gear__HoleData.boreDiameter,
                                  self._Gear__HoleData._alturaC, self._Gear__HoleData.keywayWidth, -self.__anghrot)
            #listaaguj =  DOkeyway(self._Gear__vertices, self._Gear__faces, self._Gear__normals, self.__anghrot, self.z, self.radiopoligono, self.anchoeng, self._Gear__HoleData.boreDiameter, self._Gear__HoleData._alturaC, self._Gear__HoleData.keywayWidth)

            pass
        else:
            # el else solo aplica cuando es el agujero default o no lleva.
            pass

        self._Gear__vertices = listaaguj[0]
        self._Gear__faces = listaaguj[1]
        self._Gear__normals = listaaguj[2]
        print("AAAAAAA",2*np.pi/self.__anghrot)
        #omega = len(vertijuan)
        # print(len(vertijuan)/omega)
        # Imprime en la pantalla la cantidad de dientes que tiene el engrane realmente
        # print("zzz", len(self._Gear__vertices) / len(vertieterno), len(self._Gear__faces ))
        #print("Optimizado3 en %s seconds" % (time.time() - start_time2))
        # endregion

class DoubleHelicalGear(Gear):
    def __init__(self, modulo, apdeg, ahdeg, z, anchoeng, HoleData):
        super().__init__(modulo, apdeg, z, anchoeng, HoleData)
        self.anchoeng = self.anchoeng/2
        self.ahdeg = ahdeg

        self.__ah = degToRad(ahdeg)

        # EN ESTE BLOQUE ESTÁ LA FUNCIÓN DE DIVISIÓN PARA LAS HÉLICES ORIGINAL
        # self.__seps es la variable que me indica la cantidad de self.__separaciones que tiene el engrane helicoidal (numero de engranes huecos)
        # divisor = 2
        # if 35<ahdeg>=40:
        #     divisor = 1.5
        # elif ahdeg>40:
        #     divisor = 1
        # esta es la separación para la hélice original
        # seps1 = int(lengthofhelix(self.__ah, self.anchoeng, self.__dp) / divisor)

        """
            NOTA MUY IMPORTANTE:
                        Para arreglar el bug de que no agarra el agujero en los hel y dblhel cambie la precision de la helice. Si en algun punto vuelve a fallar,
                        descomenta la solucion que ahuevo debe de jalar pero es menos eficiente.
        """

        # Si el angulo de rotación helicoidal es menor a tres el engrane se va a hacer mal (recuerda que se reservan self.__separaciones para las dos tapas)

        """ 
        Este seps1 es el OG de dblhelicoidal, pero lo cambie por lo del bug que en el helicoidal no detecta el agujero t0do el tiemposeps1 = int(helixturns(self.__ah, self.anchoeng, self._Gear__dp) * 125)
        """
        seps1 = abs(int(helixturns(self.__ah, self.anchoeng, self._Gear__dp) * 100))
        # bug fixer for hole detection of helical gears seps1 = abs(int(helixturns(self.__ah, self.anchoeng, self._Gear__dp) * 500))
        if seps1 >= 6:
            self.__seps = seps1
        else:
            self.__seps = 6

        self.__interdentsep = np.linspace((-self._Gear__alpha - self._Gear__beta), (self._Gear__alpha + self._Gear__beta),
                                          self._Gear__sep2)
        self.__sigm = sigma(self._Gear__rva, self._Gear__rb, self._Gear__dva, self._Gear__T, self._Gear__dv, self._Gear__ap)
        self.__sigsep = np.linspace(-self.__sigm / 2, self.__sigm / 2, self._Gear__sep3)

        self.__ph = np.pi * self._Gear__dv * (np.cos(self.__ah) / np.sin(self.__ah))
        self.__anchosep = self.anchoeng / self.__seps
        self.__anghrot = (self.anchoeng * 2 * np.pi) / self.__ph
        print("DHELICOIDAL", self.__anghrot)
        self.__ahsep = self.__anghrot / self.__seps

    def GetFacets(self):
        # debido a la inserción de las cuatro involutas en un mismo for, pasado de cierto punto (que lo marca quality) se requiere insertarlas en un orden de magnitud menor (ver el cuaderno)
        # o sea que cuando i llegaba la mitad de las self.__separaciones empezaba a comerse caras que no eran de las cuatro involutas
        quality = (np.ceil(self._Gear__separaciones / 2) - 1)
        # region Declaración de medida para radiopolígono
        if self._Gear__HoleData.Exists:
            self.radiopoligono = self._Gear__rvf - .5
        else:
            self.radiopoligono = 0
        # endregion

        # region Este es el GETFACETS del helicoidal
        start_time2 = time.time()
        # region Declaración de vertices, caras, tapas e inserción de la tapa inferior en el engrane
        # Comienzo a crear el engrane
        self._Gear__vertices = np.array([[0, 0, self._Gear__a1]])
        self._Gear__normals = np.array([[0, 0, 0]])
        self._Gear__faces = np.array([[0, 0, 0]])

        # k es el residuo de haber utilizado un for para crear los dientes en la primer versión
        # k es innecesario, necesita refactorizarse y ser eliminado
        k = 1
        # Involuta parte hacia arriba (que en realidad va abajo)
        xinvar = self._Gear__rb * (
                np.cos(
                    self._Gear__tetarvasep + self._Gear__angrotp + k * self._Gear__az) + self._Gear__tetarvasep * np.sin(
            self._Gear__tetarvasep + self._Gear__angrotp + k * self._Gear__az))
        yinvar = self._Gear__rb * (
                np.sin(
                    self._Gear__tetarvasep + self._Gear__angrotp + k * self._Gear__az) - self._Gear__tetarvasep * np.cos(
            self._Gear__tetarvasep + self._Gear__angrotp + k * self._Gear__az))
        # Involuta parte hacia abajo (que en realidad va arriba)
        xinvab = self._Gear__rb * (
                np.cos(
                    self._Gear__tetarvasep + self._Gear__angrotp - k * self._Gear__az) + self._Gear__tetarvasep * np.sin(
            self._Gear__tetarvasep + self._Gear__angrotp - k * self._Gear__az))
        yinvab = -self._Gear__rb * (
                np.sin(
                    self._Gear__tetarvasep + self._Gear__angrotp - k * self._Gear__az) - self._Gear__tetarvasep * np.cos(
            self._Gear__tetarvasep + self._Gear__angrotp - k * self._Gear__az))
        # Arco de adendo para el diente
        xarc = self._Gear__rva * np.cos(self._Gear__sigsep + k * self._Gear__az)
        yarc = self._Gear__rva * np.sin(self._Gear__sigsep + k * self._Gear__az)
        # Arco de dedendo interdental
        xarcf = self._Gear__rvf * np.cos(self.__interdentsep + k * self._Gear__az)
        yarcf = self._Gear__rvf * np.sin(self.__interdentsep + k * self._Gear__az)
        # Creo los vertices para el arco contiguo
        xarco = self._Gear__rvf * np.cos(self._Gear__tsep + k * self._Gear__az)
        yarco = self._Gear__rvf * np.sin(self._Gear__tsep + k * self._Gear__az)

        rotavangulo = self.__ahsep
        rotavert = np.array(
            [[np.cos(rotavangulo), -np.sin(rotavangulo), 0], [np.sin(rotavangulo), np.cos(rotavangulo), 0],
             [0, 0, 1]])

        # SECCION self._Gear__vertices
        lastijuan = len((self._Gear__vertices))

        # Diente con listas para eliminar caras
        sincarasarriba = []
        SCabajo = []
        for i in range(0, len(xinvar)):
            self._Gear__vertices = np.append(self._Gear__vertices, [[xinvar[i], yinvar[i], self._Gear__a1]], axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices, [[xinvab[i], yinvab[i], self._Gear__a1]], axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices,
                                             np.matmul([[xinvar[i], yinvar[i], self.__anchosep]], rotavert), axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices,
                                             np.matmul([[xinvab[i], yinvab[i], self.__anchosep]], rotavert), axis=0)
            # j=i*4+1
            beta = 8 * i + lastijuan - 1
            # LA REGLA PARA i parece que es (separaciones/2) -1 donde la división se redondea hacia arriba entonces para separaciones=10 -> i = 4
            # y para separaciones=11 i = 5
            if i < quality:
                # Caras de abajo
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 2, beta + 5]], axis=0)
                SCabajo.append(len(self._Gear__faces) - 1)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 2, beta + 6, beta + 5]], axis=0)
                SCabajo.append(len(self._Gear__faces) - 1)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 5, beta + 6, beta + 9]], axis=0)
                SCabajo.append(len(self._Gear__faces) - 1)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 6, beta + 10, beta + 9]], axis=0)
                SCabajo.append(len(self._Gear__faces) - 1)
                # Caras de arriba
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 3, beta + 7, beta + 4]], axis=0)
                sincarasarriba.append(len(self._Gear__faces) - 1)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 4, beta + 7, beta + 8]], axis=0)
                sincarasarriba.append(len(self._Gear__faces) - 1)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 7, beta + 11, beta + 8]], axis=0)
                sincarasarriba.append(len(self._Gear__faces) - 1)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 8, beta + 11, beta + 12]], axis=0)
                sincarasarriba.append(len(self._Gear__faces) - 1)
                # Laterales de izquierda (viendo de frente)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 5, beta + 3]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 3, beta + 5, beta + 7]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 5, beta + 9, beta + 7]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 7, beta + 9, beta + 11]], axis=0)
                # Laterales de derecha (viendo de frente)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 2, beta + 4, beta + 6]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 4, beta + 8, beta + 6]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 6, beta + 8, beta + 10]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 8, beta + 12, beta + 10]], axis=0)
            elif i == quality:
                # Caras de abajo
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 2, beta + 5]], axis=0)
                SCabajo.append(len(self._Gear__faces) - 1)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 2, beta + 6, beta + 5]], axis=0)
                SCabajo.append(len(self._Gear__faces) - 1)
                # Caras de arriba
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 3, beta + 7, beta + 4]], axis=0)
                sincarasarriba.append(len(self._Gear__faces) - 1)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 4, beta + 7, beta + 8]], axis=0)
                sincarasarriba.append(len(self._Gear__faces) - 1)
                # Laterales de izquierda (viendo de frente)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 5, beta + 3]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 3, beta + 5, beta + 7]], axis=0)
                # Laterales de derecha (viendo de frente)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 2, beta + 4, beta + 6]], axis=0)
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 4, beta + 8, beta + 6]], axis=0)
                pass

        # Arco de adendo
        self._Gear__vertices = np.append(self._Gear__vertices, [
            [self._Gear__puntor * np.cos(k * self._Gear__az), self._Gear__puntor * np.sin(k * self._Gear__az),
             self._Gear__a1]], axis=0)
        self._Gear__vertices = np.append(self._Gear__vertices,
                                         np.matmul([[self._Gear__puntor * np.cos(k * self._Gear__az),
                                                     self._Gear__puntor * np.sin(k * self._Gear__az),
                                                     self.__anchosep]],
                                                   rotavert),
                                         axis=0)
        lastijuan = len(self._Gear__vertices)
        for i in range(0, len(xarc)):
            self._Gear__vertices = np.append(self._Gear__vertices, [[xarc[i], yarc[i], self._Gear__a1]], axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices,
                                             np.matmul([[xarc[i], yarc[i], self.__anchosep]], rotavert), axis=0)
            if i < len(xarc) - 1:
                beta = i * 2 + lastijuan
                # Caras de abajo
                self._Gear__faces = np.append(self._Gear__faces, [[lastijuan - 2, beta + 2, beta]], axis=0)
                SCabajo.append(len(self._Gear__faces) - 1)
                # Caras de arriba
                self._Gear__faces = np.append(self._Gear__faces, [[lastijuan - 1, beta + 1, beta + 3]], axis=0)
                sincarasarriba.append(len(self._Gear__faces) - 1)
                # Caras frontales con la base arriba
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 2, beta + 3]], axis=0)
                # Caras frontales con base abajo
                self._Gear__faces = np.append(self._Gear__faces, [[beta, beta + 2, beta + 1]], axis=0)

        lastijuan = len(self._Gear__vertices)
        for i in range(0, len(xarco)):
            self._Gear__vertices = np.append(self._Gear__vertices, [[xarco[i], yarco[i], self._Gear__a1]], axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices,
                                             np.matmul([[xarco[i], yarco[i], self.__anchosep]], rotavert), axis=0)
            if i < len(xarc) - 1:
                beta = i * 2 + lastijuan
                # Caras de abajo
                self._Gear__faces = np.append(self._Gear__faces, [[lastijuan + 2 * len(xarco), beta + 2, beta]],
                                              axis=0)
                SCabajo.append(len(self._Gear__faces) - 1)
                # Caras de arriba
                self._Gear__faces = np.append(self._Gear__faces,
                                              [[lastijuan + 2 * len(xarco) + 1, beta + 1, beta + 3]], axis=0)
                sincarasarriba.append(len(self._Gear__faces) - 1)
                # Caras frontales con la base arriba
                self._Gear__faces = np.append(self._Gear__faces, [[beta + 1, beta + 2, beta + 3]], axis=0)
                # Caras frontales con base abajo
                self._Gear__faces = np.append(self._Gear__faces, [[beta, beta + 2, beta + 1]], axis=0)

        # Meto los centros para crear el polígono exterior (debe haber un if por si el radio es 0)
        self._Gear__vertices = np.append(self._Gear__vertices, [
            [(self.radiopoligono) * np.cos(k * self._Gear__az), (self.radiopoligono) * np.sin(k * self._Gear__az),
             self._Gear__a1]], axis=0)
        self._Gear__vertices = np.append(self._Gear__vertices,
                                         np.matmul([[(self.radiopoligono) * np.cos(k * self._Gear__az),
                                                     (self.radiopoligono) * np.sin(k * self._Gear__az),
                                                     self.__anchosep]],
                                                   rotavert), axis=0)

        # Creo el triangulo del origen rpol a las involutas
        self._Gear__faces = np.append(self._Gear__faces,
                                      [[lastijuan + 2 * len(xarco), lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 1,
                                        lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 2]], axis=0)
        SCabajo.append(len(self._Gear__faces) - 1)
        self._Gear__faces = np.append(self._Gear__faces,
                                      [[lastijuan + 2 * len(xarco) + 1, lastijuan - 2 * len(xarc) - 4 * len(xinvar),
                                        lastijuan - 2 * len(xarc) - 4 * len(xinvar) + 1]], axis=0)
        sincarasarriba.append(len(self._Gear__faces) - 1)

        lastijuan2 = len(self._Gear__vertices)
        # Selladores de los arcos contiguos con los origenes del self.radiopoligono
        # Puntos del  siguiente arco contiguo
        self._Gear__vertices = np.append(self._Gear__vertices, [[xarcf[0], yarcf[0], self._Gear__a1]], axis=0)
        self._Gear__vertices = np.append(self._Gear__vertices,
                                         np.matmul([[xarcf[0], yarcf[0], self.__anchosep]], rotavert), axis=0)

        if self._Gear__rvf < self._Gear__rb:
            # Selladores del radio de raiz a la primera involuta
            # Caras laterales
            # Derecha base abajo
            self._Gear__faces = np.append(self._Gear__faces,
                                          [[lastijuan, lastijuan - 2 * len(xarc) - 4 * len(xinvar) + 1,
                                            lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 1]], axis=0)
            # Derecha base arriba
            self._Gear__faces = np.append(self._Gear__faces,
                                          [[lastijuan + 1, lastijuan - 2 * len(xarc) - 4 * len(xinvar) + 1,
                                            lastijuan]],
                                          axis=0)
            # Izquierda
            # Base abajo
            self._Gear__faces = np.append(self._Gear__faces, [
                [lastijuan2, lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 2,
                 lastijuan - 2 * len(xarc) - 4 * len(xinvar)]],
                                          axis=0)
            # Base arriba
            self._Gear__faces = np.append(self._Gear__faces,
                                          [[lastijuan2, lastijuan - 2 * len(xarc) - 4 * len(xinvar),
                                            lastijuan2 + 1]],
                                          axis=0)

        if self.radiopoligono != 0 and self.radiopoligono != -1:
            # Selladores del radio de raiz a la primera involuta
            # Caras de abajo
            self._Gear__faces = np.append(self._Gear__faces, [
                [lastijuan + 2 * len(xarco), lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 2, lastijuan2]], axis=0)
            SCabajo.append(len(self._Gear__faces) - 1)
            # Derecha
            self._Gear__faces = np.append(self._Gear__faces,
                                          [[lastijuan + 2 * len(xarco), lastijuan,
                                            lastijuan - 2 * len(xarc) - 4 * len(xinvar) - 1]],
                                          axis=0)
            SCabajo.append(len(self._Gear__faces) - 1)

            # Caras de arriba
            # Izquierda
            self._Gear__faces = np.append(self._Gear__faces, [
                [lastijuan + 2 * len(xarco) + 1, lastijuan2 + 1, lastijuan - 2 * len(xarc) - 4 * len(xinvar)]],
                                          axis=0)
            sincarasarriba.append(len(self._Gear__faces) - 1)
            # Derecha
            self._Gear__faces = np.append(self._Gear__faces, [
                [lastijuan + 2 * len(xarco) + 1, lastijuan - 2 * len(xarc) - 4 * len(xinvar) + 1, lastijuan + 1]],
                                          axis=0)
            sincarasarriba.append(len(self._Gear__faces) - 1)

            # Creo los complementos triangulares para los arcos contiguos que siguen
            self._Gear__vertices = np.append(self._Gear__vertices,
                                             [[(self.radiopoligono) * np.cos(k * self._Gear__az + self._Gear__az),
                                               (self.radiopoligono) * np.sin(k * self._Gear__az + self._Gear__az),
                                               self._Gear__a1]],
                                             axis=0)
            self._Gear__vertices = np.append(self._Gear__vertices, np.matmul(
                [[(self.radiopoligono) * np.cos(k * self._Gear__az + self._Gear__az),
                  (self.radiopoligono) * np.sin(k * self._Gear__az + self._Gear__az), self.__anchosep]], rotavert),
                                             axis=0)

            # Triangulotes entre centros rpol
            # El de abajo
            self._Gear__faces = np.append(self._Gear__faces,
                                          [[lastijuan + 2 * len(xarco), lastijuan2 + 2, lastijuan2 - 4]], axis=0)
            SCabajo.append(len(self._Gear__faces) - 1)
            # El de arriba
            self._Gear__faces = np.append(self._Gear__faces,
                                          [[lastijuan + 2 * len(xarco) + 1, lastijuan2 - 3, lastijuan2 + 3]],
                                          axis=0)
            sincarasarriba.append(len(self._Gear__faces) - 1)

        # Appendo la primera cara(0) porque es inutil (unía el origen con el origen con el origen)
        sincarasarriba.append(0)

        # Creo SCplanas antes de que SCabajo appende el 0 para evitar que se elimine una cara extra que si es relevante
        SCplanas = []
        SCplanas.extend(sincarasarriba)
        SCplanas.extend(SCabajo)

        SCabajo.append(0)

        faceSCarriba = np.delete(self._Gear__faces, sincarasarriba, 0)
        faceSCabajo = np.delete(self._Gear__faces, SCabajo, 0)
        faceSCplanas = np.delete(self._Gear__faces, SCplanas, 0)
        # Truco para que las caras de abajo no estén
        self._Gear__faces = faceSCplanas
        facesceterno = faceSCplanas
        # Esta parte es crucial para los insertos de dientes para que se introduzcan en una matriz temporal y luego en la de verdad
        # (es meramente demostrativo porque lo que meta en la función es lo que importa)
        facejuan3 = self._Gear__faces
        vertijuan3 = self._Gear__vertices
        vertieterno = self._Gear__vertices
        faceterno = self._Gear__faces

        # endregion

        lista = conversoraBinario(self.z)
        liston = conversoraBinario(self.__seps - 2)

        # Creo las partes del diente sin tapas
        dientes = insertodientes12(liston, self.__ahsep, vertieterno, faceSCplanas, self.__anchosep)
        vertsinplanas = dientes[0]
        facesinplanas = dientes[1]

        dientes = insertodientes22(liston, self.__ahsep, vertsinplanas, facesinplanas, vertieterno, faceSCplanas,
                                   vertieterno,
                                   faceSCplanas, self.__anchosep)
        vertsinplanas = dientes[0]
        facesinplanas = dientes[1]

        # Inserto la sección media a la tapa de abajo (recuerda que la tapa de abajo ya está puesta por default en self._Gear__vertices)
        vertsinplanas = vertsinplanas + np.array([0, 0, self.__anchosep])
        vertsinplanas = rotateZ(vertsinplanas, self.__ahsep)
        zum = len(vertieterno)
        facesinplanas = facesinplanas + np.array([zum, zum, zum])
        self._Gear__vertices = np.append(self._Gear__vertices, vertsinplanas, axis=0)
        self._Gear__faces = np.append(self._Gear__faces, facesinplanas, axis=0)

        # Inserto la tapa de arriba
        vertarriba = vertieterno + np.array([0, 0, self.__anchosep * (self.__seps - 1)])
        vertarriba = rotateZ(vertarriba, self.__ahsep * (self.__seps - 1))
        zum = len(self._Gear__vertices)
        facearriba = faceSCabajo + np.array([zum, zum, zum])

        self._Gear__vertices = np.append(self._Gear__vertices, vertarriba, axis=0)
        self._Gear__faces = np.append(self._Gear__faces, facearriba, axis=0)
        # Aquí ya tengo un diente completo

        for i in range(0, len(self._Gear__faces)):
            we = self._Gear__faces[i]
            v1 = self._Gear__vertices[we[0]]
            v2 = self._Gear__vertices[we[1]]
            v3 = self._Gear__vertices[we[2]]
            normal = normal_to_surface(v1, v2, v3)
            self._Gear__normals = np.append(self._Gear__normals, [normal], axis=0)

        self._Gear__normals = np.delete(self._Gear__normals, 0, 0)
        norfinito = self._Gear__normals
        # Conversor a binerio es extremadamente importante. Crea una lista con los exponentes para 2**x de manera que controla la cantidad de dientes
        lista = conversoraBinario(self.z)

        self._Gear__normals = normales1(lista, self._Gear__az, self._Gear__normals)

        # Inserto el resto de los dientes
        vertieterno = self._Gear__vertices
        faceterno = self._Gear__faces
        facejuan3 = self._Gear__faces
        vertijuan3 = self._Gear__vertices
        lista = conversoraBinario(self.z)

        # Creo el inserto del primer grupo de dientes (recuerda que insertodientes1 afectamente la matriz final mientras que insertodientes2 no lo hace hasta el final
        dientes = insertodientes1(lista, self._Gear__az, self._Gear__vertices, self._Gear__faces, vertijuan3,
                                  facejuan3, vertieterno,
                                  faceterno)
        self._Gear__vertices = dientes[0]
        self._Gear__faces = dientes[1]
        self._Gear__normals = normales2(lista, self._Gear__az, self._Gear__normals, norfinito, norfinito)

        # Creo los insertos para el resto de los dientes
        dientes = insertodientes2(lista, self._Gear__az, self._Gear__vertices, self._Gear__faces, vertijuan3,
                                  facejuan3, vertieterno,
                                  faceterno)
        self._Gear__vertices = dientes[0]
        self._Gear__faces = dientes[1]

        MirroredVertices = MirrorX(self._Gear__vertices, mt.pi)
        MirroredNormals = MirrorX(self._Gear__normals, mt.pi)
        MirroredFaces = self._Gear__faces + np.array([len(self._Gear__vertices), len(self._Gear__vertices), len(self._Gear__vertices)])

        # Como hice un mirror las normales se invierten, por lo que es necesario cambiar una de las columnas de los extremos
        # por la de en medio (infalible)
        MirroredFaces[:, [1, 0]] = MirroredFaces[:, [0, 1]]

        self._Gear__faces = np.append(self._Gear__faces, MirroredFaces, axis=0)
        self._Gear__vertices = np.append(self._Gear__vertices, MirroredVertices, axis=0)
        self._Gear__normals = np.append(self._Gear__normals, MirroredNormals, axis=0)

        # Subo el engrane para que la cara inferior toque el suelo y no haya conflixtos con el codigo de agujeros.
        self._Gear__vertices  = self._Gear__vertices + np.array([0, 0, self.anchoeng])


        listaaguj = [self._Gear__vertices, self._Gear__faces, self._Gear__normals]
        # IMPORTANTE !! ESTO JALARIA COMO FUNCION O METODO porque DH Y DR comparten algoritmos DO solo HL es diferente
        # NOTA: HAZ QUE MATCHEEN CON LO QUE ARROJA TU UI -> HEX, SQR, CIR, KW, None
        if self._Gear__HoleData.Exists and self._Gear__HoleData.type == "HEX":
            listaaguj = DOhexagono(self._Gear__vertices, self._Gear__faces, self._Gear__normals, self.__anghrot, self.z, self.radiopoligono, 0, 2*self.anchoeng, self._Gear__HoleData.circumradius)
            pass
        elif self._Gear__HoleData.Exists and self._Gear__HoleData.type == "SQR":
            listaaguj = DOcuadrado(self._Gear__vertices, self._Gear__faces, self._Gear__normals, self.__anghrot, self.z, self.radiopoligono, 0, 2*self.anchoeng, self._Gear__HoleData.circumradius)
            pass
        elif self._Gear__HoleData.Exists and self._Gear__HoleData.type == "CIR":
            listaaguj = DOcirculo(self._Gear__vertices, self._Gear__faces, self._Gear__normals, -self.__anghrot, self.z, self.radiopoligono, 0, 2*self.anchoeng, self._Gear__HoleData.circumradius)
            pass
        elif self._Gear__HoleData.Exists and self._Gear__HoleData.type == "KW":
            #print("ES Keyway")
            listaaguj = DOkeyway(self._Gear__vertices, self._Gear__faces, self._Gear__normals, self.__anghrot, self.z, self.radiopoligono, 2*self.anchoeng, self._Gear__HoleData.boreDiameter, self._Gear__HoleData._alturaC, self._Gear__HoleData.keywayWidth)
            pass
        else:
            # el else solo aplica cuando es el agujero default o no lleva.
            pass

        self._Gear__vertices = listaaguj[0]
        self._Gear__faces = listaaguj[1]
        self._Gear__normals = listaaguj[2]


        # endregion

#endregion


