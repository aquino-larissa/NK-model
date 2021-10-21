import numpy as np
import matplotlib.pyplot as plt

def HammingDistanceMaxima():
    N = 15 
    filename = "/home/larissa/Documentos/modelo NK/Hamming Distance/graph1.txt"
    M = np.loadtxt(filename)
    print(M)
    x1 = []
    y1 = []
    x2 = []
    y2 = []
    x3 = []
    y3 = []
    for i in range (0, len(M)):
        if (M[i][0] == 10):
            x1.append(M[i][1])
            y1.append( M[i][2]/M[i][0])
        if (M[i][0] == 12):
            x2.append(M[i][1])
            y2.append( M[i][2]/M[i][0])
        if (M[i][0] == 15):
            x3.append(M[i][1])
            y3.append( M[i][2]/M[i][0])
    plt.plot(x1, y1, marker = 'o', lw = 3., ls = "-", mew = 2, color = "darkblue", label = r'$N = 10$')
    plt.plot(x2, y2, marker = 'o', lw = 3., ls = "-.", mew = 2, color = "darkorange", label = r'$N = 12$')
    plt.plot(x3, y3, marker = 'o', lw = 3., ls = ":", mew = 2, color = "green", label = r'$N = 15$')
    plt.plot(x3, f(x3), ls = "--", lw = 3., color = "black")
    plt.xlim(1,14)
    plt.xlabel(r'$K$')
    plt.ylabel(r'$\langle D_{H} \rangle/N $', fontsize = 14, fontweight ="bold")
    plt.legend(loc = "lower right")
    plt.savefig("/home/larissa/avg_hamming_distance_maxima.png")
    plt.show()
    return

def fitnessHammingDistance():
    N = 15 
    filename = "/home/larissa/Documentos/modelo NK/Hamming Distance/graph2-N15.txt"
    M = np.loadtxt(filename)
    print(M)
    x1 = []
    y1 = []
    x2 = []
    y2 = []
    x3 = []
    y3 = []
    x4 = []
    y4 = []
    for i in range (0, len(M)):
        a = 3
        if (M[i][1] == 0):
            x1.append(M[i][2])
            y1.append( M[i][a])
        if (M[i][1] == 2):
            x2.append(M[i][2])
            y2.append( M[i][a])
        if (M[i][1] == 4):
            x3.append(M[i][2])
            y3.append( M[i][a])
        if (M[i][1] == 8):
            x4.append(M[i][2])
            y4.append( M[i][a])
    plt.plot(x1, y1, marker = 'o', lw = 3., ls = "-", mew = 2, color = "darkblue", label = r'$ K = 0 $')
    plt.plot(x2, y2, marker = 'o', lw = 3., ls = "-.", mew = 2, color = "darkorange", label = r'$ K = 2 $')
    plt.plot(x3, y3, marker = 'o', lw = 3., ls = ":", mew = 2, color = "green", label = r'$ K = 4 $')
    plt.plot(x4, y4, marker = 'o', lw = 3., ls = "-", mew = 2, color = "black", label = r'$ K = 8$')
    #plt.ylim(1,1000)
    plt.xlabel(r'$D_H$')
    plt.ylabel(r'$\bar{f}(x)/f(max) $', fontsize = 14, fontweight ="bold")
    plt.legend(loc = "upper right")
    plt.savefig("/home/larissa/Documentos/modelo NK/Hamming Distance/graph2-N15.png")
    plt.show()
    return
    
fitnessHammingDistance()
#HammingDistanceMaxima()
