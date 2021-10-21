import numpy as np
import matplotlib.pyplot as plt
def avg_maxima_vs_k():
    filename = "/home/larissa/local_maxima_avg_vs_k.txt"
    M = np.loadtxt(filename)
    x1 = np.zeros(10)
    y1 = np.zeros(10)
    x2 = np.zeros(12)
    y2 = np.zeros(12)
    x3 = np.zeros(15)
    y3 = np.zeros(15)
    for i in range(1,len(M)):
        if(M[i][0]== 10):
          k = int(M[i][1])
          y1[k] = int(M[i][2])
          x1[k] = int(M[i][1])
        if(M[i][0] == 12):
          k = int(M[i][1])
          y2[k] = int(M[i][2])
          x2[k] = int(M[i][1])
        if(M[i][0]== 15):
          k = int(M[i][1])
          y3[k] = int(M[i][2])
          x3[k] = int(M[i][1])

    plt.plot(x1, y1, marker = 'o', lw = 3., ls = "-", mew = 2, color = "darkblue", label = "N = 10")
    plt.plot(x2, y2, marker = 'o', lw = 3., ls = "--", mew = 2, color = "darkorange", label = "N = 12")
    plt.plot(x3, y3, marker = 'o', lw = 3., ls = ":", mew = 2, color = "green", label = "N = 15")
    plt.xlim(1,11)
    plt.ylim(2,2000)
    plt.yscale('log')
    #plt.xscale('log')
    plt.xticks(np.arange(1, 11, step=1))
    plt.title("NK landscape with N = 12")
    plt.xlabel("K")
    plt.ylabel("average number of maxima")
    plt.legend(loc = "lower right")
    plt.savefig("/home/larissa/maxima_vs_k_sample1000.png")
    plt.show()
    return
def complexityCatastrophy():
    N = 16
    filename = "/home/larissa/Documents/nk model/graph1.txt"
    M = np.loadtxt(filename)
    print(M)
    x1 = []
    y1 = []
    x2 = []
    y2 = []
    x3 = []
    y3 = []
    y4 = []
    for i in range (0, len(M)):
        if(M[i][0] > 2):
            x1.append(M[i][1])
            y3.append(M[i][3])
            y1.append(M[i][2])
        if(M[i][0]==2):
            x2.append(M[i][1])
            y2.append(M[i][2])
            y4.append(M[i][3])
        #if(M[i][0]==8):
        #    x3.append(M[i][1])
        #    y3.append(M[i][2])

    plt.plot(x1, y1, marker = 'o', lw = 2., ls = ":", mew = 2, color = "darkblue", label = r'local maxima: $ K = N-1$')
    plt.plot(x1, y3, marker = 'o', lw = 2., ls = ":", mew = 2, color = "blue", label = r'global maximum: $  K = N-1$')

    plt.plot(x2, y2, marker = 'o', lw = 2., ls = "--", mew = 2, color = "darkorange", label = r'local maximum: $ K = 2$')
    plt.plot(x2, y4, marker = 'o', lw = 2., ls = "--", mew = 2, color = "r", label = r'global maximum: $ K = 2$')
    #plt.yscale('log')
    #plt.xscale('log')
    #plt.title("NK landscape with N = 12")
    #plt.xlim(1, 15)
    plt.xlabel(r'$N$', fontsize = 12)
    plt.ylabel(r'$\langle{\Phi}_{Max} \rangle$', fontsize = 14)
    plt.legend()
    plt.savefig("/home/larissa/Documents/nk model/complexity.png")
    plt.show()
    return

def f(x):
    y = []
    for i in range(0, len(x)):
        y.append(1/2)
    return y
def HammingDistanceMaxima():
    N = 15
    filename = "/home/larissa/maxima_hamming_distance_avg_vs_k.txt"
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

#avg_maxima_vs_k()
complexityCatastrophy()
#HammingDistanceMaxima()
