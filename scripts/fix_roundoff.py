import numpy as np
from scipy.optimize import root,show_options


def constraints(epsilons,uxs,uys,uzs):

    N = uxs.size

    b2s = uxs**2+uys**2+uzs**2

    return [np.sqrt(b2s)*np.sqrt(1+epsilons)-4.0,np.diag(np.sqrt(b2s)/np.sqrt(1+epsilons))]


def fix_roundoff(xs,ys,zs):

    uxs = xs[1:]-xs[:-1]
    uys = ys[1:]-ys[:-1]
    uzs = zs[1:]-zs[:-1]

    epsilons = uxs*0

    sol = root(constraints,epsilons,args=(uxs,uys,uzs),jac=True)

    
    uxs *= np.sqrt(1+sol.x)
    uys *= np.sqrt(1+sol.x)
    uzs *= np.sqrt(1+sol.x)

    rxs = [xs[0]]
    rys = [ys[0]]
    rzs = [zs[0]]

    for i in range(1,len(uxs)+1):

        rxs.append(rxs[-1]+uxs[i-1])
        rys.append(rys[-1]+uys[i-1])
        rzs.append(rzs[-1]+uzs[i-1])


    return np.array(rxs), np.array(rys), np.array(rzs)



if __name__ == "__main__":

    import sys
    path = sys.argv[1] + '/Documents/polymerode/scripts/'
    sys.path.append(path)
    from convertVTKtoNParrays import convertVTKtoNParrays
    import matplotlib.pyplot as plt

    proc = 1
    inputfile = f'vtkfiles/run1atom_p{proc}_100000.vtp'

    xs,ys,zs = convertVTKtoNParrays(inputfile,dname='ux')

    rxs,rys,rzs = fix_roundoff(xs,ys,zs)

    uxs = rxs[1:]-rxs[:-1]
    uys = rys[1:]-rys[:-1]
    uzs = rzs[1:]-rzs[:-1]


    
    fig = plt.figure()
    
    ax = fig.add_subplot(projection='3d')

    plt.plot(xs,ys,zs,'o')
    plt.plot(rxs,rys,rzs,'k-')
    plt.show()
