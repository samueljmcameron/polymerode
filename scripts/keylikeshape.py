import numpy as np



def keyLikeShape(bondlength,N,rotate=True):

    if (N % 2 == 0):
        beadsincircle=8
    else:
        beadsincircle=7


    phi = 2*np.pi/(beadsincircle+1)
    R = bondlength/(2*np.sin(phi/2.))

    dx = bondlength*np.cos(phi/2.)


    xs = np.zeros([N],float)
    ys = np.zeros([N],float)
    zs = np.zeros([N],float)


    offset = (N-beadsincircle)//2



    for i in range(offset):

        xs[i] = -dx
        xs[N-1-i] = dx
        ys[i] = i*bondlength
        ys[N-1-i] = i*bondlength


    dy = np.sqrt(bondlength**2-dx**2)
    yc = ys[offset-1]+bondlength - dy + R

    for i in range(offset,offset+beadsincircle):

        theta = (i+1-offset)*phi


        xs[i] = -R*np.sin(theta)
        ys[i] = -R*np.cos(theta)+yc


    if rotate:
        theta = np.arccos(1/np.sqrt(3))
        phi = np.arccos(1./np.sqrt(3)/np.sin(theta))

        # rotate the structure so that it points along (1,1,1)
        Q1 = np.array([[np.cos(phi),np.sin(phi),0],
                      [-np.sin(phi),np.cos(phi),0],
                      [0,0,1]])


        Q2 = np.array([[1,0,0],
                            [0,np.cos(theta),-np.sin(theta)],
                            [0,np.sin(theta),np.cos(theta)]])



        rot = Q2 @ Q1



        for i in range(N):
            # rotate vectors counter clockwise by the angle so
            # that middle of curve is at phi = pi/4.
            vec = np.array([xs[i],ys[i],zs[i]])

            newvec = rot @ vec


            xs[i] = newvec[0]
            ys[i] = newvec[1]
            zs[i] = newvec[2]



    return xs,ys,zs


if __name__=="__main__":

    import matplotlib.pyplot as plt
    
    def checkbonds(xs,ys,zs,bondlength):
        bxs = xs[1:]-xs[:-1]
        bys = ys[1:]-ys[:-1]
        bzs = zs[1:]-zs[:-1]

        return np.sqrt(bxs**2+bys**2+bzs**2)-bondlength


        
    bondlength=4
    N=89

    radius = 40.
    xs,ys,zs = keyLikeShape(bondlength,N)



    xs += radius/np.sqrt(3)
    ys += radius/np.sqrt(3)
    zs += radius/np.sqrt(3)



    dbs = checkbonds(xs,ys,zs,bondlength)
    index  = np.argmax(np.abs(dbs))
    print(index,dbs[index])

    fig,ax = plt.subplots(subplot_kw={'projection':'3d'})
    ax.plot(xs,ys,zs,'o-')



    # draw sphere
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = radius*np.cos(u)*np.sin(v)
    y = radius*np.sin(u)*np.sin(v)
    z = radius*np.cos(v)
    ax.plot_wireframe(x, y, z, color="r",alpha=0.5)

    L = 256.
    ax.set_box_aspect([1.0,1.0,1.0])
    ax.set_ylabel('y')
    ax.set_xlabel('x')
    ax.set_zlabel('z')
    ax.set_xlim(-L/2,L/2)
    ax.set_ylim(-L/2,L/2)
    ax.set_zlim(-L/2,L/2)

    plt.show()

    xs,ys,zs = keyLikeShape(bondlength,N,rotate=False)
    fig,ax = plt.subplots()
    ax.plot(xs,ys,'o-')
    ax.set_aspect('equal')
    plt.show()


