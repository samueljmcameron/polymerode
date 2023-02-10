import numpy as np

def convertVTKtoNParrays(fname):

    with open(fname) as vtkfile:
        while True:
            line = vtkfile.readline()
            if line.rstrip() == "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">":
                line = vtkfile.readline().rstrip()
                break
            elif (line == ""):
                print("Error: could not find data.")
                break



    words = line.split()

    Nbeads = len(words)//3

    xs = np.empty([Nbeads],float)
    ys = np.empty([Nbeads],float)
    zs = np.empty([Nbeads],float)

    bead = 0
    wi = 0

    com_calc = np.array([0,0,0],float)
    while (bead < Nbeads):


        xs[bead] = words[wi]
        ys[bead] = words[wi+1]
        zs[bead] = words[wi+2]

        com_calc[0] += xs[bead]
        com_calc[1] += ys[bead]
        com_calc[2] += zs[bead]

        bead += 1
        wi += 3



    return [xs,ys,zs]
