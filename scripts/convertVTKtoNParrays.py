import numpy as np

def convertVTKtoNParrays(fname,dname=None,dtype='Float64',ncmps=3):

    with open(fname) as vtkfile:
        while True:
            line = vtkfile.readline()
            if dname == None:
                matchline = "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
            else:
                matchline = f"<DataArray Name=\"{dname}\" type=\"{dtype}\" NumberOfComponents=\"{ncmps}\" format=\"ascii\">"
            if line.rstrip() == matchline:
                line = vtkfile.readline().rstrip()
                break
            elif (line == ""):
                print("Error: could not find data.")
                break



    words = line.split()

    Nbeads = len(words)//ncmps

    outs = []


        

    for i in range(ncmps):
        if dtype == 'Float64':
            outs.append(np.empty([Nbeads],float))
        else:
            outs.append(np.empty([Nbeads],int))

    bead = 0
    wi = 0


    while (bead < Nbeads):


        for i in range(ncmps):
            outs[i][bead] = words[wi]
            wi += 1


        bead += 1

    return outs
