# Sample functions to generate tcl files and visualize in vmd

def generateTCL_multiRotParticleMS(numparticles = 1, outfname = "testRotMulti", tclfname = "../data/vmd/out2vmd.tcl"):
    file = open(tclfname, 'w')
    file.write('set name ' + outfname + '\n')
    file.write('mol load xyz ./$name.xyz  \n \n')
    file.write('mol delrep 0 top \n')
    file.write('mol default style CPK \n')
    file.write('display resetview \n \n')

    # Define particle types (first one for bonds, second two for big and small spheres)
    file.write('mol representation CPK 0.20000 0.1 \n')
    file.write('mol selection name (type_0 and type_1) \n')
    file.write('mol color ColorID 16 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')
    
    file.write('mol representation VDW 0.20000 0.1 \n')
    file.write('mol selection name type_0 \n')
    file.write('mol color ColorID 1 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')
    
    file.write('mol representation VDW 0.10000 0.1 \n')
    file.write('mol selection name type_1 \n')
    file.write('mol color ColorID 0 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')
    
    # Define bonds within particle
    for i in range(numparticles):
        file.write('set sel [atomselect 0 "index ' + str(3*i) + ' ' + str(3*i+1) + ' ' + str(3*i+2) + '"]  \n')
        file.write('$sel setbonds {{' + str(3*i+1) + ' ' + str(3*i+2) + '} {' + str(3*i) + '} {' + str(3*i) + '}} \n')
    file.write(' \n')
    
    # Define display options
    file.write('axes location off \n')
    file.write('color Display Background white \n')
    file.write('display projection orthographic \n')
    file.write('display resize 800 800 \n')
    file.write('display nearclip set 0.0 \n')
    file.write('display depthcue off \n')
    file.write('#display cuedensity 0.20000  \n')
    file.write('#display cuemode Exp2  \n')
    file.write('display shadows on \n')
    file.write('display ambientocclusion on \n')
    file.write('display antialias on \n')
    file.write('display rendermode GLSL \n')
    file.close()

def generateTCL_dipole(frame, numparticles = 36, outfname = "odLangevinDipole", tclfname = "odLangevinDipole2vmd.tcl"):
    file = open(tclfname, 'w')
    file.write('set name ' + outfname + '\n')
    #file.write('mol new ./$name.xyz first 300 last 320 skip 2') # Load only a subset of frames.
    file.write('mol load xyz ./$name.xyz  \n \n')
    file.write('mol delrep 0 top \n')
    file.write('mol default style Licorice \n')
    file.write('display resetview \n \n')

    # Define two particle types (each on for one state)
    file.write('mol representation Licorice 0.10000 0.1 \n')
    file.write('mol selection name type_0 \n')
    file.write('mol color ColorID 1 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')

    file.write('mol representation Licorice 0.10000 0.1 \n')
    file.write('mol selection name type_1 \n')
    file.write('mol color ColorID 0 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')

    # Define bonds within particle
    for i in range(numparticles):
        file.write('set sel [atomselect 0 "index ' + str(2*i) + ' ' + str(2*i+1) + '"]  \n')
        file.write('$sel setbonds {{' + str(2*i+1) + '} {' + str(2*i) + '}} \n')
    file.write(' \n')

    # Define display options
    file.write('axes location off \n')
    file.write('color Display Background white \n')
    file.write('display projection orthographic \n')
    file.write('display resize 800 800 \n')
    file.write('display nearclip set 0.0 \n')
    file.write('display depthcue off \n')
    file.write('#display cuedensity 0.20000  \n')
    file.write('#display cuemode Exp2  \n')
    file.write('display shadows on \n')
    file.write('display ambientocclusion on \n')
    file.write('display antialias on \n')
    file.write('display rendermode GLSL \n')
    file.write('scale 1 0.5 \n')

    # Render frame and close vmd
    if (frame!=-1):
        file.write('render Tachyon filename.dat "/home/maojrs/Downloads/vmd-1.9.3/lib/tachyon/tachyon_LINUXAMD64" filename.dat '
               '-aasamples 12 -format TARGA -res 800 800 -o dipoleframe_' + str(frame).zfill(4) + '.tga \n')
        file.write('quit')
        file.close()


def generateTCL_gayBerne(numparticles = 10, outfname = "gayBerne", tclfname = "../data/vmd/gayBerne2vmd.tcl"):
    file = open(tclfname, 'w')
    file.write('set name ' + outfname + '\n')
    file.write('mol load xyz ./$name.xyz  \n \n')
    file.write('mol delrep 0 top \n')
    file.write('mol default style Licorice \n')
    file.write('display resetview \n \n')

    # Define particle types (first one for bonds, second two for big and small spheres)
    file.write('mol representation Licorice 0.50000 0.5 \n')
    file.write('mol selection name type_0 \n')
    file.write('mol color ColorID 1 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')

    # Define bonds within particle
    for i in range(numparticles):
        file.write('set sel [atomselect 0 "index ' + str(2*i) + ' ' + str(2*i+1) + '"]  \n')
        file.write('$sel setbonds {{' + str(2*i+1) + '} {' + str(2*i) + '}} \n')
    file.write(' \n')

    # Define display options
    file.write('axes location off \n')
    file.write('color Display Background white \n')
    file.write('display projection orthographic \n')
    file.write('display resize 800 800 \n')
    file.write('display nearclip set 0.0 \n')
    file.write('display depthcue off \n')
    file.write('#display cuedensity 0.20000  \n')
    file.write('#display cuemode Exp2  \n')
    file.write('display shadows on \n')
    file.write('display ambientocclusion on \n')
    file.write('display antialias on \n')
    file.write('display rendermode GLSL \n')
    file.close()


