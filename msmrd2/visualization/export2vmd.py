# Sample functions to generate tcl files and visualize in vmd

def generateTCL_multiRotParticleMS(numparticles = 1, outfname = "testRotMulti", tclfname = "../../data/vmd/out2vmd.tcl"):
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


def generateTCL_gayBerne(numparticles = 10, outfname = "gayBerne", tclfname = "../../data/vmd/gayBerne2vmd.tcl"):
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


def generateTCL_patchyParticles(numparticles = 10, outfname = "patchyParticles", tclfname = "../../data/vmd/patchyParticles2vmd.tcl"):
    file = open(tclfname, 'w')
    file.write('set name ' + outfname + '\n')
    file.write('mol load xyz ./$name.xyz  \n \n')
    file.write('mol delrep 0 top \n')
    file.write('mol default style VDW \n')
    file.write('display resetview \n \n')

    # Define particle types
    # For main particle
    file.write('mol representation VDW 0.35000 0.5 \n')
    file.write('mol selection name type_0 \n')
    file.write('mol color ColorID 23 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')
    # For patches
    file.write('mol representation VDW 0.20000 0.5 \n')
    file.write('mol selection name type_1 \n')
    file.write('mol color ColorID 3 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')

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

def generateTCL_patchyParticlesMultiColor(numparticles = 10, outfname = "patchyParticles",
                                          tclfname = "../../data/vmd/patchyParticles2vmd.tcl"):
    file = open(tclfname, 'w')
    file.write('set name ' + outfname + '\n')
    file.write('mol load xyz ./$name.xyz  \n \n')
    file.write('mol delrep 0 top \n')
    file.write('mol default style VDW \n')
    file.write('display resetview \n \n')

    # Define particle types
    # For main particle 1
    file.write('mol representation VDW 0.35000 0.5 \n')
    file.write('mol selection name type_00 \n')
    file.write('mol color ColorID 23 \n')
    file.write('mol material Opaque \n')
    file.write('mol addrep top \n \n')
    # For main particle 2
    file.write('mol representation VDW 0.35000 0.5 \n')
    file.write('mol selection name type_01 \n')
    file.write('mol color ColorID 22 \n')
    file.write('mol material Opaque \n')
    file.write('mol addrep top \n \n')
    # For patches type1
    file.write('mol representation VDW 0.20000 0.5 \n')
    file.write('mol selection name type_1 \n')
    file.write('mol color ColorID 31 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')
    # For patches type2
    file.write('mol representation VDW 0.20000 0.5 \n')
    file.write('mol selection name type_2 \n')
    file.write('mol color ColorID 32 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')

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

def generateTCL_patchyProteins(numparticles = 2, outfname = "patchyProteins", tclfname = "../../data/vmd/patchyProteins2vmd.tcl"):
    file = open(tclfname, 'w')
    file.write('set name ' + outfname + '\n')
    file.write('mol load xyz ./$name.xyz  \n \n')
    file.write('mol delrep 0 top \n')
    file.write('mol default style VDW \n')
    file.write('display resetview \n \n')

    # Define particle types
    # For main particle
    file.write('mol representation VDW 0.35000 0.5 \n')
    file.write('mol selection name type_0 \n')
    file.write('mol color ColorID 23 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')
    # For normal patches
    file.write('mol representation VDW 0.20000 0.5 \n')
    file.write('mol selection name type_1 \n')
    file.write('mol color ColorID 3 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')
    # For binding patch
    file.write('mol representation VDW 0.150000 0.5 \n')
    file.write('mol selection name type_2 \n')
    file.write('mol color ColorID 30 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')

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


def generateTCL_patchyProteinsV2(numparticles = 2, outfname = "patchyProteins", tclfname = "../../data/vmd/patchyProteins2vmd.tcl"):
    file = open(tclfname, 'w')
    file.write('set name ' + outfname + '\n')
    file.write('mol load xyz ./$name.xyz  \n \n')
    file.write('mol delrep 0 top \n')
    file.write('mol default style VDW \n')
    file.write('display resetview \n \n')

    # Define particle types
    # For main particles
    file.write('mol representation VDW 0.35000 0.5 \n')
    file.write('mol selection name type_0 \n')
    file.write('mol color ColorID 23 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')
    file.write('mol representation VDW 0.35000 0.5 \n')

    file.write('mol selection name type_1 \n')
    file.write('mol color ColorID 16 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')

    # For normal patches
    file.write('mol representation VDW 0.20000 0.5 \n')
    file.write('mol selection name type_2 \n')
    file.write('mol color ColorID 3 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')
    # For binding patch
    file.write('mol representation VDW 0.20000 0.5 \n')
    file.write('mol selection name type_3 \n')
    file.write('mol color ColorID 30 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')

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


def generateTCL_patchyProteinsMS(numparticles = 2, outfname = "patchyProteinsMS", tclfname = "../../data/vmd/patchyProteinsMS2vmd.tcl"):
    file = open(tclfname, 'w')
    file.write('set name ' + outfname + '\n')
    file.write('mol load xyz ./$name.xyz  \n \n')
    file.write('mol delrep 0 top \n')
    file.write('mol default style VDW \n')
    file.write('display resetview \n \n')

    # Define particle types
    # For main particle
    file.write('mol representation VDW 0.35000 0.5 \n')
    file.write('mol selection name type_0 \n')
    file.write('mol color ColorID 23 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')
    # For normal patches
    file.write('mol representation VDW 0.20000 0.5 \n')
    file.write('mol selection name type_1 \n')
    file.write('mol color ColorID 3 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')
    # For binding patch
    file.write('mol representation VDW 0.150000 0.5 \n')
    file.write('mol selection name type_2 \n')
    file.write('mol color ColorID 30 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')
    # For alternative particle
    file.write('mol representation VDW 0.35000 0.5 \n')
    file.write('mol selection name type_3 \n')
    file.write('mol color ColorID 15 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')
    # For alternative patch
    file.write('mol representation VDW 0.15000 0.5 \n')
    file.write('mol selection name type_4 \n')
    file.write('mol color ColorID 6 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')

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


def generateTCL_patchyProteinsMSV2(numparticles = 2, outfname = "patchyProteinsMS", tclfname = "../../data/vmd/patchyProteinsMS2vmd.tcl"):
    file = open(tclfname, 'w')
    file.write('set name ' + outfname + '\n')
    file.write('mol load xyz ./$name.xyz  \n \n')
    file.write('mol delrep 0 top \n')
    file.write('mol default style VDW \n')
    file.write('display resetview \n \n')

    # Define particle types
    # For main particle
    file.write('mol representation VDW 0.35000 0.5 \n')
    file.write('mol selection name type_0 \n')
    file.write('mol color ColorID 23 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')
    # For normal patches
    file.write('mol representation VDW 0.20000 0.5 \n')
    file.write('mol selection name type_1 \n')
    file.write('mol color ColorID 3 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')
    # For binding patch
    file.write('mol representation VDW 0.150000 0.5 \n')
    file.write('mol selection name type_2 \n')
    file.write('mol color ColorID 30 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')

    # For licorice particles
    file.write('mol representation Licorice 0.2000 0.5 \n')
    file.write('mol selection name type_3 \n')
    file.write('mol color ColorID 23 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')

    # For licorice particles alt
    file.write('mol representation VDW 0.14000 0.5 \n')
    file.write('mol selection name type_4 \n')
    file.write('mol color ColorID 30 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')

    # Define bonds within particle
    file.write('set sel [atomselect 0 "index 7 8 9 10 11 12 13"] \n')
    file.write('$sel setbonds {{8 9 10 11 12 13} {7} {7} {7} {7} {7} {7}} \n')
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


def generateTCL_pentamerTest(numparticles = 10, outfname = "pentamerTest", tclfname = "../../data/vmd/pentamerTest2vmd.tcl"):
    file = open(tclfname, 'w')
    file.write('set name ' + outfname + '\n')
    file.write('mol load xyz ./$name.xyz  \n \n')
    file.write('mol delrep 0 top \n')
    file.write('mol default style VDW \n')
    file.write('display resetview \n \n')

    # Define particle types
    # For main particle
    file.write('mol representation VDW 0.35000 0.5 \n')
    file.write('mol selection name type_0 \n')
    file.write('mol color ColorID 23 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')
    # For patches
    file.write('mol representation VDW 0.20000 0.5 \n')
    file.write('mol selection name type_1 \n')
    file.write('mol color ColorID 3 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')

    # For center
    file.write('mol representation VDW 0.10000 0.5 \n')
    file.write('mol selection name type_2 \n')
    file.write('mol color ColorID 29 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')

    # For licorice particles
    file.write('mol representation lines 2 \n')
    file.write('mol selection name type_3 \n')
    file.write('mol color ColorID 16 \n')
    file.write('mol material AOShiny \n')
    file.write('mol addrep top \n \n')

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