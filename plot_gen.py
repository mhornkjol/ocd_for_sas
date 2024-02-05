
import pyvista as pv 
import glob


data_files = ["initial000000","final000000"]

# missing DTI for 78, 127, 172, 249 

pats = [105, 175, 178,  190,  199,  215,  227,  230,  236,  105, 176,  183,  191,  205,  218,  228, 235,  240]


for pat in pats: 
    files = []
    #data files 
    pat = f'{pat:03}' 
    files.append("%s/results/initial000000.vtu" % pat) 
    files.append("%s/results/final000000.vtu" % pat) 

    #result files: first, last  
    ufiles = glob.glob("%s/results/uDTI*vtu"%pat)
    ufiles.sort()
    files.append(ufiles[0])
    files.append(ufiles[-1])
    
    print (files) 
    for f in files:  
        
        u = pv.read(f)
        u_slice = u.slice(normal=(0,0,1))

        plotter = pv.Plotter(off_screen=True)
        plotter.camera_position = [(0,0,300), (0, 0, 0), (0, 0, 0)]
#        plotter.add_mesh_slice_orthogonal(u)
        plotter.add_mesh(u_slice)
        plotter.update_scalar_bar_range([0,300])
        plotter.screenshot("pics/%s.png" % (f.replace("/", "").replace(".", "")))
