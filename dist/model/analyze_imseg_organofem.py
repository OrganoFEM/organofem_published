

#general functions
import os
import os.path
from os.path import exists as file_exists
import sys

sys.path.append('../../../') #to load organofem functions
import organofem
from organofem.dist.io.ImageHandling.IO import SpatialImage, imread, imsave

#image analysis specifac functions
from scipy import ndimage as nd
import numpy as np

### variables ####
##################
image_to_analyse="seg_organoid.tif" # in .mha or .tif format
vsize=0.32
# define the z step used during acquisition to transform voxel size in z
z_interpolation =1 
# used to reduce image size to test the code and decreasing the calculation time
subsample=1 # for no subsamble value = 1


#Surface of Contacts
def __slices_dilation(slices, maximum=[np.inf, np.inf, np.inf]):
    return tuple([slice(max(0, s.start-1), min(s.stop+1, maximum[i])) for i, s in enumerate(slices)])

def slices_dilation(slices, maximum=[np.inf, np.inf, np.inf], iterations=1):
    for i in range(iterations):
        slices=__slices_dilation(slices, maximum)
    return slices    

def process_contact_surface(seg_file_format):
    t=1
    seg=imread(seg_file_format)
    bboxes=nd.find_objects(seg)
    cells=np.unique(seg)
    contact_surf={}
    seg[seg==0]=1
    labels=cells[2:]
    barycenters_tmp=dict(zip(t*10**4+labels, nd.center_of_mass(np.ones_like(seg), seg, labels)))
    barycenters={k:(v[0]*2, v[1]*2, v[2]*2) for k, v in barycenters_tmp.items()}
    for c in cells:
        if c!=1 and c!=0:
            print("processing cell id: ",c)
            bb=slices_dilation(bboxes[c-1], iterations=2)
            tmp=seg[bb]
            bounds=tmp[(nd.binary_dilation(tmp==c) & (tmp!=c))]
            contact_surf[c]=dict(zip(np.unique(bounds), nd.sum(np.ones_like(bounds), bounds, np.unique(bounds))))
    surfaces={}
    for k, v in contact_surf.items():
        surfaces[t*10**4+k]={t*10**4+key:val for key, val in v.items()}
    return contact_surf, surfaces, barycenters

def process_multiple_contacts(seg_file_format):
    print("Process Surface of Contact")
    contact_surf, surfaces, barycenters, inertia_axis, eig_values = {}, {}, {}, {}, {}
    cs, s, b = process_contact_surface(seg_file_format)
    contact_surf.update(cs)
    surfaces.update(s)
    barycenters.update(b)
    return surfaces,barycenters,contact_surf


# define voxel size properly if subsampling is required
vsize=vsize/subsample

#make a save of a correctly writted tif 
try:
    image_to_analyse=image_to_analyse.replace(".mha", ".tif")
    imsave(image_to_analyse,imread(image_to_analyse))
except:
    pass


rescaled_image=image_to_analyse.replace(".tif", "_res.mha")


if file_exists(rescaled_image)==True:
    D_interpolated=imread(rescaled_image)
    os.chdir('../../test/tmp')

else:
    seg1=imread(image_to_analyse)
    os.chdir('../../test/tmp')

    Nx, Ny, Nz = seg1.shape[0],seg1.shape[1],seg1.shape[2]
    
    import scipy.interpolate as interpolate
    
    # Create a uniformly spaced grid
    Mx, My, Mz =int(float(seg1.shape[0]*subsample)),  int(float(seg1.shape[1]*subsample)),  int(float(seg1.shape[2]*subsample*z_interpolation))
    
    x, y, z = seg1.shape
    X, Y, Z = np.mgrid[:x,:y,:z]
    print("initial shape")
    print(Nx, Ny, Nz)
    X, Y, Z = np.broadcast_arrays(X, Y, Z)
    print("resampling factor used: "+str(subsample))
    print("new shape with resampling and isotropic resolution")
    print(Mx, My, Mz)
    # Create a uniformly spaced grid
    xi = np.linspace(X.ravel().min(), X.ravel().max(), Mx)
    yi = np.linspace(Y.ravel().min(), Y.ravel().max(), My)
    zi = np.linspace(Z.ravel().min(), Z.ravel().max(), Mz)
    X_uniform, Y_uniform, Z_uniform = (
        xi[:, None, None], yi[None, :, None], zi[None, None, :])
    
    print("begin of interpolation")
    
    D_interpolated = interpolate.griddata(
        (X.ravel(), Y.ravel(), Z.ravel()),
        seg1.ravel(),
        (X_uniform, Y_uniform, Z_uniform),
        method='nearest')
    
    print("end of interpolation")
    D_interpolated+=2 # to avoid problem with background at 1 and no measurments of stats in "process_multiple_contacts"
    imsave(rescaled_image, SpatialImage(D_interpolated, voxelsize=(vsize*subsample,vsize*subsample,vsize*subsample)).astype(np.uint16))


#calcul des volumes
unique, counts = np.unique(D_interpolated.ravel(), return_counts=True)
volumes=np.asarray((unique, counts*vsize*vsize*vsize)).T
#np.savetxt("volumes.csv", volumes, fmt="%i", delimiter=",") # possibility to save values in individual .csv files


import pandas as pd
Stats=pd.DataFrame(volumes[:,1],columns=['volume'])

id_ext=Stats['volume'].nlargest(2).keys()[0]
id_lum=Stats['volume'].nlargest(2).keys()[1]

print(Stats)
print(Stats['volume'].nlargest(2))

#calcul des surfaces et barycentres from ASTEC
surfaces,barycenters,contact_surf=process_multiple_contacts(rescaled_image)

id_ext=Stats['volume'].nlargest(2).keys()[0]
id_lum=Stats['volume'].nlargest(2).keys()[1]
id_ext=list(contact_surf.keys())[id_ext]
id_lum=list(contact_surf.keys())[id_lum]
print("matrigel id found: ",id_ext)
print("lumen id found: ",id_lum)

apical_surf={}
apical_lateral_surf={}
neibhoors_number={}

for x in contact_surf:
    surfL=0
    surfA=0
    nlat=0
    for i,j in contact_surf[x].items():
        if i== id_lum:
            surfA+=j
        elif i== id_ext:
            pass
        else:
            surfL+=j
            if j > 10:
                nlat+=1
    if surfA==0:
        apical_lateral_surf[x]=0
    elif surfL != 0:
        apical_lateral_surf[x]=surfA/surfL*2
    else:
        apical_lateral_surf[x]=0
    apical_surf[x]=surfA*vsize*vsize
    neibhoors_number[x]=nlat

basal_surf={}
basal_lateral_surf={}
for x in contact_surf:
    surfL=0
    surfB=0
    for i,j in contact_surf[x].items():
        if i== id_ext:
            surfB+=j
        elif i== id_lum:
            pass        
        else:
            surfL+=j
    if surfB==0:
        basal_lateral_surf[x]=0
    elif surfL != 0:
        basal_lateral_surf[x]=surfB/surfL*2
    else:
        basal_lateral_surf[x]=0

    basal_surf[x]=surfB*vsize*vsize

from sklearn.neighbors import KDTree
bary = np.array(list(barycenters.values()))
tree= KDTree(bary)
nearst_dist, nearest_ind = tree.query(bary, k=5)
nearst_dist = nearst_dist
near_list=[]

for x in nearst_dist:
    mean_near=np.mean(x)*vsize
    near_list.append(mean_near)
near_list.append(np.mean(near_list))
near_list.append(np.mean(near_list))

Stats["spacing"]=near_list

id_ext=Stats['volume'].nlargest(1).keys()[0]
#Stats = Stats.drop(labels=[id_ext], axis=0)
Stats["apical"]=apical_surf.values()
Stats["basal"]=basal_surf.values()
Stats["apical-lateral"]=apical_lateral_surf.values()
Stats["basal-lateral"]=basal_lateral_surf.values()
Stats["neibhoors_number"]=neibhoors_number.values()


id_lum=Stats['volume'].nlargest(1).keys()[0]
Stats = Stats.drop(labels=[id_lum], axis=0)

new_csv=image_to_analyse.replace(".tif","_stats.csv")
new_csv2=image_to_analyse.replace(".tif","_stats_post.csv")

Stats.to_csv(new_csv)

Stats = Stats[Stats.apical > 1] #do not take in accout cells that do not touch lumen (error of segmentation)
Stats = Stats[Stats.basal > 1] #do not take in accout cells that do not touch matrigel (error of segmentation)
Stats = Stats[Stats.volume > 500] #do not take tiny cells (error of segmentation)

Stats["basal_apical"]=Stats.basal/Stats.apical
Stats.to_csv(new_csv2)
print("Final stats saved in: "+str(new_csv2))
