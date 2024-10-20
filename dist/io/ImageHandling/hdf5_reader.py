from h5py import File
import numpy as np
#from openalea.spatial_image import SpatialImage

def imread_h5(path):
    """Read an hdf5 file

    :Parameters:
     - `filename` (str) - name of the file to read
    """
    data = File(path)['Data']
    im_out = np.zeros(data.shape, dtype=data.dtype)
    data.read_direct(im_out)
    return SpatialImage(im_out.transpose(2, 1, 0))
    
    
