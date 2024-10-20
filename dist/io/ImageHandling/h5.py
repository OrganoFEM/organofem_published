from h5py import File
#from spatial_image import SpatialImage
import numpy as np

def read_h5(path):
    """Read an hdf5 file

    :Parameters:
     - `filename` (str) - name of the file to read
    """
    data = File(path.replace('\\', ''), 'r')['Data']
    im_out = np.zeros(data.shape, dtype=data.dtype)
    data.read_direct(im_out)
    return SpatialImage(im_out.transpose(2, 1, 0))