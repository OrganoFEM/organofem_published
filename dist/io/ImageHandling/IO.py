from os.path import exists, splitext, split as psplit, expanduser as expusr
import os, fnmatch

# from scipy.misc import imsave as _imsave
# this is not used (GM)
from struct import pack,unpack,calcsize
from pickle import dumps,loads
import numpy as np

from .spatial_image import SpatialImage
from .inrimage import read_inrimage, write_inrimage
from .metaimage import read_metaimage, write_metaimage
from .tif import read_tif, write_tif # OLD VERSION LIBTIFF depreciated
from .h5 import read_h5

def imread (filename, dimension=3) :
    """Reads an image file completely into memory.

    It uses the file extension to determine how to read the file. It first tries
    some specific readers for volume images (Inrimages, TIFFs, LSMs, NPY) or falls
    back on PIL readers for common formats if installed.

    In all cases the returned image is 3D (2D is upgraded to single slice 3D).
    If it has colour or is a vector field it is even 4D.

    :Parameters:
     - `filename` (str)

    :Returns Type:
        |SpatialImage|
    """
    program='imread'
    filename = expusr(filename)

    if not os.path.isfile(filename) and os.path.isfile(filename+".gz"):
        filename=filename+".gz"
        print ("%s: Warning: path to read image has been changed to %s."%(program,filename))

    if not os.path.isfile(filename) and os.path.isfile(filename+".zip"):
        filename=filename+".zip"
        print ("%s: Warning: path to read image has been changed to %s."%(program,filename))

    if not exists(filename) :
        raise IOError("The requested file do not exist: %s" % filename)

    root, ext = splitext(filename)
    ext = ext.lower()
    if ext == ".gz":
        root, ext = splitext(root)
        ext = ext.lower()
    if ext == ".inr":
        return read_inrimage(filename)
    if ext == ".mha":
        return read_metaimage(filename)
    elif ext in [".tif", ".tiff"]:
        return read_tif(filename)
    elif ext in [".npz", ".npy"]:
        return load(filename)
    elif ext in [".h5", ".hdf5"]:
        return read_h5(filename)


def imsave(filename, img):
    """Save a |SpatialImage| to filename.

    .. note: `img` **must** be a |SpatialImage|.

    The filewriter is choosen according to the file extension. However all file extensions
    will not match the data held by img, in dimensionnality or encoding, and might raise `IOError`s.

    For real volume data, Inrimage and NPY are currently supported.
    For |SpatialImage|s that are actually 2D, PNG, BMP, JPG among others are supported if PIL is installed.

    :Parameters:
     - `filename` (str)
     - `img` (|SpatialImage|)
    """

    assert isinstance(img, SpatialImage)
    # -- images are always at least 3D! If the size of dimension 3 (indexed 2) is 1, then it is actually
    # a 2D image. If it is 4D it has vectorial or RGB[A] data. --
    filename = expusr(filename)
    head, tail = psplit(filename)
    head = head or "."
    if not exists(head):
        raise IOError("The directory do not exist: %s" % head)

    root, ext = splitext(filename)

    is2D = img.shape[2] == 1
    ext = ext.lower()
    if ext == ".gz":
        root, ext = splitext(root)
        ext = ext.lower()
    if ext == ".inr":
        write_inrimage(filename, img)
    if ext == ".mha":
        write_metaimage(filename, img)
    elif ext in [".npz", ".npy"]:
        save(filename, img)
    elif ext in [".tiff", ".tif"]:
        write_tif(filename, img)
