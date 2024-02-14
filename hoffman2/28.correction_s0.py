### apply flatfield correction to s0 image

import numpy as np
from scipy.ndimage import zoom
import zarr

from matplotlib import pyplot as plt
import seaborn as sns
import tifffile

import logging
logging.basicConfig(level=logging.INFO)

def transform(images, f, d=None):
    if d is not None:
        images_transformed = (np.clip(images-d[np.newaxis],0,None))/f[np.newaxis]
    else:
        images_transformed = images/f[np.newaxis]
        
    return images_transformed

def transform_pipe(images_raw, dx, ff, upscale=8):
    """
    """
    ff_s = zoom(ff, upscale).astype(np.float16)
    dx_s = zoom(dx, upscale).astype(np.float16)
    images_s_transformed = transform(images_raw.astype(np.float16), ff_s, d=dx_s).astype(np.uint16)
    return images_s_transformed

f_dx = '/u/home/f/f7xiesnm/project-zipursky/easifish/results/sparse06_r1c0_flatfield/darkfield_rb.tiff' 
f_ff = '/u/home/f/f7xiesnm/project-zipursky/easifish/results/sparse06_r1c0_flatfield/flatfield_rb.tiff' 

path = "/u/scratch/f/f7xiesnm/sparse06/dataset.n5"
outpath = "/u/home/f/f7xiesnm/project-zipursky/data/hold/sparse06/r1_test_flatfield_v8.n5"

scale   = 's0' #'s0'
upscale = 16    # 16

zarr_data = zarr.open(store=zarr.N5Store(path), mode='r')
# slots_h1 = list(zarr_data.keys())
# slots_h1 = list(zarr_data.keys())[23:]
# slots = [f'/{h1}/timepoint0/{scale}' for h1 in slots_h1]
i_s = [4, 5, 6, 7, 8, 9]
slots = [f'/setup{i}/timepoint0/{scale}' for i in i_s]
logging.info(f"Number of slots: {len(slots)}")
logging.info(f"{slots}")

dx = tifffile.imread(f_dx).astype(np.float32)
ff = tifffile.imread(f_ff).astype(np.float32)

n5_root = zarr.open_group(store=zarr.N5Store(outpath), mode='r+')
for slot in slots: # this can be parallelized
    logging.info(slot)

    images_s0_handle  = zarr_data[slot]
    attributes = images_s0_handle.attrs.asdict()
    attributes['pixelResolution'] = [0.23, 0.23, 0.42] # useful later

    # get
    logging.info('retrieving the image...')
    images_s0_raw = images_s0_handle[...]

    # transform
    logging.info('transforming the image...')
    a = transform_pipe(images_s0_raw, dx, ff, upscale=upscale)

    # save
    logging.info('saving the image...')
    dataset = n5_root.require_dataset(
        slot,
        data=a,
        shape=a.shape,
        chunks=images_s0_handle.chunks, # (64, 128, 128),
        dtype=images_s0_handle.dtype,
        compressor=images_s0_handle.compressor,  # GZip(level=1),
        )
    dataset.attrs.update(**attributes)