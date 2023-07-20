import os
import sys
import numpy as np
import pandas as pd
from glob import glob
from os.path import abspath, dirname
from skimage.measure import regionprops
from scipy.spatial import cKDTree

def spot_counts(lb, spot_dir, s=[0.92,0.92,0.84], verbose=True):
    """
    Returns spot counts for each ROI. 
    
    lb: segmenntation mask
    spot_dir: Accepts 3 different data types. 
        a) Folder where extracted spot centroid position is stored for batch processing
        b) single .txt file for spot extraction in single channel
        c) numpy arrays with spot info  
    s: pixel size for segmentation mask, default to 0.92µm in x, y and 0.84µm in z. 
    
    """
    if type(spot_dir)==str:
        if os.path.isdir(spot_dir):
            fx = sorted(glob(spot_dir+"/*.txt"))
            lb_id = np.unique(lb[lb != 0])
            z, y, x = lb.shape
            count = pd.DataFrame(np.empty([len(lb_id), 0]), index=lb_id)
            for f in fx:
                print("Load:", f)
                r = os.path.basename(f).split('/')[-1]
                r = r.split('.')[0]
                spot = np.loadtxt(f, delimiter=',')
                n = len(spot)
                rounded_spot = np.round(spot[:, :3]/s).astype('int')
                df = pd.DataFrame(np.zeros([len(lb_id), 1]),
                                index=lb_id, columns=['count'])

                for i in range(0, n):
                    if np.any(np.isnan(spot[i,:3])):
                        print('NaN found in {} line# {}'.format(f, i+1))
                    else:
                        if np.any(spot[i,:3]<0) or np.all(np.greater(rounded_spot[i], [x, y, z])):
                            if verbose:
                                print('Point outside of fixed image found in {} line# {}'.format(f, i+1))
                        else:
                            try:
                                # if all non-rounded coord are valid values (none is NaN)
                                Coord = np.minimum(rounded_spot[i], [x, y, z])
                                idx = lb[Coord[2]-1, Coord[1]-1, Coord[0]-1]
                                if idx > 0 and idx <= len(lb_id):
                                    df.loc[idx, 'count'] = df.loc[idx, 'count']+1
                            except Exception as e:
                                print('Unexpected error in {} line# {}: {}'.format(f, i+1, e))
                count.loc[:, r] = df.to_numpy()
        else:
            lb_id = np.unique(lb[lb != 0])
            z, y, x = lb.shape
            count = pd.DataFrame(np.empty([len(lb_id), 0]), index=lb_id)
            print("Load:", spot_dir)
            r = os.path.basename(spot_dir).split('/')[-1]
            r = r.split('.')[0]
            spot = np.loadtxt(spot_dir, delimiter=',')
            n = len(spot)
            print(n)
            rounded_spot = np.round(spot[:, :3]/s).astype('int')
            df = pd.DataFrame(np.zeros([len(lb_id), 1]),
                            index=lb_id, columns=['count'])

            for i in range(0, n):
                if np.any(np.isnan(spot[i,:3])):
                    print('NaN found in {} line# {}'.format(f, i+1))
                else:
                    if np.any(spot[i,:3]<0) or np.all(np.greater(rounded_spot[i], [x, y, z])):
                        if verbose:
                            print('Point outside of fixed image found in line# {}'.format(i+1))
                            print(rounded_spot[i], (x,y,z))
                    else:
                        try:
                            # if all non-rounded coord are valid values (none is NaN)
                            Coord = np.minimum(rounded_spot[i], [x, y, z])
                            idx = lb[Coord[2]-1, Coord[1]-1, Coord[0]-1]
                            if idx > 0 and idx <= len(lb_id):
                                df.loc[idx, 'count'] = df.loc[idx, 'count']+1
                        except Exception as e:
                            print('Unexpected error in line# {}: {}'.format(i+1, e))
            count.loc[:, r] = df.to_numpy()
    else:
        lb_id = np.unique(lb[lb != 0])
        z, y, x = lb.shape
        count = pd.DataFrame(np.empty([len(lb_id), 0]), index=lb_id)
        spot = spot_dir.copy()
        n = len(spot)
        rounded_spot = np.round(spot[:, :3]/s).astype('int')
        df = pd.DataFrame(np.zeros([len(lb_id), 1]),
                        index=lb_id, columns=['count'])

        for i in range(0, n):
            if np.any(np.isnan(spot[i,:3])):
                print('NaN found in line# {}'.format(i+1))
            else:
                if np.any(spot[i,:3]<0) or np.all(np.greater(rounded_spot[i], [x, y, z])):
                    if verbose:
                        print('Point outside of fixed image found in line# {}'.format(i+1))
                else:
                    try:
                        # if all non-rounded coord are valid values (none is NaN)
                        Coord = np.minimum(rounded_spot[i], [x, y, z])
                        idx = lb[Coord[2]-1, Coord[1]-1, Coord[0]-1]
                        if idx > 0 and idx <= len(lb_id):
                            df.loc[idx, 'count'] = df.loc[idx, 'count']+1
                    except Exception as e:
                        print('Unexpected error in line# {}: {}'.format(i+1, e))
        count.loc[:, 'spot count'] = df.to_numpy()

    return(count)

def rm_lipofuscin(channel_1, channel_2, radius=0.69):
    """
    
    returns real FISH spots in channel_1 and channel_2 after removing identified lipofuscin spots.
    Autofluorescence lipofuscin spots are identified using two FISH channels. 
    Spots appearing in both channels at the same position are identified as autofluorescence spots. 
    
    channel_1: spots detected from first FISH channel (4 columns with x, y, z position (µm) and spot integrated intensity) 
    channel_2: spots detected from second FISH channel (4 columns with x, y, z position (µm) and spot integrated intensity) 
    radius: maximum distance between two spots, in µm. 
    """
    neighbor_radius   = radius
    kdtree_c0 = cKDTree(channel_1[:,:3])
    kdtree_c1 = cKDTree(channel_2[:,:3])
    neighbors = kdtree_c0.query_ball_tree(kdtree_c1, neighbor_radius)

    no_neighbors = 0
    one_neighbor = 0
    more_neighbors = 0
    max_neighbors = 0
    for nnn in neighbors:
        if len(nnn) == 0: no_neighbors += 1
        if len(nnn) == 1: one_neighbor += 1
        if len(nnn) > 1 : more_neighbors += 1
        if len(nnn) > max_neighbors: max_neighbors = len(nnn)

    # print(no_neighbors, one_neighbor, more_neighbors, max_neighbors)

    neighbors_num = np.array([len(x) for x in neighbors]).sum()
    pAind = np.empty(neighbors_num, dtype=np.uint32)
    pBind = np.empty(neighbors_num, dtype=np.uint32)

    p_ind = 0
    for c0_ind, nnn in enumerate(neighbors):
        if len(nnn) == 0:
            continue
        for c1_ind in nnn:
            pAind[p_ind]  = c0_ind
            pBind[p_ind]  = c1_ind
            p_ind += 1

    lipo_c0 = channel_1[pAind]
    lipo_c1 = channel_2[pBind]

    true_pos_c0 = np.delete(channel_1, pAind, axis=0)
    true_pos_c1 = np.delete(channel_2, pBind, axis=0)

    return (true_pos_c0, true_pos_c1, pAind, pBind)


def spot_counts_v2(lb, spot_dir, s, 
                   verbose=True,
                   remove_emptymask=True,
                   ):
    """
    Returns spot counts for each ROI. 
    
    lb: segmenntation mask (integer 3d image/matrix)
    spot_dir: Accepts 3 different data types. 
        a) Folder where extracted spot centroid position is stored for batch processing
        b) single .txt file for spot extraction in single channel
        c) numpy arrays with spot info  
    s: pixel size for segmentation mask, default to 
                       s=[0.92,0.92,0.84], 
    
    """
    # directory
    if isinstance(spot_dir, str) and os.path.isdir(spot_dir):
        lb_id = np.unique(lb[lb!=0]) # exclude 0
        lb_id = np.hstack([[0], lb_id]) # include 0

        counts_all = pd.DataFrame(index=lb_id)
        fx = sorted(glob(spot_dir+"/spots_*.txt"))
        for f in fx:
            print("Load:", f)
            r_c = os.path.basename(f).split('.')[0]
            spot = np.loadtxt(f, delimiter=',')
            res = spot_counts_worker(lb, spot, lb_id=lb_id, 
                                     remove_emptymask=remove_emptymask, 
                                     verbose=verbose,
                                     )
            counts_all[r_c] = res 
        
        if remove_emptymask:
            counts_all = counts_all.iloc[1:]
        return counts_all

    # txt file
    if isinstance(spot_dir, str) and os.path.isfile(spot_dir):
        print("Load:", spot_dir)
        r = os.path.basename(f).split('.')[0]
        spot = np.loadtxt(spot_dir, delimiter=',')
        res = spot_counts_worker(lb, spot, 
                                remove_emptymask=remove_emptymask, 
                                verbose=verbose,
                                )
        return res

    # numpy array
    if isinstance(spot_dir, np.ndarray):
        res = spot_counts_worker(lb, spot, 
                                remove_emptymask=remove_emptymask, 
                                verbose=verbose,
                                )
        return res

def spot_counts_worker(lb, spots, s, 
                       lb_id=None,
                       remove_noncell=True, 
                       selected_roi_list=None,
                       ):
    """
    Returns spot counts for each ROI. 
    
    lb: segmenntation mask (integer 3d image/matrix)
    spots: spots info (numpy array)
    s: pixel size for segmentation mask, default to 
                       s=[0.92,0.92,0.84], 
    
    """
    # numpy array
    z, y, x = lb.shape
    if lb_id is None:
        lb_id = np.unique(lb[lb!=0]) # exclude 0
        lb_id = np.hstack([[0], lb_id]) # include 0
    res = pd.Series(np.zeros(len(lb_id)), index=lb_id) # counts per roi (including 0) 

    # spot info
    # remove nan
    n = len(spots)
    spot = spots[~np.any(np.isnan(spots[:,:3]), axis=1),:3]
    print(f"removed {n-len(spot)} spots due to nan")

    # um to pixel
    spot = np.round(spot[:,:3]/s).astype('int')
    spot = spot-1  # why??? shift everything by 1 pixel - should be fine

    # remove outside range
    spot = spot[~np.any(spot<0, axis=1)]
    spot = spot[~(spot[:,0]>=x)]
    spot = spot[~(spot[:,1]>=y)]
    spot = spot[~(spot[:,2]>=z)]
    print(f"{len(spot):,}/{n:,} spots in range {(x,y,z)}")

    # get index and counts
    idx = lb[spot[:,2], spot[:,1], spot[:,0]]
    rois, counts = np.unique(idx, return_counts=True)
    res.loc[rois] = counts

    # remove noncell
    if remove_noncell and lb_id[0] == 0:
        res = res.iloc[1:] # 0 means not inside any mask
    
    # report only those selected cells
    if selected_roi_list is not None:
        res = res.loc[selected_roi_list]
    return res

def remove_bleed_thru_spots(ref_dots, query_dots, epsilon=3):
    """query dots that are within epsilon distance of any ref_dots will be flagged
    epsilon distance has the unit of query/ref dots, usually in um (post-expansion)
    """
    from scipy.spatial import cKDTree
    index = cKDTree(ref_dots[:,:3])
    nn_dists, nn_ids = index.query(query_dots[:,:3], k=1)
    
    cond = nn_dists < epsilon
    print(f"{cond.sum()}/{len(cond)} = {100*cond.sum()/len(cond):.1f}% removed")

    return query_dots[~cond], query_dots[cond]

def get_unit_spot_intn(path_spot, vox, xymin=250, xymax=1500, zmin=150, zmax=650):
    """estimate unit spot intensity
    
    vox=[0.92,0.92,0.84]
    """
    # (x, y, z, I)
    spot = np.loadtxt(path_spot, delimiter=',')
    
    ## remove spots on edges (eliminate false detection)
    spot[:,:3]=spot[:,:3]/vox  # convert from um (2x) to pixel unit
    spot = spot[np.logical_and(spot[:,0]<=xymax, spot[:,0]>xymin)]
    spot = spot[np.logical_and(spot[:,1]<=xymax, spot[:,1]>xymin)]
    spot = spot[np.logical_and(spot[:,2]<= zmax, spot[:,2]> zmin)]   
    
    ## assign the most frequent intensity as the single-spot-intensity
    spot_int = spot[:,3]
    spot_int = spot_int[spot_int!=-8.0] # ???
    n,b=np.histogram(spot_int, bins=5000)
    unit_intn = b[np.argmax(n)]
    
    return unit_intn

def get_spot_counts_from_intn(path_intn, path_spot, roi_meta, vox):
    """estimate spot_counts from cell_intensities; estimate unit-spot intensity first
    """
    unit_intn = get_unit_spot_intn(path_spot, vox) # get unit intn
    
    cell_int = pd.read_csv(path_intn, sep=',', index_col=0)
    cell_int = cell_int.reindex(roi_meta.index) ## only include intact ROIs###

    vec_mean = cell_int['mean_intensity'].values
    vec_area = roi_meta['area'].values

    # background - this is to assume most cells do not express this gene
    n,b = np.histogram(vec_mean, bins=1000)
    bg = b[np.argmax(n)]    
    
    # count
    vec_count = np.clip(vec_mean - bg, 0, None)*vec_area/unit_intn
    return vec_count