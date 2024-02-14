import xml.etree.ElementTree as ET
import numpy as np
import glob
import os
import json
import shutil


def generate_xml_s1(f_src, f_dst):
    """
    """
    tree = ET.parse(f_src)
    root = tree.getroot()
    
    ## show
    # get viewSetups - voxelSize and image size
    for viewsetup in root[1][1].findall('ViewSetup'):
        print(viewsetup.find('voxelSize').find('size').text)
        print(viewsetup.find('size').text)
    
    # get viewRegistrations - ViewTransform
    for viewreg in root.find('ViewRegistrations').findall('ViewRegistration'):
        for viewtrans in viewreg.findall('ViewTransform'):
            # print(viewtrans.find('Name').text)
            print(viewtrans.find('affine').text)
            
    print('---')
    
    ## modify
    for viewsetup in root[1][1].findall('ViewSetup'):
        a = viewsetup.find('voxelSize').find('size').text
        b = viewsetup.find('size').text
    
        a = np.array(a.split(' ')).astype(float)
        a[:2] = a[:2]*2
        a = " ".join(a.astype(str).tolist())
        
        b = np.array(b.split(' ')).astype(int)
        b[:2] = b[:2]/2
        b = " ".join(b.astype(str).tolist())
    
        print(a)    
        print(b)    
        viewsetup.find('voxelSize').find('size').text = a
        viewsetup.find('size').text = b
        
    for viewreg in root.find('ViewRegistrations').findall('ViewRegistration'):
        for viewtrans in viewreg.findall('ViewTransform'):
            transform_type = viewtrans.find('Name').text
            transform_vals = viewtrans.find('affine').text
            a = transform_vals
            a = np.array(a.split(' ')).astype(float)
            
            if transform_type == 'Translation':
                a[3] = a[3]/2
                a[7] = a[7]/2
                
            elif transform_type == 'calibration':
                a[10] = a[10]/2
            else:
                raise ValueError('transform_type should be Translation or calibration only')
                
            a = " ".join(a.astype(str).tolist())
            print(a)    
            viewtrans.find('affine').text = a
            
    tree.write(f_dst, encoding='utf-8')
    
def generate_setup_attributes_s0_s1(f_n5):
    """
    # update attributes
    # 1. if attributes_s0 doesn't exist - copy from attributes
    # 2. use attrbutes_s0 to generate attributes_s1
    # 3. copy attributes_s1 as attributes
    """
    
    for view_setup in sorted(glob.glob(os.path.join(f_n5, 'setup*'))):
        f_attr    = os.path.join(view_setup, 'attributes.json')
        f_attr_s0 = os.path.join(view_setup, 'attributes_s0.json')
        f_attr_s1 = os.path.join(view_setup, 'attributes_s1.json')
        
        if not os.path.isfile(f_attr):
            raise ValueError("attributes.json doesn't exist")
            
        if not os.path.isfile(f_attr_s0):
            print('copy attributes.json as attributes_s0.json')
            shutil.copyfile(f_attr, f_attr_s0)
            
        if not os.path.isfile(f_attr_s1):
            print('generate attributes_s1.json from attributes_s0.json')
            with open(f_attr_s0, 'r') as fi:
                data = json.load(fi)
            
                a = data['downsamplingFactors']
                print('befor:', a)
                a = np.array(a).astype(int)
                a[:,:2] = (a[:,:2]/2).astype(int)
                a = a.tolist()
                print('after:', a)
                data['downsamplingFactors'] = a
                
            with open(f_attr_s1, 'w') as fo:
                json.dump(data, fo)
    return

def set_attributes(f_n5, level):
    """
    """
    for view_setup in sorted(glob.glob(os.path.join(f_n5, 'setup*'))):
        f_attr    = os.path.join(view_setup, 'attributes.json')
        f_attr_s0 = os.path.join(view_setup, 'attributes_s0.json')
        f_attr_s1 = os.path.join(view_setup, 'attributes_s1.json')

        if level == 's0':
            f_src = f_attr_s0
        elif level == 's1':
            f_src = f_attr_s1
        else:
            raise ValueError('level must be s0 or s1')
            
        if os.path.isfile(f_attr_s0) and os.path.isfile(f_attr_s1):
            if os.path.isfile(f_attr):
                os.remove(f_attr)
            shutil.copyfile(f_src, f_attr)
        else:
            raise ValueError('attributes_s0.json and attributes_s1.json must both exist')

if __name__ == '__main__':
    
    f_n5      = '/u/home/f/f7xiesnm/project-zipursky/data/00_tiled/sparse06_r1.n5/'
    f_xml_old = '/u/home/f/f7xiesnm/project-zipursky/data/00_tiled/lt185_r2.xml'
    f_xml_new = '/u/home/f/f7xiesnm/project-zipursky/data/00_tiled/lt185_r2_autos1.xml'

    generate_xml_s1(f_xml_old, f_xml_new)
    generate_setup_attributes_s0_s1(f_n5)
    set_attributes(f_n5, 's1')