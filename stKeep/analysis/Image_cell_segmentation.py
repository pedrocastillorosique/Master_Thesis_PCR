# -*- coding: utf-8 -*-
"""

@author: chunman zuo
@modified_by: Pedro Castillo
"""

import numpy as np
import pandas as pd
import json
import cv2
from pathlib import Path
import os
import time
import sys
import tifffile as tiff

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import stKeep

print('Start processing cell segmentation')

start =  time.time()
parser =  stKeep.parameter_setting()
args = parser.parse_args()
imageSeg_dir    =  args.inputPath

image_loc_in    = pd.read_csv(imageSeg_dir + "spatial/tissue_positions_list.csv", header= None)
image_loc_in.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
image_loc_in = image_loc_in.set_index('barcode')

sh_fname        = os.path.join(imageSeg_dir, f"spatial/{args.jsonFile}")
scale_fname     = os.path.join(imageSeg_dir, "spatial/scalefactors_json.json")

image = cv2.imread(imageSeg_dir + "/spatial/tissue_hires_image.png")

with open(sh_fname, 'r') as f:
    sh_json = json.load(f)

with open(scale_fname, 'r') as f:
    scale_json = json.load(f)


if type(image_loc_in['in_tissue'].iloc[1]) == int:
    image_loc_in = image_loc_in.loc[image_loc_in['in_tissue'] == 1, ['pxl_row_in_fullres', 'pxl_col_in_fullres']].astype(float) * scale_json['tissue_hires_scalef']
else: 
    image_loc_in = image_loc_in.loc[image_loc_in['in_tissue'] == "1", ['pxl_row_in_fullres', 'pxl_col_in_fullres']].astype(float) * scale_json['tissue_hires_scalef']
    
check_or_not    = [False] * len(image_loc_in.index)
cell_type_dict  = {}	
region_pro      = 1e-5
spot_size       = int((scale_json['spot_diameter_fullres'] * scale_json['tissue_hires_scalef']) / 2) # Size of the spot in the original size image

if args.annotation_method == 'labelme':
    shape = (sh_json['imageWidth'], sh_json['imageHeight'])

    for sh in sh_json['shapes']:
        geom   = np.squeeze(np.array(sh['points']))
        pts    = geom.astype(int)
        mask   = np.zeros(shape, dtype= int)
        cv2.fillPoly(mask, [pts], 1)
        temp_celltype = sh['label']
        print(temp_celltype)
        count_status  = -1

        for barcode, imagerow, imagecol in zip(image_loc_in.index, image_loc_in["pxl_row_in_fullres"], image_loc_in["pxl_col_in_fullres"]):
            count_status   =+ 1
            imagerow_down  = imagerow - spot_size / 2
            imagerow_up    = imagerow + spot_size / 2
            imagecol_left  = imagecol - spot_size / 2
            imagecol_right = imagecol + spot_size / 2
            spot  = np.array([[imagecol_left, imagerow_up],[imagecol_left, imagerow_down], 
                            [imagecol_right,imagerow_down], [imagecol_right,imagerow_up]])
            pts   = spot.astype(int)
            mask1 = np.zeros(shape)
            cv2.fillPoly(mask1, [pts], 1)
            mask1 = mask1.astype(int)
            mask2 = mask + mask1

            if np.sum(mask2>1)/np.sum(mask1>0) > region_pro:
                if check_or_not[count_status] is False:
                    if cell_type_dict.__contains__(temp_celltype):
                        temp_sopt_list = cell_type_dict[temp_celltype]
                        temp_sopt_list.append(barcode)
                    else:
                        temp_spo_l = []
                        temp_spo_l.append(barcode)
                        cell_type_dict[temp_celltype] = temp_spo_l
                    check_or_not[count_status] = True

elif args.annotation_method == 'QuPath':
    shape = image.shape[:2]

    for idx, feature in enumerate(sh_json):
        props = feature["properties"]
        geometry = feature["geometry"]

        temp_celltype = props.get("name", "unknown")
        geom_type = geometry["type"]
        coords = geometry["coordinates"]

        print(f"\nAnnotation {idx + 1}")
        print(f"Label: {temp_celltype}")
        print(f"Geometry Type: {geom_type}")

        if geom_type == "Polygon":
            pts = np.array(coords[0], dtype=int)
        elif geom_type == "MultiPolygon":
            pts = np.array(coords[0][0], dtype=int)
        else:
            print("Geometry type not supported: ", geom_type)
            continue

        mask = np.zeros(shape, dtype=int)
        cv2.fillPoly(mask, [pts], 1)

        count_status = -1 # This will iterate through each barcode 

        for barcode, imagerow, imagecol in zip(image_loc_in.index, image_loc_in["pxl_row_in_fullres"], image_loc_in["pxl_col_in_fullres"]):
            count_status += 1

            mask1 = np.zeros(shape)
            cv2.circle(mask1, ((round(imagecol), round(imagerow))), spot_size , (1), -1)
            mask1 = mask1.astype(int)
            mask2 = mask + mask1
        
            # if count_status > 1000:
            #     tiff.imsave("mask1.tiff", mask1)
            #     tiff.imsave("mask2.tiff", mask2)
            #     tiff.imsave("mask.tiff", mask)
                
            if np.sum(mask2 > 1) / np.sum(mask1 > 0) > region_pro: # Check when the number of pixels cropped is greater than the 50% of the background 
                if not check_or_not[count_status]:
                    if temp_celltype in cell_type_dict:
                        cell_type_dict[temp_celltype].append(barcode)
                    else:
                        cell_type_dict[temp_celltype] = [barcode]
                    check_or_not[count_status] = True
else: 
      raise Exception("Non valid annotation method")


cell_types = []
cell_names = []
for cell_type in cell_type_dict.keys():
    temp_barcodes = cell_type_dict[cell_type]
    cell_names.extend(temp_barcodes)
    cell_types.extend([cell_type] * len(temp_barcodes))

remain_cells = list(set(image_loc_in.index.tolist()) - set(cell_names))
cell_names.extend(remain_cells)
cell_types.extend(['cluster'+str((len(cell_type_dict)+1))] * len(remain_cells)) # That correspond to those lower than 0.5 in the intersection spot/annotation
df = {'Cell': cell_names, 'Layer': cell_types}

data_frame   = pd.DataFrame(df)
unique_labels = data_frame['Layer'].unique()
label_to_cluster = {label: idx + 1 for idx, label in enumerate(unique_labels)}
data_frame['Cluster'] = data_frame['Layer'].map(label_to_cluster)
data_frame.to_csv(imageSeg_dir + '/Image_cell_segmentation.txt', sep='\t', index=False) 

duration     = time.time() - start
print('Finish training, total time is: ' + str(duration) + 's' )