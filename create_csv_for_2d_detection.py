import pandas as pd
import pydicom
import os
import glob
import numpy as np
from tqdm import tqdm
import pydicom
import math

# input local csv
LOCAL_CSV = "liver_ct_mass_export_271020.csv"
# output formatted csv
DETECTION_CSV = "liver_ct_mass_detection_interp_data_export_271020.csv"
# input raw data folder
RAW_LIVER_CT_FOLDER = "/media/datnt/data/liver-ct-raw-mass/"


def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1], idx-1
    else:
        return array[idx], idx

def num_slices_between(a,b):
    return max(abs(a-b)-1, 0)

df = pd.read_csv()

row_data = []
for enum, s in enumerate(tqdm(df.studyUid.unique())):
    if not os.path.exists(RAW_LIVER_CT_FOLDER + f"{s}/"):
        continue
    else:
        series_data = {}
        for f in glob.glob(RAW_LIVER_CT_FOLDER + f"{s}/*.*"):
            data = pydicom.read_file(f, force=True)
            slice_pos = data.get('ImagePositionPatient')
            spacing = data.get('PixelSpacing')
            series_data[f.split('/')[-1][:-4]] = [float(slice_pos[0]), float(slice_pos[1]), float(slice_pos[2]), float(spacing[0]), float(spacing[1])]
        series_data = {k: v for k, v in sorted(series_data.items(), key=lambda item: item[1])}
        name_data = np.array([f for f in series_data.keys()])
        pos_data = np.array([f for f in series_data.values()])
        dfs = df[df["studyUid"] == s]
        for ts in dfs.createTimestamp.unique():
            dfsts = dfs[dfs["createTimestamp"] == ts]
            idx_with_box = []
            x_min_vals = []
            x_max_vals = []
            y_min_vals = []
            y_max_vals = []
            interpolated_x_min_vals = []
            interpolated_x_max_vals = []
            interpolated_y_min_vals = []
            interpolated_y_max_vals = []
            interpolated_idx_with_box = []
            for zp in dfsts.z_pos.unique():
                dfstszp = dfsts[dfsts["z_pos"] == zp]
                _, idx = find_nearest(pos_data[:,2], zp)
                x_vals = ((dfstszp.x_pos.values - pos_data[idx,0]) / pos_data[idx,3]).astype(np.int)
                y_vals = ((dfstszp.y_pos.values - pos_data[idx,1]) / pos_data[idx,4]).astype(np.int)                
                idx_with_box.append(idx)
                x_min_vals.append(min(x_vals))
                x_max_vals.append(max(x_vals))
                y_min_vals.append(min(y_vals))
                y_max_vals.append(max(y_vals))
            if len(idx_with_box) == 1:
                continue
            x_min_vals = [x_min_vals[i] for i in np.argsort(idx_with_box)]
            x_max_vals = [x_max_vals[i] for i in np.argsort(idx_with_box)]
            y_min_vals = [y_min_vals[i] for i in np.argsort(idx_with_box)]
            y_max_vals = [y_max_vals[i] for i in np.argsort(idx_with_box)]
            idx_with_box = [idx_with_box[i] for i in np.argsort(idx_with_box)]            
            for i in range(len(idx_with_box) - 1):
                num_pad_slice = num_slices_between(idx_with_box[i], idx_with_box[i+1])
                interpolated_idx_with_box.append(idx_with_box[i])
                interpolated_x_min_vals.append(x_min_vals[i])
                interpolated_x_max_vals.append(x_max_vals[i])
                interpolated_y_min_vals.append(y_min_vals[i])
                interpolated_y_max_vals.append(y_max_vals[i])
                for ii in range(num_pad_slice):
                    interpolated_idx_with_box.append(idx_with_box[i] + (ii+1) * (idx_with_box[i+1] - idx_with_box[i]) / (num_pad_slice+1))
                    interpolated_x_min_vals.append(x_min_vals[i] + (ii+1) * (x_min_vals[i+1] - x_min_vals[i]) / (num_pad_slice+1))
                    interpolated_x_max_vals.append(x_max_vals[i+1] + (ii+1) * (x_max_vals[i+1] - x_max_vals[i]) / (num_pad_slice+1))
                    interpolated_y_min_vals.append(y_min_vals[i] + (ii+1) * (y_min_vals[i+1] - y_min_vals[i]) / (num_pad_slice+1))
                    interpolated_y_max_vals.append(y_max_vals[i] + (ii+1) * (y_max_vals[i+1] - y_max_vals[i]) / (num_pad_slice+1))
            interpolated_idx_with_box.append(idx_with_box[i+1])
            interpolated_x_min_vals.append(x_min_vals[i+1])
            interpolated_x_max_vals.append(x_max_vals[i+1])
            interpolated_y_min_vals.append(y_min_vals[i+1])
            interpolated_y_max_vals.append(y_max_vals[i+1])
                
            for iii in range(len(interpolated_idx_with_box)):
                row_data.append([s, name_data[int(interpolated_idx_with_box[iii])], int(interpolated_x_min_vals[iii]), int(interpolated_x_max_vals[iii]), int(interpolated_y_min_vals[iii]), int(interpolated_y_max_vals[iii])])

pd.DataFrame(data=row_data, columns=["studyUid", "imageUid", "xmin", "xmax", "ymin", "ymax"]).drop_duplicates().to_csv(DETECTION_CSV, index=False)                