import os
import pandas as pd
import numpy as np
import shutil
import csv

def read_positions_csv(filename):
    f = open(filename)
    csvreader = csv.reader(f)
    positionsx = []
    positionsy = []
    positionsz = []

    for line in csvreader:
        #print(line)
        try:
            if "." in line[0]:
                positionsx.append(float(line[0]))
                positionsy.append(float(line[1]))
                positionsz.append(float(line[2]))
        except:
            pass
    positions=np.array([positionsx,positionsy,positionsz])
    #print(positions)
    return positions

def write_storm_csv(
    filename, points, coords=["x", "y", "z"], split_by=None, **csv_args
):
    """
    Saves a point cloud array in the storm format
    """
    columns = ["frame", "x [nm]", "y [nm]", "z [nm]", "uncertainty_xy", "uncertainty_z"]
    points = points.dropna()
    storm_points = pd.DataFrame(np.zeros((points.shape[0], 6)), columns=columns)
    storm_points[["x [nm]", "y [nm]", "z [nm]"]] = points[coords].values
    storm_points["frame"] = 1
    storm_points[["uncertainty_xy", "uncertainty_z"]] = 2.1
    # tab separated values are faster and more portable than excel
    if split_by is None:
        if not filename.endswith(".csv"):
            filename = filename + ".csv"
        storm_points.to_csv(filename, **csv_args)
    elif split_by in points.columns():
        storm_points[split_by] = points[split_by]
        # separated files by the column split_by
        storm_points.groupby(split_by).apply(
            lambda df: df.to_csv(
                "{}_{}.csv".format(filename, df[split_by].iloc[0]), **csv_args
            )
        )

def eptm_to_csv(ept, path=None):
    """
    Saves epithelium data in separates csv files
    """
    if path==None:
        path = os.getcwd()
        csv_path = str(path)+"/csv/"
    else:
        csv_path= str(path)
    if os.path.exists(csv_path)==False:
        os.mkdir(csv_path)
        print ("creation of the directory %s " % csv_path)
    else:
        shutil.rmtree(csv_path)
        os.mkdir(csv_path)	
    for key in ept.data_names:
        df=getattr(ept, "{}_df".format(key))
        df.to_csv(('%s%s.csv'% (csv_path, key)), index=False)
    print("Properties "+ str(ept.data_names)+ " are saved in csv format")
