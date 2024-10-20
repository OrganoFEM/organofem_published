#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import odbAccess
import sys
from os import path
import numpy as np
# Get odb file ================================================================
modelodb = str(sys.argv[-1]) # Last argument must be the path to get the odb file.
odbdir=path.dirname(modelodb)
odb=odbAccess.openOdb(modelodb)
# Settings ====================================================================
my_step='Step-1'
my_var='U'
## Check nodes outside instance  =======================================================
#instanced_node=[]
#for my_node in range(len(odb.steps[my_step].frames[0].fieldOutputs[my_var].values)):
#    if odb.steps[my_step].frames[0].fieldOutputs[my_var].values[my_node].instance  is not None:
#        instanced_node.append(my_node)
# Check dimension  ====================================================================
if len(odb.steps[my_step].frames[0].fieldOutputs['U'].componentLabels)==3:
    # Get displacement ====================================================================
    U=np.empty((3*len(odb.steps[my_step].frames[0].fieldOutputs[my_var].values),0),dtype='float64')
    for my_frame in range(len(odb.steps[my_step].frames)):
        Utmp=np.array([],dtype='float64')
        Vtmp=np.array([],dtype='float64')
        Wtmp=np.array([],dtype='float64')
        frame = odb.steps[my_step].frames[my_frame]
        disps = frame.fieldOutputs[my_var]
        for my_node in range(len(disps.values)):
            Utmp=np.append(Utmp,disps.values[my_node].dataDouble[0])
            Vtmp=np.append(Vtmp,disps.values[my_node].dataDouble[1])
            Wtmp=np.append(Wtmp,disps.values[my_node].dataDouble[2])
        U=np.column_stack([U,np.hstack([Utmp,Vtmp,Wtmp])])
    """ U shape ===================================================================
	        Step0  Step1  Step2  ...
	     N0  U       U      U    ...
	     N1  U       U      U    ...
	         ...    ...    ...
	     N0  V       V      V    ...
	     N1  V       V      V    ...
	         ...    ...    ...
	     N0  W       W      W    ...
	     N1  W       W      W    ...
	         ...    ...    ...
    """
elif len(odb.steps[my_step].frames[0].fieldOutputs['U'].componentLabels)==2:
    # Initialize Utmp and Vtmp ====================================================
    Utmp=np.zeros((len(odb.steps[my_step].frames[0].fieldOutputs[my_var].values),len(odb.steps[my_step].frames)))
    Vtmp=np.array(Utmp) # copy Utmp (and not markers of Utmp!)
    # Get displacement ============================================================
    for my_frame in range(len(odb.steps[my_step].frames)):
        frame = odb.steps[my_step].frames[my_frame]
        disps = frame.fieldOutputs[my_var]
        for my_node in range(len(disps.values)):
            Utmp[my_node,my_frame]=disps.values[my_node].dataDouble[0]
            Vtmp[my_node,my_frame]=disps.values[my_node].dataDouble[1]
    U=np.vstack([Utmp,Vtmp])
    """ U shape ===================================================================
            Step0  Step1  Step2  ...
         N0  U       U      U    ...
         N1  U       U      U    ...
             ...    ...    ...
         N0  V       V      V
         N1  V       V      V
             ...    ...    ...
	"""
else:
    raise ValueError('Dimension mismatch !')

# Write data in file ==========================================================
if int(sys.version[0]) == 3:
    np.save(path.join(odbdir,'disp-generate.npy'),U)
elif int(sys.version[0]) == 2:
    import pickle
    pickle.dump(U, open(path.join(odbdir,'disp-generate.npy'), 'wb'))