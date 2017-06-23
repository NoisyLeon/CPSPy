import vmodel
import copy
import matplotlib.pyplot as plt
import numpy as np
model=vmodel.Model1d(modelindex=2)

model.read_axisem_bm('/home/leon/code/axisem/SOLVER/MESHES/prem_aniso_10s/1dmodel_axisem.bm')