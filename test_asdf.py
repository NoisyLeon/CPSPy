#!/usr/bin/env python
# -*- coding: utf-8 -*-
import multiprocessing
from asdf import AsdfFile
import numpy as np
import pyasdf

# def process(i):
#     print i
#     dset= pyasdf.ASDFDataSet("example"+str(i)+".asdf")
#     obspy.read()
#     # dset.add_wave
#     # tree = {
#     #     'author': 'John Doe'+str(i),
#     #     'my_array': np.random.rand(8, 8)
#     #   }
#     # ff = AsdfFile(tree)
#     # ff.write_to("example"+str(i)+".asdf")
# 
# p = multiprocessing.Pool(processes=6)
# p.map(process, xrange(10))

from asdf import AsdfFile
import numpy as np

tree = {
  'author': 'John Doe',
  'my_array': np.random.rand(8, 8)
}
ff = AsdfFile(tree)
ff.write_to("example.asdf")