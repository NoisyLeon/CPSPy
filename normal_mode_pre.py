#!/usr/bin/env python
import cpsfile
dfile=cpsfile.DistFile();
dfile.addEqualDist(0, 10, 200)
dfile.write('dfile')
