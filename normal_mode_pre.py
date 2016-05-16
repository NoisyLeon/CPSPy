#!/usr/bin/env python
import cpspy
dfile=cpspy.DistFile();
dfile.addEqualDist(100, 100, 80)
dfile.write('dfile')
