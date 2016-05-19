import stations

SLst=stations.StaLst()
SLst.ReadDistFile('dfile')
SLst.WriteStaList('stations.lst')

SLst2=stations.StaLst()
# SLst.ReadDistFile('dfile')
SLst2.ReadStaList('stations.lst')
inv=SLst.GetInventory()