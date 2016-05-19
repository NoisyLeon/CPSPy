import os, shutil
import numpy as np
import obspy
import time
import obspy.geodetics

class StaInfo(object):
    """
    An object contains a station information several methods for station related analysis.
    -----------------------------------------------------------------------------------------------------
    General Parameters:
    stacode      - station name
    network     - network
    chan          - channels for analysis
    lon, lat, z           - position for station
    distance    - epicentral distance
    -----------------------------------------------------------------------------------------------------
    """
    def __init__(self, stacode=None, network='CPS', lon=None, lat=None, z=0, distance=None):

        self.stacode=stacode
        self.network=network
        self.lon=lon
        self.lat=lat
        self.z=z
        self.distance=distance
        return
    
class StaLst(object):
    """
    An object contains a station list(a list of StaInfo object) information several methods for station list related analysis.
        stations: list of StaInfo
    """
    def __init__(self,stations=None):
        self.stations=[]
        if isinstance(stations, StaInfo):
            stations = [stations]
        if stations:
            self.stations.extend(stations)

    def __add__(self, other):
        """
        Add two StaLst with self += other.
        """
        if isinstance(other, StaInfo):
            other = StaLst([other])
        if not isinstance(other, StaLst):
            raise TypeError
        stations = self.stations + other.stations
        return self.__class__(stations=stations)

    def __len__(self):
        """
        Return the number of Traces in the StaLst object.
        """
        return len(self.stations)

    def __getitem__(self, index):
        """
        __getitem__ method of obspy.Stream objects.
        :return: Trace objects
        """
        if isinstance(index, slice):
            return self.__class__(stations=self.stations.__getitem__(index))
        else:
            return self.stations.__getitem__(index)

    def append(self, station):
        """
        Append a single StaInfo object to the current StaLst object.
        """
        if isinstance(station, StaInfo):
            self.stations.append(station)
        else:
            msg = 'Append only supports a single StaInfo object as an argument.'
            raise TypeError(msg)
        return self
    
    def ReadDistFile(self, distfname, network='CPS', lon0=0., lat0=0., lat=0.):
        """
        Read Distance file 
        DIST DT NPTS T0 VRED
        """
        InArr=np.loadtxt(distfname)
        distArr=InArr[:,0]
        dist, azi, baz= obspy.geodetics.gps2dist_azimuth(lon0, lat0, lat, 1.)
        dist=dist/1000.
        for i in np.arange(distArr.size):
            distance=distArr[i]
            lon=distance/dist
            stacode='B%03d' %(i+1)
            self.stations.append( StaInfo(stacode=stacode, network=network, lon=lon, lat=lat, z=0, distance=distance) )        
        return
    
    def ReadStaList(self, stafile):
        """
        Read Sation List from a txt file
        stacode network lon lat z distance
        """        
        with open(stafile, 'r') as f:
            Sta=[]
            for lines in f.readlines():
                lines=lines.split()
                stacode=lines[0]
                network=lines[1]
                lon=float(lines[2])
                lat=float(lines[3])
                z=float(lines[4])
                distance=float(lines[5])
                netsta=network+'.'+stacode
                if Sta.__contains__(netsta):
                    index=Sta.index(netsta)
                    if abs(self[index].lon - lon) >0.01 and abs(self[index].lat - lat) and abs(self[index].z-z)>0.01:
                        raise ValueError('Incompatible Station Location:' + netsta+' in Station List!')
                    else:
                        print 'Warning: Repeated Station:' +netsta+' in Station List!'
                        continue
                Sta.append(netsta)
                self.append(StaInfo (stacode=stacode, network=network, lon=lon, lat=lat, z=z, distance=distance) )
        return
    
    def WriteStaList(self, stafile):
        """
        Read Sation List from a txt file
        stacode network lon lat z distance
        """
        with open(stafile, 'w') as f:
            for sta in self.stations:
                tempstr='%s %s %g %g %g %g  \n' %( sta.stacode, sta.network, sta.lon, sta.lat, sta.z, sta.distance)
                f.writelines(tempstr)
        # f.close()
        return
    
    def AddSingle(self, lon, lat, stacode, network='CPS',  z=0., lon0=0., lat0=0.):
        dist, azi, baz= obspy.geodetics.gps2dist_azimuth(lon0, lat0, lat, lon)
        dist=dist/1000.
        self.append(StaInfo (stacode=stacode, network=network, lon=lon, lat=lat, z=z, distance=distance) )
        return
    
    def GetInventory(self, outfname=None, chans=['ZVF'],  source='CUCPS'):
        chandict={'ZDD': '01', 'RDD': '02', 'ZDS': '03', 'RDS': '04', 'TDS': '05', 'ZSS': '06', 'RSS': '07', 'TSS': '08', 'ZEX': '09', 'REX': '10', 'ZVF': '11', 'RVF': '12',
                      'ZHF': '13', 'RHF': '14', 'THF': '15'}
        stations=[]
        total_number_of_channels=len(chans)
        site=obspy.core.inventory.util.Site(name='01')
        creation_date=obspy.core.utcdatetime.UTCDateTime(0)
        for sta in self.stations:
            channels=[]
            for chan in chans:
                channel=obspy.core.inventory.channel.Channel(code=chan, location_code=chandict[chan], latitude=sta.lat, longitude=sta.lon, elevation=sta.z, depth=0.0)
                channels.append(channel)
            station=obspy.core.inventory.station.Station(code=sta.stacode, latitude=sta.lat, longitude=sta.lon, elevation=sta.z,
                    site=site, channels=channels, total_number_of_channels = total_number_of_channels, creation_date = creation_date)
            stations.append(station)
        network=obspy.core.inventory.network.Network(code=sta.network, stations=stations)
        networks=[network]
        inv=obspy.core.inventory.inventory.Inventory(networks=networks, source=source)
        if outfname!=None:
            inv.write(outfname, format='stationxml')
        return inv

    
