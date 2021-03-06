#!/usr/bin/env python
# -*- conding: utf-8 -*-

import os
import numpy as np
import datetime
import configparser
from . import letkf

class LETKF_core(object):
    
    required_instance_vals = ['nLon', 'east', 'assimE', 'nLat', 'res', 'patch', 'assimS', 'west', 'assimW', 'north', 'assimN', 'ensMem', 'south', 'undef']
    instance_vals_num = len(required_instance_vals)

    def __init__(self,configPath=None):
        """
            initial settings.
        """
        if type(configPath) == str and os.path.exists(configPath):
            # Read configuraion file.
            print("Read variables from configuration...")
            config = configparser.ConfigParser()
            config.read(configPath)

            self.assimN = float(config.get("assimilation","assimN"))
            self.assimS = float(config.get("assimilation","assimS"))
            self.assimE = float(config.get("assimilation","assimE"))
            self.assimW = float(config.get("assimilation","assimW"))
            self.patch  = int(config.get("assimilation","patch"))
            self.ensMem = int(config.get("assimilation","ensMem"))
            self.nLon   = int(config.get("model","nLon"))
            self.nLat   = int(config.get("model","nLat"))
            self.res    = float(config.get("model","res"))
            self.north  = float(config.get("model","north"))
            self.south  = float(config.get("model","south"))
            self.east   = float(config.get("model","east"))
            self.west   = float(config.get("model","west"))
            self.undef  = float(config.get("observation","undef"))

            print("############Check instance variables############")
            self.__showProperties()
            print("##############")


    def letkf(self,ensembles,observation,obserr,ocean,excGrids):
        """
        Data Assimilation with Local Ensemble Transformed Kalman Filter
        inputs:
            ensembles: numpy.ndarray([nLat,nLon,eNum]): ensemble simulation
            observation: numpy.ndarray([nLat,nLon]): gridded observation with observed or undef values
            ocean: numpy.ndarray([nLat,nLon]): gridded ocean mask
            excGrids: numpy.ndarray([nLat,nLon]): grids to exclude
        """
        # check all instance variables are set.
        self.__checkInstanceVals()
        # check shapes are correspond with instance correctly
        eNum = ensembles.shape[0]
        if eNum != self.ensMem:
            raise IOError("Specified emsemble member %d is not match with the passed array shape %d." % (self.ensMem,eNum))
        for data in [ensembles,observation,ocean,excGrids]:
            self.__checkLatLonShapes(data)
        # main letkf
        xa = letkf.letkf(ensembles,observation,obserr,ocean,excGrids,self.patch,self.ensMem,self.assimE,self.assimW,self.assimN,self.assimS,self.east,self.west,self.north,self.south,self.res,self.undef)
    
        return xa


    def __checkInstanceVals(self):
        
        keys = self.__dict__.keys()
        if len(keys) != self.instance_vals_num:
            nonDefined = set(self.required_instance_vals) - set(keys)  
            raise IOError("Not all instance variables defined. %s" % str(list(nonDefined)))


    def __checkLatLonShapes(self,array):

        if len(array.shape) == 2:
            lat_shape = array.shape[0]
            lon_shape = array.shape[1]
        elif len(array.shape) == 3:
            lat_shape = array.shape[1]
            lon_shape = array.shape[2]
        if lat_shape != self.nLat:
            raise IOError("Specified latitude number %d is not match with the passed array shape %d." % (self.nLat,lat_shape))
        if lon_shape != self.nLon:
            raise IOError("Specified longitude number %d is not match with the passed array shape %d." % (self.nLon,lon_shape))


    def __showProperties(self):
        
        for key in self.__dict__.keys():
            print("%s:%s" % (key,self.__dict__[key]))


if __name__ == "__main__":

    chunk = LETKF_core("./config.ini")
