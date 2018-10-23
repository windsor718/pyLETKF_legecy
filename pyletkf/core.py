#!/usr/bin/env python
# -*- conding: utf-8 -*-

import os
import numpy as np
import datetime
import configparser
import argparse
import subprocess
import src.extentions.sendSlack as slack
from multiprocessing import Pool
from multiprocessing import Process

class LETKF(object):
    
    def __init__(self,configPath="./config.ini"):
        """
            initial settings.
        """
        # Read configuraion file.
        self.config = configparser.SafeConfigParser()
        self.config.read(configPath)

        self.assimN = self.config.get("assimilation","assimN")
        self.assimS = self.config.get("assimilation","assimS")
        self.assimE = self.config.get("assimilation","assimE")
        self.assimW = self.config.get("assimilation","assimW")
        self.patch  = self.config.get("assimilation","patch")
        self.ensMem = self.config.get("assimilation","ensMem")
        self.errexp = self.config.get("assimilation","errExp")

        self.srcDir = self.config.get("model","srcDir")
        self.nLon   = self.config.get("model","nLon")
        self.nlat   = self.config.get("model","nLat")
        self.res    = self.config.get("model","res")

        self.sysDir = self.config.get("system","sysDir")

        self.slack  = self.config.getboolean("log","slackNotify")
        
        if self.slack:
            self.url = self.config.get("log","url")
            slackSender = slack.sendSlack(self.url)

        # Version
        string      = "LETKF-Python-API-Core:[%s]"
        version     = self.config.get("develop","version")
        developer   = self.config.get("develop","developer")
        lastUpdated = self.config.get("develop","lastUpdated")
        print "="*80
        print string % (version)
        print "Developers: %s" % (developer)
        print "Last updated: %s" % (lastUpdated)
        print "="*80


        # argparse option setting
        parser      = argparse.ArgumentParser(description="LETKF-Python-API-Core")
        parser.add_argument("-c","--compile",action="store_true",default=False,help="Compile Option; Use this to compile source codes.")
        self.pArgs  = parser.parse_args()


        #create namelist for the fortran code.
        self.createNameList()


    def main(self):
        """
            main controller
        """

        # <-- Preparation

        if self.pArgs.compile == True:
            print "compile option on."
            flag = compile_func()
           
        initial()

        # -->


        # <-- send notification if true

        if self.slack == True:
            user = "LETKF-Python-API-Core"
            text = "\"Data Assimilation is in progress.\""
            slack.progress(user,text)

        # -->


        # <-- Main Data Assimilation
        
        flag = self.dataAssim()

        # -->

        
        # <-- send notification if true

        if self.slack == True:
            user = "LETKF-Python-API-Core"
            text = "\"Data Assimilation is successfully finished.\""
            slack.success(user,text)

        # -->
        
        return 0

        
    def dataAssim(self):

        try:
            path = os.path.join(self.sysDir,"src/letkf/main")
            subprocess.check_call([path])
        except(subprocess.CalledProcessError):
            if self.slack == True:
                user = "LETKF-Python-API-Core"
                text = "\"Subprocess Error @ main\""
                slack.failed(user,text)
            raise subprocess.CalledProcessError

        return 0 


    def compileFunc(self):
        
        print "Compiling function ON."
        crntDir = os.getcwd()
        os.changedirs(os.path.join(sysDir,"src"))
        subprocess.call(["make","clean"])
        subprocess.call(["make","all"])
        os.changedirs(crntDir)

        return 0


    def createNameList(self):

        namelist    = ["assimN", "assimS", "assimE", "assimW", "ensMem",\
                       "errexp", "srcDir", "nLon",\
                       "nlat", "res"]
        namelistDct = {"assimN":self.assimN, "assimS":self.assimS, "assimE":self.assimE, "assimW":self.assimW, "ensMem":self.ensMem,\
                       "errexp":self.errexp, "srcDir":self.srcDir, "nLon":self.nLon,\
                       "nlat":self.nlat, "res":self.res}

        path = os.path.join(self.sysDir,"src/letkf/namelist.txt")
        with open(path,"w") as f:
            f.write("&input\n")
            for name in namelist:
                f.write("%s=%s\n"%(name,namelistDct[name]))
            f.write("/\n!\n")
        
        return 0


if __name__ == "__main__":

    chunk = LETKF()
    chunk.createNameList()
