# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
PCAS driver for software syncing the Nikon DS-Qi2 Detector with Oxford & SRI X-ray sources.
Adapted from varianSync.py. Requires sscan and scanProgress IOC's to be running.

__author__   =   Alireza Panna
__status__   =   Stable
__version__  =   1-3
__to-do__    =   Need to test usb-6501 to see how fast we can poll camera output lines and how long to wait after
                 sending trigger inputs to the camera. Currently after sending trigger we wait 25ms. Since the camera acquisition
                 starts at rising transition (0->1)  We effectively wait 50ms between exposures. 
__date__     =   4/1/2016
"""
"""
CHANGELOG:
05/04/2016       (AP) all input records are postfixed by _RBV, 
                      fixed minor bug in polling of FIRING_RBV record for xray, update to version 1-1, 
                      remove unused imports.
06/15/2016       (AP) fixed bug in DOC enum string. ON corresponds to value 1 and OFF corresponds to value 0. 
                      Added new record $(P):REPORT to write scan report and image statistics file for scans. 
06/26/2016       (AP) added all abort method (and corresponding record $(P):ABORT) for master abort of
                      all pvs. 
06/27/2016       (AP) Clean exit on keyboard interrupt (andrews method ported from cpisync), minor bugs were fixed.
                      Bug in statistics plot still needs fixing! (plot does not refresh)
07/28/2016       (AP) Fix bug in statistics plot. Now for each scan, past stats are destroyed by deleting the axis
                      object and creating a new one. Remove all SCAN_POSITIONER_PVS since this is not necessary for the
                      virtual sync circuit to function.
08/02/2016       (AP) Fix bug in reporting method. When the ioc is launched with $(P)REPORT set, during the first scan 
                      statistics weren't being saved. This is now fixed. Also fixed an off by one error in the calculation
                      of average/stddev/min/max shutter times for the camera. Scan report now saves proper accumulated 
                      x-ray exposure values for scans. It is derived from the x-ray driver pv's ($(P)MAS_RBV, $(P)MS_RBV). 
                      Removed setting value record for all pv's in this ioc during initialization, 
                      since autosave already restores previously saved values.     
08/04/2016       (AP) Added readback pv's related to acquisition timing, more debug prints added. 
08/08/2016       (AP) Set interactive plot off (ioff()), incase we are not running  daemon, the figure window should not
                      display
01/11/2017       (AP) Remove scanProgress dependency, fix bug in scan abort callback, scan report file title 
                      updated, exception added if dutycycle list is empty after scan, scan abort sequence now runs in a seperate
                      thread rather than from the callback function itself, implemented soft sync (rename hardSync() to sync()
                      as a fallback incase DAQ is disconnected or not working. Update version to 1-3.
"""
from pcaspy import SimpleServer, Driver
from epics import *
from PyDAQmx import *
import os, sys, datetime, psutil, time
import threading
import numpy as np
import matplotlib.pyplot as plt
import inspect

sys.path.append(os.path.realpath('../utils'))
import epicsApps

EXPERIMENT = 'HPFI:'
DAQ_NAME                      = 'cpiSync'
NIKON_EXPOSE                  = DAQ_NAME + '/port2/line0' # Output from DAQ to Qi2 requests exposure (0->1 starts exposure)
NIKON_TRIGGERREADY            = DAQ_NAME + '/port2/line1' # High when Qi2 ready to expose
NIKON_EXPOSEOUT               = DAQ_NAME + '/port2/line2' # High when it is exposing
NIKON_INPUTS                  = NIKON_TRIGGERREADY + "," + NIKON_EXPOSEOUT
# Main IOC records
XRAY_IOC                      = EXPERIMENT  + 'OXFORD:xray:'
SCAN_IOC                      = EXPERIMENT + 'SCAN:'
MOTOR_IOC                     = EXPERIMENT + 'KOHZU:'
DET_IOC                       = EXPERIMENT + 'Qi2:'
# Nikon Ks PV's
NIKON_ACQUIRE                 = PV(DET_IOC + 'cam1:Acquire', callback = True)
NIKON_ACQUIRE_TIME_RBV        = PV(DET_IOC + 'cam1:AcquireTime_RBV', callback = True)
NIKON_FULL_FILENAME_RBV       = PV(DET_IOC + 'TIFF1:FullFileName_RBV', callback = True)
NIKON_FILEPATH                = PV(DET_IOC + 'TIFF1:FilePath', callback = False)
NIKON_IMAGE_MODE              = PV(DET_IOC + 'cam1:ImageMode', callback = False)
NIKON_TRIGGER_MODE            = PV(DET_IOC + 'cam1:NikonTriggerMode', callback = False)
# Nikon AD statistics PV's
NIKON_STATS_MEAN              = PV(DET_IOC + 'Stats1:MeanValue_RBV', callback = False)
NIKON_STATS_SIGMA             = PV(DET_IOC + 'Stats1:Sigma_RBV', callback = False)
# sscan PV's
SCAN_DETECTOR_1               = PV(SCAN_IOC + 'scan1.T1PV', callback = False)
SCAN_ABORT                    = PV(SCAN_IOC + 'AbortScans.PROC', callback = True)

SCAN_1_NPTS                   = PV(SCAN_IOC + 'scan1.NPTS')
SCAN_2_NPTS                   = PV(SCAN_IOC + 'scan2.NPTS')  
SCAN_3_NPTS                   = PV(SCAN_IOC + 'scan3.NPTS')
SCAN_4_NPTS                   = PV(SCAN_IOC + 'scan4.NPTS')  

SCAN_1_BUSY                   = PV(SCAN_IOC + 'scan1.BUSY')
SCAN_2_BUSY                   = PV(SCAN_IOC + 'scan2.BUSY')
SCAN_3_BUSY                   = PV(SCAN_IOC + 'scan3.BUSY')
SCAN_4_BUSY                   = PV(SCAN_IOC + 'scan4.BUSY')
SCAN_DIM                      = PV(SCAN_IOC + 'ScanDim.VAL')

# Constant numpy arrays for setting digital outputs high or low on NI-DAQs
LOW  = np.zeros((1,), dtype=np.uint8)
HIGH = np.ones((1,), dtype=np.uint8)
# while loop and time.sleep poll time 1ms
POLL_TIME = 0.001
# If DOC is ON (1) save motor pv's.
MOTOR_IOC_LIST = [
                    MOTOR_IOC + 'm1',  MOTOR_IOC + 'm2',  MOTOR_IOC + 'm3', \
                    MOTOR_IOC + 'm4',  MOTOR_IOC + 'm5',  MOTOR_IOC + 'm6', \
                    MOTOR_IOC + 'm7',  MOTOR_IOC + 'm8',  MOTOR_IOC + 'm9', \
                    MOTOR_IOC + 'm10', MOTOR_IOC + 'm11', MOTOR_IOC + 'm12',\
                    MOTOR_IOC + 'm13', MOTOR_IOC + 'm14', MOTOR_IOC + 'm15',\
                    MOTOR_IOC + 'm16', MOTOR_IOC + 'm17', MOTOR_IOC + 'm18',\
                    MOTOR_IOC + 'm19'
                 ]
# pcas records
prefix = EXPERIMENT + 'NikonKsSync:'
pvdb = {

    'XSYNC'                     : {'type'  : 'enum',
                                   'enums' : ['NONE', 'SRI', 'OXFORD'],},
    'XSYNC_RBV'                 : {'type'  : 'enum',
                                   'enums' : ['NONE', 'SRI', 'OXFORD'],
                                   'scan'  : 1,},
    'TRIGGER'                   : {'asyn'  : True},
    'SEND_TRIGGER'              : {'asyn'  : True},
    'TRIGGER_READY_RBV'         : {'asyn'  : True},
    'EXPOSING_RBV'              : {'asyn'  : True},
    'DOC'                       : {'type'  : 'enum',
                                   'enums' : ['OFF', 'ON'],},      
    'REPORT'                    : {'type'  : 'enum',
                                   'enums' : ['OFF', 'ON'],},
    'ABORT'                     : {'asyn'  : True},
    'DUTYCYCLE_RBV'             : {'prec'  : 6,
                                   'unit'  : '%'},
    'ACQUIRE_RBV'               : {'prec'  : 6,
                                   'unit'  : 's'},
    'IDLE_RBV'                  : {'prec'  : 6,
                                   'unit'  : 's'},
    'TRIGGER_ACQUIRE_DELAY_RBV' : {'prec'  : 6,
                                   'unit'  : 's'},
    'SYNC_MODE_RBV'             : {'type'  : 'enum',
                                   'enums' : ['SOFT', 'HARD'],},
}
pvdb.update(epicsApps.pvdb)

class myDriver(Driver):
    """
    PCAS Driver
    """
    def  __init__(self):
        """
        Constructor to set initial parameters polling threads, callbacks etc
        """
        super(myDriver, self).__init__()
        # set high priority for this process
        self.setProcessPriority()
        self.iocStats()
        self.pulseDelayAvrg = []
        # callback PV's
        NIKON_FULL_FILENAME_RBV.add_callback(self.checkDoc)    
        SCAN_ABORT.add_callback(self.scanAbortCallback) 
        self.xrayReadyTime          = 0
        self.shutter_close          = 0
        self.shutter_open           = 0
        self.idleTime               = 0
        self.trigger_not_ready_rbv  = 0
        self.total_npts             = 0
        self.scanStartTime          = 0
        self.elapsedScanTime        = 0
        self.total_cpts             = 0
        self.dutyCycleList          = []
        self.shutterTimeList        = []
        self.mean_stats             = []
        self.sigma_stats            = []
        # Create the figure object
        plt.ioff()
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(1, 1, 1)
        # Set up DO task to trigger Qi2 expose
        self.written = int32()
        try:
            self.qi2TriggerTask = TaskHandle()
            DAQmxCreateTask("",byref(self.qi2TriggerTask))
            DAQmxCreateDOChan(self.qi2TriggerTask, NIKON_EXPOSE, "", DAQmx_Val_ChanForAllLines)
            DAQmxStartTask(self.qi2TriggerTask)     
            # Set up polling for the 2 Qi2 Inputs
            self.qi2PollTask = TaskHandle()
            DAQmxCreateTask("",byref(self.qi2PollTask))
            DAQmxCreateDIChan(self.qi2PollTask, NIKON_INPUTS , "", DAQmx_Val_ChanForAllLines)
            DAQmxStartTask(self.qi2PollTask)
            # Start the Qi2 polling thread
            self.xid = threading.Thread(target=self.pollQi2,args=())
            self.xid.start()
            self.daq_init_flag = 1
        except:
            self.daq_init_flag = 0
            print str(datetime.datetime.now())[:-3], "In function: " + inspect.stack()[0][3] + \
                  ", Failed to create DAQ task. Check to see if DAQ is connected"
            pass
        print self.daq_init_flag
        # Set scan detector PV to NikonSync hard trigger PV on init.
        SCAN_DETECTOR_1.put(EXPERIMENT + 'NikonKsSync:TRIGGER')
        # handle autosave related stuff
        epicsApps.buildRequestFiles(prefix, pvdb.keys(), os.getcwd())
        epicsApps.makeAutosaveFiles()
        print '############################################################################'
        print '## ADNIKONKS PCAS IOC Online $Date:' + str(datetime.datetime.now())[:-3]
        print '############################################################################'

    def setProcessPriority(self):
        """
        Sets high priority for this pcas process
        """
        try:
            sys.getwindowsversion()
        except:
            isWindows = False
        else:
            isWindows = True
        p = psutil.Process(os.getpid())
        if isWindows:
            try:
                p.set_nice(psutil.HIGH_PRIORITY_CLASS)
            except:
                print str(datetime.datetime.now())[:-3], "In function: " + inspect.stack()[0][3] + \
                      ", Failed setting high priority for this process"
        else:
            try:
                os.nice(-10)
            except IOError:
                print str(datetime.datetime.now())[:-3], "In function: " + inspect.stack()[0][3] + \
                      ", Could not set high priority"

    def iocStats(self):
        """
        Sets the iocAdmin related records
        """
        self.start_time = datetime.datetime.now()
        self.setParam('ENGINEER', 'Alireza Panna')
        self.setParam('LOCATION', 'B1D521D DT-ASUSWIN102')
        self.setParam('RECORD_CNT', len(pvdb.keys()))
        self.setParam('APP_DIR1', str(os.getcwd()))
        self.setParam('UPTIME', str(self.start_time))
        self.setParam('PARENT_ID', os.getpid())
        self.setParam('HEARTBEAT', 0)

    def read(self, reason):
        """
        pcaspy native read method
        """
        format_time = ""
        global XRAY_IOC
        if reason == 'UPTIME':
            format_time = datetime.datetime.now() - self.start_time
            value  = str(format_time).split(".")[0] 
        elif reason == 'TOD':
            value = str(datetime.datetime.now().strftime("%m/%d/%Y %H:%M:%S"))
        elif reason == 'HEARTBEAT':
            value = self.getParam('HEARTBEAT') + 1
            self.setParam('HEARTBEAT', value)
        elif reason == 'XSYNC_RBV':
            value = self.getParam('XSYNC')
            if value == 2:
                XRAY_IOC = EXPERIMENT + 'OXFORD:xray:'
            elif value == 1:
                XRAY_IOC = EXPERIMENT + 'SRI:xray:'
            elif value == 0:
                XRAY_IOC = ''
        else: 
            value = self.getParam(reason)
        return value

    def write(self, reason, value):
        """
        pcaspy native write method
        """
        if reason == "XSYNC":
           self.setParam(reason, value)
        elif reason == 'TRIGGER' and value == 1:
            if self.daq_init_flag == 1:
                NIKON_IMAGE_MODE.put(2)   # set image mode to continous
                NIKON_TRIGGER_MODE.put(1) # set to hard mode (this allows sending trigger through the daq)
                NIKON_ACQUIRE.put(1)
            elif self.daq_init_flag == 0:
                NIKON_IMAGE_MODE.put(2)   # set image mode to continous 
                NIKON_TRIGGER_MODE.put(2) # set to soft mode 
                NIKON_ACQUIRE.put(1)
            self.sync_thread = threading.Thread(target = self.sync)
            self.sync_thread.daemon = True
            self.sync_thread.start()
        elif reason == 'SEND_TRIGGER':
            if self.daq_init_flag == 1:
                self.fid = threading.Thread(target = self.sendTrigger, args = (self.qi2TriggerTask, value))
                self.fid.daemon = True
                self.fid.start()
        elif reason == "DOC":
            self.setParam(reason, value)
        elif reason == "REPORT":
            if value == 1:
                print str(datetime.datetime.now())[:-3], 'POST-SCAN REPORT WILL BE SAVED'
            else:
                print str(datetime.datetime.now())[:-3], 'POST-SCAN REPORT WILL NOT BE SAVED'
            self.setParam(reason, value)
        elif reason == "DUTYCYCLE_RBV":
            self.setParam(reason, value)
        elif reason == "IDLE_RBV":
            self.setParam(reason, value)
        elif reason == "ACQUIRE_RBV":
            self.setParam(reason, value)
        elif reason == "TRIGGER_ACQUIRE_DELAY_RBV":
            self.setParam(reason, value)
        elif reason == "ABORT" and value == 1:
            self.aid = threading.Thread(target = self.allAbort)
            self.aid.daemon = True
            self.aid.start()
        elif reason == "SYNC_MODE_RBV":
            if self.daq_init_flag == 1:
                value = 1
            elif self.daq_init_flag == 0:
                value = 0
            self.setParam(reason, value)
        self.setParam(reason, value)
        self.updatePVs() 

    def pollQi2(self):
        """
        Polls Nikon DS-Qi2 input lines to check for exposing and trigger ready lines
        """
        oldval = np.array([0,0], dtype=np.uint8)
        DAQmxReadDigitalLines(self.qi2PollTask, 1, 1, 0,  oldval, 2, None, None, None)
        while True:
            newval = np.array([0,0], dtype=np.uint8)
            DAQmxReadDigitalLines(self.qi2PollTask, 1, 1, 0,  newval, 2, None, None, None) 
            self.write("TRIGGER_READY_RBV", newval[0])
            self.write("EXPOSING_RBV", newval[1])
            self.updatePVs()
            if newval[0] != oldval[0]:
                if newval[0] == 1: # trigger is ready
                    self.idleTime = self.trigger_not_ready_rbv - self.shutter_close
                    self.setParam('IDLE_RBV', self.idleTime)
                    print str(datetime.datetime.now())[:-3], 'DETECTOR END EXPOSURE'
                    self.shutter_close = time.clock()
                    self.shutterTime = self.shutter_close - self.shutter_open
                    self.setParam('ACQUIRE_RBV', self.shutterTime)
                    print str(datetime.datetime.now())[:-3], 'DETECTOR ACQUIRE TIME:', \
                    self.shutterTime,'s'
                if newval[0] == 0: # trigger is not ready
                    self.trigger_not_ready_rbv = time.clock()
            if newval[1] != oldval[1]:
                if newval[1] == 1: # acquisition started/acquiring
                    print str(datetime.datetime.now())[:-3], 'IDLE TIME:', \
                    self.idleTime,'s'
                    print str(datetime.datetime.now())[:-3], 'DETECTOR START EXPOSURE'
                    self.shutter_open = time.clock()
                    print str(datetime.datetime.now())[:-3], 'TRIGGER_READY->ACQUIRING DELAY:', \
                    self.shutter_open - self.trigger_not_ready_rbv,'s'
                    self.setParam('TRIGGER_ACQUIRE_DELAY_RBV', self.shutter_open - self.trigger_not_ready_rbv)
                if newval[1] == 0: # acquisition ended/not acquiring
                    self.notExposing = time.clock()
                    print str(datetime.datetime.now())[:-3], 'NOT ACQUIRING->NEXT TRIGGER_READY DELAY:', \
                    self.notExposing - self.shutter_close,'s'
            oldval = newval
            self.updatePVs()

    def sendTrigger(self, DAQtaskName, value):
        """
        Andrew's compact code to send trigger to Qi2. Used to be setDigiOut
        """
        if type(value) != np.ndarray:
            if value == 0:
                value = LOW
            else:
                value = HIGH
        self.processDAQstatus(DAQmxWriteDigitalLines(DAQtaskName,1,1,10.0,DAQmx_Val_GroupByChannel,value,self.written,None))  

    def processDAQstatus(self, errorcode):
        if errorcode != 0:
            print str(datetime.now())[:-3], "NI-DAQ error! Code:", errorcode

    def startXray(self):
        """
        Returns when the x-ray is on and outputting x-rays at set values
        """
        if self.getParam('XSYNC') == 2: # x-ray sync is set to oxford/nova
            if caget(XRAY_IOC + 'STATUS_RBV') == 5: # make sure x-ray is not in fault mode.
                print str(datetime.datetime.now())[:-3], 'X-ray is in fault mode!'
                return
            else:
                if caget(XRAY_IOC + 'STATUS_RBV') == 0: # xray is warming
                    print str(datetime.datetime.now())[:-3], 'Waiting for warm up to finish!'
                    return
                if caget(XRAY_IOC + 'STATUS_RBV') == 1: # if xray in standby mode i.e not outputting 
                    # turn on x-ray -> this is output mode
                    caput(XRAY_IOC + 'ON', '1')
                    # wait for x-ray to reach set points
                    while (caget(XRAY_IOC + 'FIRING_RBV') != 1): # this record is sampled at 10Hz in the db
                        time.sleep(0.01)
                    print str(datetime.datetime.now())[:-3], 'X-ray is outputting at set points'
                elif caget(XRAY_IOC + 'STATUS_RBV') == 3 or caget(XRAY_IOC + 'STATUS_RBV') == 2: # if xray is pulsing or outputting
                    caput(XRAY_IOC + 'PULSE_MODE', '0')
                    # switch from pulse to output mode...
                    while (caget(XRAY_IOC + 'FIRING_RBV') != 1):
                        time.sleep(0.01)
                    print str(datetime.datetime.now())[:-3], 'X-ray is outputting at set points'
        elif self.getParam('XSYNC') == 1: # x-ray sync is set to sri
            caput(XRAY_IOC + 'ON', '1') 
            print str(datetime.datetime.now())[:-3], 'X-ray is outputting at set points'
        elif self.getParam('XSYNC') == 0:
            print str(datetime.datetime.now())[:-3], 'Acquiring dark image!'
            return

    def stopXray(self):
        """
        Stops the x-ray flux
        """
        if self.getParam('XSYNC') == 0:
            return
        else:
            caput(XRAY_IOC + 'ON', '0')
        print str(datetime.datetime.now())[:-3], 'X-ray is off'

    def afterScan(self):
         self.stopXray()
         self.setParam('TRIGGER', 0)
         self.elapsedScanTime = str(datetime.timedelta(0, time.clock() - self.scanStartTime))
         if (self.getParam('REPORT') == 1):
             # the first two points don't make sense
             try:
                 self.dutyCycleList.pop(0)
                 self.dutyCycleList.pop(1)
             except:
                 print str(datetime.datetime.now())[:-3], 'In function: ' + inspect.stack()[0][3] + \
                       ', self.dutyCycleList is empty'
                 pass
             self.reportStatistics()
             self.scanReport()
         self.shutterTime   = []
         self.dutyCycleList = []
         self.fig.delaxes(self.ax)
         self.ax = self.fig.add_subplot(1, 1, 1)
         self.scanStartTime = 0
         self.total_cpts = 0
        
    def sync(self):
        """
        Synchronizes the x-ray source exposure (master clock) with the detector shutter (slave)
        """
        # sequence for sync trigger
        if self.scanStartTime == 0:
            self.scanStartTime = time.clock()
        self.startXray()
        if self.daq_init_flag == 1:
            self.write('SEND_TRIGGER', 0)
            time.sleep(0.025)
            # check to see if we are ready to expose
            while self.getParam('TRIGGER_READY_RBV') != 1:
                time.sleep(POLL_TIME)
            self.updatePVs()
            self.write('SEND_TRIGGER', 1)
            time.sleep(0.025)
            while self.getParam('EXPOSING_RBV') == 1:
                time.sleep(POLL_TIME)
            self.updatePVs()
        elif self.daq_init_flag == 0:
            print str(datetime.datetime.now())[:-3], 'DETECTOR START EXPOSURE'
            caput(DET_IOC + 'cam1:SoftTrigger', "1")
            time.sleep(NIKON_ACQUIRE_TIME_RBV.get())
            print str(datetime.datetime.now())[:-3], 'DETECTOR END EXPOSURE'
        # scan logic start
        if SCAN_1_BUSY.get() or SCAN_2_BUSY.get() or SCAN_3_BUSY.get() or SCAN_4_BUSY.get():
            self.total_cpts += 1
            caput(SCAN_IOC + 'scan1.WAIT', 0)
            if self.daq_init_flag == 1:
                self.shutterTimeList.append(self.shutterTime)
            elif self.daq_init_flag == 0:
                self.shutterTimeList.append(NIKON_ACQUIRE_TIME_RBV.get())
            if self.total_cpts > 1 and self.daq_init_flag == 1:
                self.dutyCycleList.append(self.shutterTime/(self.idleTime + self.shutterTime))
                self.setParam('DUTYCYCLE_RBV',100 *self.shutterTime/(self.idleTime + self.shutterTime))
            # Grab the image statistics right after qi2 is finished exposing.
            if (self.getParam('REPORT') == 1):
                self.mean_stats.append(NIKON_STATS_MEAN.get(as_string=False))
                self.sigma_stats.append(NIKON_STATS_SIGMA.get(as_string=False))
            # check if all images from scan are done. 
            if SCAN_DIM.get() == 4:
                self.total_npts = SCAN_1_NPTS.get()*SCAN_2_NPTS.get()*SCAN_3_NPTS.get()*SCAN_4_NPTS.get()
            elif SCAN_DIM.get() == 3:
                self.total_npts = SCAN_1_NPTS.get()*SCAN_2_NPTS.get()*SCAN_3_NPTS.get()
            elif SCAN_DIM.get() == 2:
                self.total_npts = SCAN_1_NPTS.get()*SCAN_2_NPTS.get()
            elif SCAN_DIM.get() == 1:
                self.total_npts = SCAN_1_NPTS.get()
            if (self.total_cpts >= self.total_npts):
                self.afterScan()
            self.updatePVs()
        # single exposure
        else:
            self.stopXray()
            self.setParam('TRIGGER', 0)
        self.updatePVs()

    def scanReport(self):
        """
        Writes sscan info to text file
        """
        with open(NIKON_FILEPATH.get(as_string = True) + os.sep + 'scanReport_' + time.strftime("%Y%m%d-%H%M%S") + '.txt', 'w+') as f:
            f.write('# ' + EXPERIMENT.split(':')[0] + ' Scan Report auto-generated on ' + str(datetime.datetime.now()) + '\n')
            f.write('No. Images in Scan: ' + str(self.total_npts) + '\n')
            f.write('Total Scan time: (HH:MM:SS) ' + self.elapsedScanTime + '\n')
            if self.daq_init_flag == 1:
                f.write('Camera Acquire time per image during scan: ' + str(np.nanmean(self.shutterTime)) +  \
                        ' +/- ' +  str(np.nanstd(self.shutterTime)) + ' s/image ' + 'Min:' + \
                        str(np.nanmin(self.shutterTime)) + ' s/image ' + 'Max: ' + str(np.nanmax(self.shutterTime)) + ' s/image\n')
                f.write('Camera Duty Cyle during scan: ' + str(100*np.nanmean(self.dutyCycleList)) +  \
                        ' +/- ' +  str(100*np.nanstd(self.dutyCycleList)) + ' % ' + 'Min:' + \
                        str(100*np.nanmin(self.dutyCycleList)) + ' % ' + 'Max: ' + str(100*np.nanmax(self.dutyCycleList)) + ' %\n')
            if self.getParam('XSYNC') != 0:
                f.write('Accumulated X-ray Exposure time during scan: ' + str(caget(XRAY_IOC + 'MS_RBV', as_string=False)/60000.) + str(' mins\n'))
                f.write('Accumulated X-ray Exposure mAs during scan: ' +  caget(XRAY_IOC + 'MAS_RBV', as_string=True) + str(' mAs\n'))
    
    def reportStatistics(self):
        """
        Writes image statistics info to text file, useful for alignment purposes during scan
        """
        norm_sigma = []
        with open(NIKON_FILEPATH.get(as_string = True) + os.sep + 'imageStatistics_' + time.strftime("%Y%m%d-%H%M%S") + '.txt', 'w+') as f:
            f.write('# VPFI Image Statistics auto-generated on ' + str(datetime.datetime.now()) + '\n')
            f.write('Mean' + '\t' + 'Sigma' + '\t' + 'Normalized Sigma' + '\n')
            try:
                for m, s in zip(self.mean_stats, self.sigma_stats):
                    norm_sigma.append(s/m)
                    f.write(str(m) + '\t' + str(s) + '\t' + str(s/m) + '\n')
            except:
                pass
        self.ax.plot(norm_sigma, 'ro-', label='Measured', linewidth=1.25)
        self.ax.set_ylabel(r'Normalized Stdev (A.U.)', fontsize=22)
        self.ax.set_xlabel(r'Image No.', fontsize=22)
        self.ax.grid(which='major', axis='x', linewidth=0.75, linestyle='--', color='gray')
        self.ax.grid(which='major', axis='y', linewidth=0.75, linestyle='--', color='gray')
        self.ax.legend(loc='upper left')
        self.fig.savefig(NIKON_FILEPATH.get(as_string = True) + os.sep + 'imageStatistics_' + \
                         time.strftime("%Y%m%d-%H%M%S") + ".jpg", bbox_inches='tight', dpi=80)
        plt.close(self.fig)
        f.close()
        self.mean_stats = []
        self.sigma_stats = []
        norm_sigma = []
        print str(datetime.datetime.now())[:-3], 'Report log Successful.'
    
    def saveParams(self):
        """
        Daemon thread that saves a text file with the same name as 
        the image name. The file contains motor position readback values.
        """
        filePath = caget(DET_IOC + 'TIFF1:FilePath_RBV', as_string=True)
        fullFileName = caget(DET_IOC + 'TIFF1:FullFileName_RBV', as_string=True)
        if fullFileName == '':
            return 
        else:
            self.val = []
            for pvs in MOTOR_IOC_LIST:
                try:
                    self.val.append(caget(pvs + '.RBV'))
                except:
                    pass
            fileName = (fullFileName.split('\\')[-1]).split('.')[0]        
            f = open(filePath + fileName + '.txt', 'w')
            # Save relevant PVs in PV_LIST
            for i, j in zip(MOTOR_IOC_LIST, self.val):
                f.write(i + " - " + str(j) + '\n')
            f.close()
            print str(datetime.datetime.now())[:-3], 'Document successful.'
            self.save_params_thread = None

    def checkDoc(self, **kw):
        """
        Callback function for NikonKs FullFileName_RBV PV. When image is saved,
        this function will start a thread to save parameter file if DOC is ON.
        """
        if self.getParam('DOC') == 1:
            self.save_params_thread = threading.Thread(target = self.saveParams, args=())
            self.save_params_thread.daemon = True
            self.save_params_thread.start()
        else:
            return

    def scanAbortCallback(self, **kw):
        """
        Start the scan abort thread
        """
        self.scan_abort_thread = threading.Thread(target = self.scanAbortSeq, args=())
        self.scan_abort_thread.daemon = True
        self.scan_abort_thread.start()
    
    def scanAbortSeq(self):
        """
        Start the scan abort sequence
        """
        print str(datetime.datetime.now())[:-3], 'START SCAN ABORT SEQUENCE!'
        self.stopXray()
        self.setParam('TRIGGER', 0)
        if self.getParam('REPORT') == 1:
             self.mean_stats         = []
             self.sigma_stats        = []
             self.shutterTime        = []
             self.dutyCycleList      = []
             self.fig.delaxes(self.ax)
             self.ax = self.fig.add_subplot(1, 1, 1)
        print str(datetime.datetime.now())[:-3], 'END SCAN ABORT SEQUENCE!'
        self.scan_abort_thread = None
    
    def allAbort(self):
        """
        Master all abort sequence. Similar to scan abort callback function, but in addition
        stops any moving motors and disables camera acquisition
        """
        print str(datetime.datetime.now())[:-3], 'START MASTER ABORT SEQUENCE!'
        # stop x-ray, motor, disable camera acquisition and clear report lists.
        SCAN_ABORT.put(1)
        for pvs in MOTOR_IOC_LIST:
            if caget(pvs + '.DMOV') == 0:
                caput(pvs + '.STOP', 1)
        NIKON_ACQUIRE.put(0)
        print str(datetime.datetime.now())[:-3], 'END MASTER ABORT SEQUENCE!'
        self.aid = None

if __name__ == '__main__':
    server = SimpleServer()
    server.createPV(prefix, pvdb)
    # create the driver object
    driver = myDriver()
    # Process CA transactions
    while True:
        try:
            server.process(0.0001)
        except KeyboardInterrupt:
            try:
                sys.exit(0)
            except SystemExit:
                os._exit(0)
