###############################################
# Functions to read in pickled scope traces
# and perform standard measurements.
# 
# Will calculate: Integrated area, rise time,
# fall time, pulse width, peak voltage and
# the jitter on two signals. 
# 
# Author: Ed Leming
# Date: 17/03/2015
################################################
import pickle
import utils
import time
import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import array

def readPickleChannel(file, channel_no, correct_offset=True):
    """Read data set as stored in pickle file"""
    # Make sure file path is correct format
    if file[-4:] == ".pkl":
        file = file[0:-4]
    ### READ Pickle File ###
    wave = utils.PickleFile(file, 4)
    wave.load()
    xRaw = wave.get_meta_data("timeform_%i" % channel_no)
    yRaw = wave.get_data(channel_no)
    # Correct for trigger offset in timestamps
    x = xRaw - xRaw[0]
    # Count how many pulses saved the file
    count = 0
    for i in yRaw:
        count = count + 1
    ### Make 2D array of pulse y values ###
    y = np.zeros( (count, len(xRaw)) )
    for i, ent in enumerate(yRaw):
        if correct_offset == True:
            y[i, :] = ent  - np.mean(ent[0:20])
        else:
            y[i, :] = ent  
    return x,y

def positive_check(y):
    if np.mean(y[1,:]) > 0:
        return True
    else:
        return False

def calcArea(x,y):
    """Calc area of pulses"""
    trapz = np.zeros( len(y[:,0]) )
    for i in range(len(y[:,0])):
        trapz[i] = np.trapz(y[i,:],x)
    return np.mean(trapz), np.std(trapz)

def interpolate_threshold(x, y, thresh, rise=True, start=0):
    """Calculate the threshold crossing using a linear interpolation"""
    if rise == True:
        index_high = np.where( y > thresh )[0][start]
    else:
        index_high = np.where( y < thresh )[0][start]
    index_low = index_high - 1
    dydx = (y[index_high] - y[index_low])/(x[1]-x[0])
    c = -dydx * x[index_low]
    return ((thresh - c) / dydx)

def calcRise(x,y):
    """Calc rise time of pulses"""
    rise = np.zeros( len(y[:,0]) )
    f = positive_check(y)
    if f == True:
        for i in range(len(y[:,0])):
            m = max(y[i,:])
            lo_thresh = m*0.1
            hi_thresh = m*0.9
            low = interpolate_threshold(x, y[i,:], lo_thresh)
            high = interpolate_threshold(x, y[i,:], hi_thresh)
            rise[i] = high - low
        return np.mean(rise), np.std(rise)
    else: 
        for i in range(len(y[:,0])):
            m = min(y[i,:])
            lo_thresh = m*0.1
            hi_thresh = m*0.9
            low = interpolate_threshold(x, y[i,:], lo_thresh, rise=False)
            high = interpolate_threshold(x, y[i,:], hi_thresh, rise=False)
            rise[i] = high - low
        return np.mean(rise), np.std(rise)

def calcFall(x,y):
    """Calc fall time of pulses"""
    fall = np.zeros( len(y[:,0]) )
    f = positive_check(y)
    if f == True:
        for i in range(len(y[:,0])):
            m = max(y[i,:])
            m_index = np.where(y[i,:] == m)[0][0]
            lo_thresh = m*0.1
            hi_thresh = m*0.9
            low = interpolate_threshold(x, y[i,m_index:], lo_thresh, rise=False)
            high = interpolate_threshold(x, y[i,m_index:], hi_thresh, rise=False)
            fall[i] = low - high
        return np.mean(fall), np.std(fall)
    else:
        for i in range(len(y[:,0])):
            m = min(y[i,:])
            m_index = np.where(y[i,:] == m)[0][0]
            lo_thresh = m*0.1
            hi_thresh = m*0.9
            low = interpolate_threshold(x, y[i,m_index:], lo_thresh)
            high = interpolate_threshold(x, y[i,m_index:], hi_thresh)
            fall[i] = low - high
        return np.mean(fall), np.std(fall)
        
def calcWidth(x,y):
    """Calc width of pulses"""
    width = np.zeros( len(y[:,0]) )
    f = positive_check(y)
    if f == True:
        for i in range(len(y[:,0])):
            m = max(y[i,:])
            thresh = m*0.5
            index_1 = interpolate_threshold(x, y[i,:], thresh)
            index_2 = interpolate_threshold(x, y[i,:], thresh, start=-1)
            width[i] = index_2 - index_1
        return np.mean(width), np.std(width)
    else:
        for i in range(len(y[:,0])):
            m = min(y[i,:])
            thresh = m*0.5
            index_1 = interpolate_threshold(x, y[i,:], thresh, rise=False)
            index_2 = interpolate_threshold(x, y[i,:], thresh, rise=False, start=-1)
            width[i] = index_2 - index_1
        return np.mean(width), np.std(width)

def calcPeak(x,y):
    """Calc min amplitude of pulses"""
    peak = np.zeros( len(y[:,0]) )
    f = positive_check(y)
    if f == True:
        for i in range(len(y[:,0])):
            peak[i] = max(y[i,:])
        return np.mean(peak), np.std(peak)
    else:
        for i in range(len(y[:,0])):
            peak[i] = min(y[i,:])
        return np.mean(peak), np.std(peak)

def calcSinglePeak(pos_check, y_arr):
    """Calculate peak values for single trace inputs can be positive or negative."""
    if pos_check == True:
        m = max(y_arr)
    else:
        m = min(y_arr)
    return m

def calcJitter(x1, y1, x2, y2):
    """Calc jitter between trig and signal using CFD"""
    p1 = positive_check(y1)
    p2 = positive_check(y2)
    times = np.zeros(len(y1[:,0]))
    for i in range(len(y1[:,0])):
        m1 = calcSinglePeak(p1, y1[i,:])
        m2 = calcSinglePeak(p2, y2[i,:])
        time_1 = interpolate_threshold(x1, y1[i,:], 0.1*m1, rise=p1)
        time_2 = interpolate_threshold(x2, y2[i,:], 0.1*m2, rise=p2)
        times[i] = time_1 - time_2
    return np.mean(times), np.std(times), np.std(times)/np.sqrt(2*len(y1[:,0]))

def printParams(x,y, name):
    """Calculate standard parameters and print to screen"""
    area, areaStd = calcArea(x,y)
    rise, riseStd = calcRise(x,y)
    fall, fallStd = calcFall(x,y)
    width, widthStd = calcWidth(x,y)
    mini, miniStd = calcPeak(x,y)

    print "\n%s:" % name
    print "--------"
    print "Area \t\t= %1.2e +/- %1.2e Vs" % (area, areaStd)
    print "Fall time \t= %1.2f +/- %1.2f ns" % (fall*1e9, fallStd*1e9)
    print "Rise time \t= %1.2f +/- %1.2f ns" % (rise*1e9, riseStd*1e9)
    print "Width \t\t= %1.2f +/- %1.2f ns" % (width*1e9, widthStd*1e9)
    print "Peak \t\t= %1.2f +/- %1.2f V" % (mini, miniStd)

def plot_eg_pulses(x,y,n,title=None,fname=None,show=False):
    """Plot example pulses""" 
    plt.figure()
    for i in range(n):
        plt.plot(x*1e9,y[i,:])
    if title == None:
        plt.title( "Example pulses")
    else:
        plt.title(title)
    plt.xlabel("Time (ns)")
    plt.ylabel("Amplitude (V)")
    if fname is not None:
        plt.savefig(fname)
    if show == True:
        plt.show()
    plt.clf()
