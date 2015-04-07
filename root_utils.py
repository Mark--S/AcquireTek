#!/usr/bin/env python
#
# root_utils.py
#
# Utility functions to convert waveforms to the high energy physics root format and save the data.
# -- Further functions added to generate historgrams containing results of standard measurements to 
# saved waveforms - much functionality inhereted from calc_utils.py (Ed).
#
# Author P G Jones - 04/06/2013 <p.g.jones@qmul.ac.uk> : First revision
#        Ed Leming - 07/03/2015 <e.leming09@googlemail.com> : Second revision
#################################################################################################### 
import ROOT
import calc_utils as calc
import numpy as np

def waveform_to_hist(timeform, waveform, data_units, title="hist"):
    """ Pass a tuple of dataforms and data units.
    Loaded values are in divs, must use scalings to convert to correct units if desired."""
    histogram = ROOT.TH1D("data", title, len(timeform), timeform[0], timeform[-1])
    histogram.SetDirectory(0)
    for index, data in enumerate(waveform):
        histogram.SetBinContent(index + 1, data)
    histogram.GetXaxis().SetTitle(data_units[0])
    histogram.GetYaxis().SetTitle(data_units[1])
    return histogram

def plot_area(x, y, name):
    """Calc area of pulses"""
    area, areaErr = calc.calcArea(x,y)
    bins = np.arange((area-8*areaErr)*1e9, (area+8*areaErr)*1e9, (areaErr/5)*1e9)
    hist = ROOT.TH1D("%s" % name,"%s" % name, len(bins), bins[0], bins[-1])
    hist.SetTitle("Pulse integral")
    hist.GetXaxis().SetTitle("Integrated area (V.ns)")
    for i in range(len(y[:,0])-1):
        hist.Fill(np.trapz(y[i,:],x)*1e9)
    return hist, area, areaErr

def plot_rise(x, y, name):
    """Calc and plot rise time of pulses"""
    rise, riseErr = calc.calcRise(x,y)
    bins = np.arange((rise-8*riseErr)*1e9, (rise+8*riseErr)*1e9, (riseErr/5.)*1e9)
    hist = ROOT.TH1D("%s" % name,"%s" % name, len(bins), bins[0], bins[-1])
    hist.SetTitle("Rise time")
    hist.GetXaxis().SetTitle("Rise time (ns)")
    f = calc.positive_check(y)
    if f == True:
        for i in range(len(y[:,0])-1):
            m = max(y[i,:])
            lo_thresh = m*0.1
            hi_thresh = m*0.9
            low = calc.interpolate_threshold(x, y[i,:], lo_thresh)
            high = calc.interpolate_threshold(x, y[i,:], hi_thresh)
            hist.Fill((high - low)*1e9)
    else:
        for i in range(len(y[:,0])-1):
            m = min(y[i,:])
            lo_thresh = m*0.1
            hi_thresh = m*0.9
            low = calc.interpolate_threshold(x, y[i,:], lo_thresh, rise=False)
            high = calc.interpolate_threshold(x, y[i,:], hi_thresh, rise=False)
            hist.Fill((high - low)*1e9)
    return hist, rise, riseErr

def plot_fall(x, y, name):
    """Calc and plot fall time of pulses"""
    fall, fallErr = calc.calcFall(x,y)
    bins = np.arange((fall-8*fallErr)*1e9, (fall+8*fallErr)*1e9, (fallErr/5.)*1e9)
    hist = ROOT.TH1D("%s" % name,"%s" % name, len(bins), bins[0], bins[-1])
    hist.SetTitle("Fall time")
    hist.GetXaxis().SetTitle("Fall time (ns)")
    f = calc.positive_check(y)
    if f == True:
        for i in range(len(y[:,0])-1):
            m = max(y[i,:])
            m_index = np.where(y[i,:] == m)[0][0]
            lo_thresh = m*0.1
            hi_thresh = m*0.9
            low = calc.interpolate_threshold(x, y[i,m_index:], lo_thresh, rise=False)
            high = calc.interpolate_threshold(x, y[i,m_index:], hi_thresh, rise=False)
            hist.Fill((low - high)*1e9)
    else:
        for i in range(len(y[:,0])-1):
            m = min(y[i,:])
            m_index = np.where(y[i,:] == m)[0][0]
            lo_thresh = m*0.1
            hi_thresh = m*0.9
            low = calc.interpolate_threshold(x, y[i,m_index:], lo_thresh)
            high = calc.interpolate_threshold(x, y[i,m_index:], hi_thresh)
            hist.Fill((low - high)*1e9)
    return hist, fall, fallErr

def plot_peak(x, y, name):
    """Plot pulse heights for array of pulses"""
    peak, peakErr = calc.calcPeak(x,y)
    bins = np.arange((peak-8*peakErr), (peak+8*peakErr), (peakErr/5.))
    hist = ROOT.TH1D("%s" % name,"%s" % name, len(bins), bins[0], bins[-1])
    hist.SetTitle("Pulse hieght")
    hist.GetXaxis().SetTitle("Pulse height (V)")
    f = calc.positive_check(y)
    if f == True:
        for i in range(len(y[:,0])-1):
            hist.Fill(max(y[i,:]))
    else:
        for i in range(len(y[:,0])-1):
            hist.Fill(min(y[i,:]))
    return hist, peak, peakErr

def plot_jitter(x1, y1, x2, y2, name):
    """Calc and plot jitter of pulse pairs"""
    sep, jitter, jittErr = calc.calcJitter(x1, y1, x2, y2)
    bins = np.arange((sep-8*jitter)*1e9, (sep+8*jitter)*1e9, (jitter/5.)*1e9)
    hist = ROOT.TH1D("%s" % name,"%s" % name, len(bins), bins[0], bins[-1])
    hist.SetTitle("Jitter between signal and trigger out")
    hist.GetXaxis().SetTitle("Pulse separation (ns)")
    p1 = calc.positive_check(y1)
    p2 = calc.positive_check(y2)
    for i in range(len(y1[:,0])-1):
        m1 = calc.calcSinglePeak(p1, y1[i,:])
        m2 = calc.calcSinglePeak(p2, y2[i,:])
        time_1 = calc.interpolate_threshold(x1, y1[i,:], 0.1*m1, rise=p1)
        time_2 = calc.interpolate_threshold(x2, y2[i,:], 0.1*m2, rise=p2)
        hist.Fill((time_1 - time_2)*1e9)
    return hist, jitter, jittErr

def fit_gauss(hist):
    """Fit generic gaussian to histogram"""
    f = ROOT.TF1("f1","gaus")
    f.SetLineColor(1)
    p = hist.Fit(f, "S")

    # Write to canvas
    #stats = c1.GetPrimitive("stats")
    #stats.SetTextColor(1)
    #c1.Modified(); c1.Update()

    return f.GetParameters(), f.GetParErrors()

def print_hist(hist, savename, c):
    """Function to print histogram to png"""
    c.Clear()
    hist.Draw("")
    c.Update()
    stats = c.GetPrimitive("stats")
    stats.SetTextSize(0.04)
    c.Modified(); c.Update()
    c.Print("%s" % savename, "pdf")
