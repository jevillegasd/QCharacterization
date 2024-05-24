#Additional functions for basic qubit parametrization with vna and local oscilators
#Original Fucntions developed by: Lucas July 2023

import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt, row_stack
import scipy.interpolate as spy
from scipy import odr
from scipy.optimize import curve_fit
from E5080B_driver import E5080B_driver

#Spectrometry one-liner
def meas_spectrum(instr: E5080B_driver, f_c :float= 7.33e9, f_span:float = 80, power = -60, sleep_t = 3, npoints:int = 1001, data_format = 'MA', avg_count = 1000, meas = 'S21', overule_power = False, *arg, **kwargs): 
    from time import sleep
    instr.set_data_format(data_format)
    instr.set_center_frequency(f_c)
    instr.set_span_frequency(f_span)
    instr.set_power(power,overule_power)
    instr.set_average_count(avg_count)
    instr.set_sweep_npoints(npoints )

    instr.reset_average()
    # instr.send(":CALC1:MEAS1:CORR:EDEL:TIME %2.4e"%delay)

    

    state= "+3\n"
    while state  == "+3\n":
        sleep(0.02)
        try:
            state = instr.query("STAT:QUES:INT:MEAS1:COND?")
        except Exception:
            print(Exception)
        
    sleep(0.1)
    f = instr.get_freq_array()
    instr.start_rf()
    sleep(sleep_t)
    mag, phi = instr.get_data(meas)
    return f, mag, phi

def lorentzian( x, x0, a, gam ):
    return a/(2*np.pi) * gam / ( (gam/2)**2 + ( x - x0 )**2)

def Q2(f,mag, plot_f = False):
    '''
    Loaded quality factor, fit to a conventional Lorentzian using the scipy signal fitting tools
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.signal import find_peaks, peak_widths

    thresh_top    = np.median(mag) + 2 * np.std(mag)
    thresh_bottom = np.median(mag) - 2 * np.std(mag)

    peaks_idx, _   = find_peaks(mag,  height =  thresh_top)
    valley_idx, _ = find_peaks(-mag, height = -thresh_bottom)

    f_c = f[peaks_idx]
    df = (f[1]-f[0])
    peaks_w = peak_widths(mag, peaks_idx, rel_height=0.2)
    fwhm = peaks_w[0]*df
    Q = f_c[0]/fwhm[0]
    if plot_f:
        plt.plot(f,mag)
    return f_c[0], Q, fwhm[0]


def Q(f, mag, plot_f = False, format='MA'):
    '''
    Loaded quality factor, fit to a conventional Lorentzian
    '''
    from scipy.optimize import curve_fit
    from scipy.signal import find_peaks

    num_samp = round(len(f)/10)

    if format == 'DB':
        mag = 10**(mag/20)

    mag_min = mag.min()
    mag_lin = mag -mag.min()

    min_w = int(len(f)/100)

    peaks, _ = find_peaks(mag_lin, width=min_w, prominence=mag_lin.max()/2)
    flip = False
    if len(peaks) == 0:
        zero = mag_lin.max()
        mag_lin =zero -mag_lin
        peaks, _ = find_peaks(mag_lin, width=min_w, prominence=mag_lin.max()/2)
        flip = True
        idx_fmax = peaks[0]
    else:
        idx_fmax = peaks[0]

    max_v = mag_lin.max() # Normalization constant
    yData = mag_lin / max_v #Normalize
    xData = range(len(f))

    param, _ = curve_fit(lorentzian, 
                         xData[idx_fmax-int(num_samp/2):idx_fmax+int(num_samp/2)], 
                         yData[idx_fmax-int(num_samp/2):idx_fmax+int(num_samp/2)])
    
    mag_fit = [ lorentzian(x, *param)*max_v for x in xData ]

    if flip:
        mag_fit = zero-mag_fit
    mag_fit = mag_fit+mag_min

    if format == 'DB':
        mag_fit = 20*np.log10(mag_fit)


    f_c = f[idx_fmax]
    df = f[1]-f[0]
    FWHM = abs(df*param[-2]) # Gamma = param[-1]
    Q = f_c/FWHM
    if plot_f:
        plt.plot(f,mag)
        plt.plot(f,mag_fit)
        pass
    return f_c, Q, FWHM

def find_max(freqs,mag):
    return freqs[np.where(mag==np.max(mag))][0]

def find_min(freqs,mag):
    return freqs[np.where(mag==np.min(mag))][0]

def find_flux_resonators_frequencies_max(currents,freqs,mags):
    resonator_frequencies = np.ndarray(len(currents))
    for i in range(len(currents)):
        y = find_max(freqs,mags[i])
        resonator_frequencies[i] = y
    return resonator_frequencies

def find_flux_resonators_frequencies_min(currents,freqs,mags):
    resonator_frequencies = np.ndarray(len(currents))
    for i in range(len(currents)):
        y = find_min(freqs,mags[i])
        resonator_frequencies[i] = y
    return resonator_frequencies

def find_rough_estimate_Q(freqs,mags,format = "DB",width = 20):
    '''MA: linear DB: Log'''

    #print(format)
    if format != "MA" and format != "DB":
        raise ValueError("Wrong format. It should be MA or DB")

    if format == "MA":
        mags = 20*np.log10(mags);

    half = -3;

    max = np.max(mags);

    mags = mags - max -half;

    f0 = find_max(freqs,mags)
    idx_f0 = np.where(mags == 3)[0][0]


    a_min = 3
    a_idx  = 0
    for i in range(int(idx_f0-width/2),idx_f0):
        aux = np.abs(mags[i])
        if aux < a_min:
            a_min = aux
            a_idx  = i

    
    b_min = 3
    b_idx  = 0
    for i in range(idx_f0,int(idx_f0+width/2)):
        aux = np.abs(mags[i])
        if aux < b_min:
            b_min = aux
            b_idx  = i

    
    delta = freqs[b_idx]-freqs[a_idx]

    #print(b_idx)
    #print(a_idx)
    #print(idx_f0)

    return f0/delta
    
    

def fitCPW(f,Q,f0,Qc,floor):
    """ Complex fit using Lorentzian curve and error term.

    """
    
    S21 = (1-(2*Q/Qc)/(1+2j*Q*(f-f0)/f0))*floor

    return S21

def fitCPW_real(f,Q,f0,Qc,floor):
    """ Complex fit using Lorentzian curve and error term.
    """
    
    return fitCPW(f,Q,f0,Qc,floor).real

µs = 1e-6
ns = 1e-9

# Wrap phase______________
def unwrap_d(Phase):
    """
    Unwrap Phase in Degree"""
    return np.unwrap(Phase * np.pi / 180)


def wrap(Phase):
    """ Wrap back the phase 
    """
    return (Phase + np.pi) % (2 * np.pi) - np.pi


def wrap1(Phase):
    """ Wrap back the phase using arctan2
    """
    return np.arctan2(np.sin(Phase), np.cos(Phase))

    


def dBmtoV(Amp):
    """Convert dBm to Volts
    Input Amplitude in dBm
    """
    c = 10 * np.log10(20)
    return np.power(10, (Amp - c) / 20)


def VtodBm(V):
    """Convert Volts to dBm
    Input Voltage in Volts
    """
    c = 10 * np.log10(20)
    return 20 * np.log10(V) + c

# Convert dBm to W

def WtodB(W):
    return 10 * np.log10(W)

def dBtoW(Amp):
    """Convert dB to Linear Scale 
    Imput Amp in dB
    """
    return np.power(10, (Amp) / 10)

def dBmtoW(Amp):
    """Convert dBm to Watts
    Imput Amp in dBm
    """
    return np.power(10, (Amp - 30) / 10)

def toRadians(degree):
    return degree*np.pi/180

def toDegree(radians):
    return radians*180/np.pi

def toComplex(mag,phase):
    return mag*np.exp(1j*phase)

def plot_sparam(freqs, Sparam):
    plt.plot(freqs,Sparam.real)
    plt.plot(freqs,Sparam.imag)

def plot_sparam_abs(freqs, Sparam):
    plt.plot(freqs,np.absolute(Sparam))
    
def plot_sparam_angle(freqs, Sparam):
    plt.plot(freqs,np.angle(Sparam))
    
def plot_sparam_polar(Sparam):
    fig = plt.figure(figsize=(10,10))
    ax = fig.gca()
    ax.scatter(Sparam.real,Sparam.imag)
    
    return ax

def remove_cable_delay(freqs, Sparam, delay):
    return Sparam*np.exp(2*np.pi*1j*delay*freqs)

def find_center(Sparam,x0,y0):
    """ Finds the center of a circle. Give the lists of X and Y, also a estimate of the center. Use odr method. """
    X = Sparam.real
    Y = Sparam.imag
    
    # any point on the plot, which must be a circle
    x1 = X[0]
    y1 = Y[0]

    def calc_R(xc, yc):
        """ calculate the distance of each 2D points from the center (xc, yc) """
        return sqrt((x1-xc)**2 + (y1-yc)**2)

    x_m = x0
    y_m = y0

    # initial guess for parameters
    R_m = calc_R(x_m, y_m).mean()
    beta0 = [ x_m, y_m, R_m]


    # method_3 = "odr"

    def f_3(beta, x):
        """ implicit definition of the circle """
        return (x[0]-beta[0])**2 + (x[1]-beta[1])**2 -beta[2]**2

    lsc_data  = odr.Data(row_stack([X, Y]),y=1)#,we=error,wd=error)
    lsc_model = odr.Model(f_3, implicit=True)
    lsc_odr   = odr.ODR(lsc_data, lsc_model, beta0)
    lsc_out   = lsc_odr.run()

    xc_3, yc_3, R_3 = lsc_out.beta
    Ri_3 = calc_R(xc_3, yc_3)
    #    residu_3 = sum((Ri_3 - R_3)**2)
    x0 = lsc_out.beta[0]
    y0 = lsc_out.beta[1]
    Zc = complex(x0,y0)
    
    return Zc,R_3

def rotate_Sparam(Zc,Sparam):
    return (Zc - Sparam)*np.exp(-1j*np.angle(Zc))+np.abs(Zc)

def fit_peak(X,Y,n):
    from scipy.optimize import curve_fit
    def func(x, a, b, c):
        return a*(x-b)**2+c

    peak = np.where(Y == np.min(Y))[0][0]
    begin = peak-n
    end = peak +n
    n_interpolate = 41
    xnew = np.linspace(X[begin], X[end-1], num=n_interpolate, endpoint=True)
    popt, pcov = curve_fit(func, X[begin:end], Y[begin:end],p0=[1,X[peak],Y[peak]])
    ynew= func(xnew,*popt)
    
    return xnew,ynew

def S2PParser(filename,starting_line = 5):
    '''S2P files parser. Useful to parse the data returned by our network analyser'''
    file = open(filename,"r")
    lines = file.readlines()
    file.close()
    info = ''.join(lines[:5])
    
    Freqs = []
    S11s_dB = []
    S11s_angle = []
    S21s_dB = []
    S21s_angle = []
    S12s_dB = []
    S12s_angle = []
    S22s_dB = []
    S22s_angle = []


    for line in lines[starting_line:]:
        data = line.split('\n')[0].split('\t')

        Freqs.append(float(data[0]))

        S11s_dB.append(float(data[1]))
        S11s_angle.append(float(data[2]))

        S21s_dB.append(float(data[3]))
        S21s_angle.append(float(data[4]))

        S12s_dB.append(float(data[5]))
        S12s_angle.append(float(data[6]))

        S22s_dB.append(float(data[7]))
        S22s_angle.append(float(data[8]))
    
    return (info,Freqs,S11s_dB,S11s_angle,S21s_dB,S21s_angle,S12s_dB,S12s_angle,S22s_dB,S22s_angle)

def bigplot(X,Y,xlabel,ylabel,title):
    '''Helper function to get east bigplots'''
    fig,ax = plt.subplots(figsize=(14,7))
    ax.plot(X,Y,linewidth=3)
    ax.set_xlabel(xlabel,fontsize=20)
    ax.set_ylabel(ylabel,fontsize=20)

    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)

    ax.set_title(title,fontsize=20)
    
    return fig,ax


def YfromS(S11,S12,S21,S22,gamma):
    det1_S = (1+S22)*(1+S11) - (1+S21)*(1+S12)
    
    Y11 = sqrt(gamma)*( (1-S11)*(1+S22)+(1-S12)*(-1-S21) )/det1_S*sqrt(gamma)
    Y12 = sqrt(gamma)*( (1-S11)*(-1-S12)+(1-S12)*(1+S11) )/det1_S*sqrt(gamma)
    Y21 = sqrt(gamma)*( (1-S21)*(1+S22)+(1-S22)*(-1-S21) )/det1_S*sqrt(gamma)
    Y22 = sqrt(gamma)*( (1-S21)*(-1-S12)+(1-S22)*(1+S11) )/det1_S*sqrt(gamma)
    
    return Y11, Y12, Y21, Y22
    
def get_plot_str(params, filename=""):
    if filename == "":
        filename = params['test_type']+"_"+params['qubit_name']+"_"+params['test_date']
    header = "import numpy as np \n"
    for key in params.keys():
        header += "# " + key +': \t'+ str(params[key]) + "\n"
    
    if params['test_type'] == 'twotone_fpsweep' or params['test_type'] == 'sweep_resonatorquality':
        plot_str = header+\
        "data = np.load('"+filename+".npz')\n\
        fig = plt.figure(figsize=(10,7))\n\
        ax = fig.gca()\n\
        plt.pcolor(data['axis2'],data['freqs'],20*np.log10(np.abs(data['mags'])))\n\
        ax.tick_params(labelsize=20)\n\
        ax.set_xlabel('Qubit Frequency (MHz)',fontsize=20)\n\
        ax.set_clabel('S21 (dB)',fontsize=20)\n\
        ax.set_title('"+params['qubit_name']+" "+params['test_type']+"',fontsize=16)\n\
        plt.show()"
    else:
        plot_str = header+\
        "data = np.load('"+filename+".npz')\n\
        fig = plt.figure(figsize=(10,7))\n\
        ax = fig.gca()\n\
        plt.plot(data['freqs'],20*np.log10(np.abs(data['mags'])))\n\
        ax.tick_params(labelsize=20)\n\
        ax.set_xlabel('Qubit Frequency (MHz)',fontsize=20)\n\
        ax.set_ylabel('S21 (dB)',fontsize=20)\n\
        ax.set_title('"+params['qubit_name']+" "+params['test_type']+"',fontsize=16)\n\
        plt.show()"
    #print(plot_str)
    return plot_str

def punch_out(instr, f_c=7.33e9, hp = -10, lp = -65, f_span = 30, meas = 'S21',overule_power=True):
    f, m,_ =  meas_spectrum(instr, f_c = f_c, f_span = f_span, power =hp, sleep_t = 5, npoints = 1000, data_format = 'MA', meas = meas, overule_power=overule_power)
    f_c0, Qf0,_ = Q(f,m, True)
    f, m , _ =  meas_spectrum(instr, f_c = f_c, f_span = f_span, power =lp, sleep_t = 30, npoints = 1000, data_format = 'MA',meas = meas)
    f_c1, Qf1,_ = Q(f,m, True)

    df = f_c1 - f_c0
    return df