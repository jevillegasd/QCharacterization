"""
This code fits given zdata for transmission resonators. 
David Eslava Sabaté @ Qilimanjaro 22/03/2022 - david.eslava@qilimanjaro.tech

"""

from calendar import c
from cmath import log10
from fileinput import filename
from pyexpat import model
from queue import Empty
from re import A
from tkinter import CENTER
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import lmfit
from scipy.stats import linregress
import matplotlib.gridspec as gridspec
import os
import h5py

#File postprocessing
def read_s2p(file):
    df = pd.read_csv(file, delimiter='\t', skiprows=5,names=['frequency', 'S11_dB', 'S11_ang', 'S21_dB', 'S21_ang', 'S12_dB', 'S12_ang', 'S22_dB', 'S22_ang'])
    return df

def read_txt(file):
    df = pd.read_csv(file, delim_whitespace=True,  header=None, skiprows=5).T
    df.rename(columns={0: "frequency", 1: "s21"}, inplace=True)
    df['frequency'] = df['frequency'].astype(float)
    df['s21'] = df['s21'].astype(complex)
    return df

def read_hdf5(file, fc_inx, pwr_inx ):
    """
    From a hdf5 file creates a dataframe with the dataframe that contains frequency, magnitude in dB and angle in degrees.
    Input: 
        file (string), file path.
        fc_inx (int), index of the frequency center that is wanted to access. 
        pwr_inx (int), index of the power that is wanted to access. 
    Output: 
        df (DataFrame), dataframe that contains "frequency", "S21_dB" and "S21_ang".

    """
    rawdata = h5py.File(file, 'r')
    rawdata.keys()
    xparams = rawdata['Step config']
    xparams.keys()
    stepinfo = xparams['VNA - # of points']
    npoints = int(stepinfo['Step items'][0][2]) #0 is numbers (the only value of the array) and the second one is the number of points. 
    stepinfo['Step items']
    span = xparams['VNA - Span']['Step items'][0][2]
    stepsize = span/(npoints-1)
    def extract(a): #always accessing the second item because is stored as shown for npoints.
        r = []
        for x in a:
            r.append(x[2])
        return np.array(r)
    centerfreqs = extract(xparams[ 'VNA - Center frequency']['Step items'][:]) #extract center frequencies from parameters array
    trace = rawdata['Traces']['VNA - S21']
    trace = np.moveaxis(np.moveaxis(trace,0,-1),0,1) #relabeling if you have a order of frequency, phase, amplitudes change the order. 
    #moving the first one to the last and the second to the first on the other upper array
    trace.shape
    trace.shape #first element is the number of power points, second the number of traces in this case 2 for I and Q and last element the number of points.
    #54 is from 9 central frequencies * 6 different powers.
    data = trace.reshape(9,6,2,-1) #-1 is because is overspecified and is defined by the other 3 parameters 
    #data[fr,P,IQ,Voltages]
    params = data[:,:,:,0]
    S21 = data[:,:,0,:]+1j*data[:,:,1,:] #array of the measured data
    X2 = params[:,0,0] #array of frequencies. This would have been the same as [:,0,0,0]. size 9
    X3 = params[0,:,0] # array of powers. size 6
    # generate frequency array for each trace
    X1 = np.zeros((len(X2),len(X3),npoints))
    for x2 in range(len(X2)): #this to go over the 54 total traces.
        for x3 in range(len(X3)):
            X1[x2,x3] = np.linspace(-span/2,span/2,npoints) + centerfreqs[x2] #constructing the frequency array for each line trace.
    #create dataframe df (so we are compatible with read_s2p and read_txt)
    s21_complex = np.array(S21[fc_inx, pwr_inx].T,dtype=np.complex128)
    df = pd.DataFrame(
    {
        "frequency": X1[fc_inx, pwr_inx],
        "S21_dB": 10*np.log10(np.abs(s21_complex)), #in dB!
        "S21_ang": np.angle(s21_complex,True) #in degrees!
    }
)
    return df

def file2data(file, tau, c, theta):
    """
    From a file gets zdata and freq to do the fitting. File types supported are .txt, .s2p and 
    Input: 
        file (string), file path.
        tau (float), the electrical delay in ns.
        c (float), the gain.
        theta (float), the rotation angle.
    Output: 
        zdata (ndarray(np.complex128)), complex numbers with R and I values of the measurement
        zdata (ndarray(np.complex128)), complex numbers with R and I values of the measurement after phase unwind.
        zdata_r (ndarray(np.complex128)), complex numbers with R and I values of the measurement after phase unwind, escale and rotation.
        freq (ndarray), x-axis frequency values.

    """
    if 'txt' in file:
        df = read_txt(file)
        s21_complex = df.s21
    elif 's2p' in file:
        df = read_s2p(file)
        voltage = 10**(df.S21_dB/20)
        s21_complex = np.array(voltage*np.exp(1j*df.S21_ang*np.pi/180.),dtype=np.complex128) #dtype = # 128-bit complex floating-point number
    #---trying linear magnitude IN VOLTS!
    freq = df.frequency.values
    zdata = np.array(s21_complex*np.exp(1j*freq*tau*1e-9*2*np.pi),dtype=np.complex128) # This is the unwinded data.Add the electrical delay as exp^(j·w·tau)
    zdata_r = c*zdata*np.exp(1j*theta) #rotated zdata
    return s21_complex, zdata, zdata_r, freq

#general fit
def lorentzian_fit(zdata, freq, plot = True, report = True, modelinfo = True): 
    """
    Fits a lorentzian into given zdata magnitude.
    Input: 
        zdata (ndarray(np.complex128)), complex numbers with R and I values of the measurement.
        freq (ndarray), x-axis frequency values.
        plot (Boolean), when true plots the data with the fit
        report (Boolean), when true stores the plot in a .pdf file
        modelinfo (Boolean), when true prints the fit report as well as the initial guess
    Output: 
        f0 (float), fitted center frequency of the resonator.
        BW (float), 3dB bandwidth of the resonator.
        Q (float), Loaded Quality factor of the resonator.
        offset (float), offset 

    """
    frequency = freq
    data = np.abs(zdata)
    voltage_min_i = np.argmin(data)
    voltage_max_i = np.argmax(data)
    voltage_max = np.max(data)
    voltage_min = np.min(data)

    if voltage_min*0.85 <= data[0] <= voltage_min*1.15:
        peak = 'max'
    else:
        peak = 'min'

    #mathematical model of a Lorentzian
    def resonator_peak(frequency,amplitude,center,sigma,offset):
        return (amplitude/np.pi) * (sigma/((frequency-center)**2 + sigma**2) + offset)
    model_Q = lmfit.Model(resonator_peak)

    #Guess parameters
    #to guess center
    if peak == 'max':
        guess_center = frequency[np.argmax(data)] #Argmax = Returns the indices of the maximum values along an axis.
    else:
        guess_center = frequency[np.argmin(data)] #Argmin = Returns the indices of the minimum values along an axis.

    model_Q.set_param_hint('center',value=guess_center,vary=True)

    #to guess the sigma
    guess_sigma = 50e3

    model_Q.set_param_hint('sigma',value=guess_sigma,
                                vary=True)
    #to guess the amplitude 
    #http://openafox.com/science/peak-function-derivations.html
    #Better sigma estimation is a better amplitude estimation
    if peak == 'max':
        guess_amp = voltage_max*guess_sigma*np.pi
    else: 
         guess_amp = -voltage_min*guess_sigma*np.pi       

    model_Q.set_param_hint('amplitude',value=guess_amp,
                                vary=True)
    #to guess the offset
    guess_offset = voltage_min
    model_Q.set_param_hint('offset',value=guess_offset,
                            vary=True)
    #guessed parameters
    guess_parameters = model_Q.make_params()

    #fit the model with the data and guessed parameters
    fit_res = model_Q.fit(data=data,frequency=frequency,params=guess_parameters)

    #get the values for output.
    f0 = fit_res.best_values['center']
    BW = (fit_res.best_values['sigma']*2)
    Q = abs(f0/BW)
    offset = fit_res.best_values['offset']

    if plot == True:
        #plot the fitted curve
        dummy_frequencies = np.linspace(np.amin(frequency),np.amax(frequency),101)
        fit_fine = resonator_peak(dummy_frequencies,**fit_res.best_values)
        fig = plt.figure(figsize=(8, 6))
        ax = plt.subplot()
        ax.plot(frequency,20*np.log(data),'o',label='Data') #plot data in dB
        ax.plot(dummy_frequencies,20*np.log(fit_fine),'r-', label='Lorentzian fit') #plot data fit in dB
        ax.set_ylabel('S21 (dB)')
        ax.set_xlabel('Frequency (GHz)')
        ax.legend()
        print("\nf0: {:.4f} GHz\n"
              "\nQ: {:.4f}\n" 
              "\nBW: {:.4f} MHz\n"
              "\noffset: {:.4f} \n".format(f0/1e9, Q, BW/1e6, offset))
        if report == True:
            fig.savefig('Lorentzian fit.pdf',format='pdf')

    if modelinfo == True:
        fig,ax = plt.subplots(1,1,figsize=(8,3))
        print("\nGuess_parameters:\n")
        print(guess_parameters)
        print("\nFit report:\n")
        print(fit_res.fit_report())
        print("\nFit best values:\n")
        print(fit_res.best_values)
        fit_res.plot_fit(show_init=True)

    return f0, BW, Q, offset

def lorentzian_fit_typeB(zdata, freq, plot = True, report = True, modelinfo = True): 
    """
    Fits a lorentzian into given zdata magnitude.
    Input: 
        zdata (ndarray(np.complex128)), complex numbers with R and I values of the measurement.
        freq (ndarray), x-axis frequency values.
        plot (Boolean), when true plots the data with the fit
        report (Boolean), when true stores the plot in a .pdf file
        modelinfo (Boolean), when true prints the fit report as well as the initial guess
    Output: 
        f0 (float), fitted center frequency of the resonator.
        BW (float), 3dB bandwidth of the resonator.
        Q (float), Loaded Quality factor of the resonator.
        offset (float), offset 

    """
    frequency = freq
    data = np.abs(zdata)
    voltage_min_i = np.argmin(data)
    voltage_max_i = np.argmax(data)
    voltage_max = np.max(data)
    voltage_min = np.min(data)

    peak = 'max'

    #mathematical model of a Lorentzian
    def resonator_peak(frequency,amplitude,center,sigma,offset):
        return (amplitude/np.pi) * (sigma/((frequency-center)**2 + sigma**2) + offset)
    model_Q = lmfit.Model(resonator_peak)

    #Guess parameters
    #to guess center
    if peak == 'max':
        guess_center = frequency[np.argmax(data)] #Argmax = Returns the indices of the maximum values along an axis.
    else:
        guess_center = frequency[np.argmin(data)] #Argmin = Returns the indices of the minimum values along an axis.

    model_Q.set_param_hint('center',value=guess_center,vary=True)

    #to guess the sigma
    guess_sigma = 3e6

    model_Q.set_param_hint('sigma',value=guess_sigma,
                                vary=True)
    #to guess the amplitude 
    #http://openafox.com/science/peak-function-derivations.html
    #Better sigma estimation is a better amplitude estimation
    if peak == 'max':
        guess_amp = voltage_max*guess_sigma*np.pi
    else: 
         guess_amp = -voltage_min*guess_sigma*np.pi       

    model_Q.set_param_hint('amplitude',value=guess_amp,
                                vary=True)
    #to guess the offset
    guess_offset = 0 
    model_Q.set_param_hint('offset',value=guess_offset,
                            vary=True)
    #guessed parameters
    guess_parameters = model_Q.make_params()

    #fit the model with the data and guessed parameters
    fit_res = model_Q.fit(data=data,frequency=frequency,params=guess_parameters)

    #get the values for output.
    f0 = fit_res.best_values['center']
    BW = (fit_res.best_values['sigma']*2)
    Q = abs(f0/BW)
    offset = fit_res.best_values['offset']

    if plot == True:
        #plot the fitted curve
        dummy_frequencies = np.linspace(np.amin(frequency),np.amax(frequency),101)
        fit_fine = resonator_peak(dummy_frequencies,**fit_res.best_values)
        fig = plt.figure(figsize=(8, 6))
        ax = plt.subplot()
        ax.plot(frequency,20*np.log(data),'o',label='Data') #plot data in dB
        ax.plot(dummy_frequencies,20*np.log(fit_fine),'r-', label='Lorentzian fit') #plot data fit in dB
        ax.set_ylabel('S21 (dB)')
        ax.set_xlabel('Frequency (GHz)')
        ax.legend()
        print("\nf0: {:.4f} GHz\n"
              "\nQ: {:.4f}\n" 
              "\nBW: {:.4f} MHz\n"
              "\noffset: {:.4f} \n".format(f0/1e9, Q, BW/1e6, offset))
        if report == True:
            fig.savefig('Lorentzian fit.pdf',format='pdf')

    if modelinfo == True:
        fig,ax = plt.subplots(1,1,figsize=(8,3))
        print("\nGuess_parameters:\n")
        print(guess_parameters)
        print("\nFit report:\n")
        print(fit_res.fit_report())
        print("\nFit best values:\n")
        print(fit_res.best_values)
        fit_res.plot_fit(show_init=True)

    return f0, BW, Q, offset

#Type B fittings

def circle_fit_transmission(zdata, freq, guess_BW = 0, guess_kc = 20e6, guess_offset = 0.0025, plot = True, report = True, modelinfo = True):
    """
    Fits a TRANSMISSION resonator into given zdata (complex).
    Input: 
        zdata (ndarray(np.complex128)), complex numbers with R and I values of the measurement.
        freq (ndarray), x-axis frequency values.
        guess_BW (float), guess of the BW.
        guess_kc (int), guess of the parameter kc/2*pi in MHz (this is the coupling)
        guess_offset (int), guess of the parameter offset
        plot (Boolean), when true plots the data with the fit
        report (Boolean), when true stores the plot in a .pdf file
        modelinfo (Boolean), when true prints the fit report as well as the initial guess
    Output: 
        f0 (float), fitted center frequency of the resonator.
        BW (float), 3dB bandwidth of the resonator.
        Q (float), Loaded Quality factor of the resonator.
        Qi (float), Internal Quality factor of the resonator.
        kl (float), kl =BW expressed as angular frequency (direct result of the fit). This kl is kl/2*pi to express it in MHz.
        kc (float), coupling expressed as angular frequency (direct result of the fit). This kc is kc/2*pi to express it in MHz.
        offset (float), offset 

    """
    frequency = freq
    data = zdata
    voltage_max_i = np.argmax(np.abs(data))

    def linear_resonator(frequency, center, kc, kl,offset,amp):
        w0=center*2*np.pi
        return amp*((2*kc)/(kl-2*1j*(frequency*2*np.pi-w0))+offset)
    model_complex = lmfit.Model(linear_resonator)

    #Guess parameters
    #to guess center
    guess_center = frequency[voltage_max_i] #Argmax = Returns the indices of the maximum values along an axis.
    model_complex.set_param_hint('center',value=guess_center,vary=True)
    
    #guess kl
    if guess_BW == 0:
        fmin = min(frequency)
        fmax = max(frequency)
        Q_min = 0.1 * (guess_center/(fmax-fmin))  # (minimum Q is a 0.1 than the measured (maximum) BW) assume the user isn't trying to fit just a small part of a resonance curve
        delta_f = np.diff(frequency)  # assume f is sorted
        min_delta_f = delta_f[delta_f > 0].min()
        Q_max = guess_center/min_delta_f  # (max Q is the minimum BW which is the step)assume data actually samples the resonance 
        Q_guess = np.sqrt(Q_min*Q_max)  # geometric mean, why not?
        guess_BW = guess_center/Q_guess

    guess_kl = guess_BW*2*np.pi 
    model_complex.set_param_hint('kl',value=guess_kl,vary=True, min= 0, max= 600e6)

    #guess kc
    model_complex.set_param_hint('kc',value=guess_kc*2*np.pi,vary=True, min= 0)

    #guess offset
    model_complex.set_param_hint('offset',value=guess_offset,vary=True)

    #guessed parameters
    guess_parameters = model_complex.make_params()

    #fit the model with the data and guessed parameters
    fit_res = model_complex.fit(data=data,frequency=frequency,params=guess_parameters)

    #get the values for the output.
    f0 = fit_res.best_values['center']
    kl = fit_res.best_values['kl']
    kc = fit_res.best_values['kc']
    BW = kl/(2*np.pi)
    Q = f0/BW
    Qc = f0/kc
    Qi = 1/(1/Q-1/Qc)
    offset = fit_res.best_values['offset']

    fit_s21 = model_complex.eval(params=fit_res.params, frequency=frequency)
    guess_s21 = model_complex.eval(params=guess_parameters, frequency=frequency)

    if plot == True:
            #plot the fitted curves (magnitude, phase and RI)
            fig = plt.figure(figsize=(12, 10))
            gs = gridspec.GridSpec(3, 2)
            gs.update(wspace=0.2)
            ax1 = plt.subplot(gs[0, 0], )
            ax2 = plt.subplot(gs[0, 1])
            ax3 = plt.subplot(gs[1:, :])

            #plot RI
            ax3.plot(data.real,data.imag,'o', label='Data') #add the data
            ax3.plot(fit_s21.real,fit_s21.imag, 'r-', label='best fit')
            ax3.plot(guess_s21.real, guess_s21.imag,'--', label='initial fit')
            ax3.legend()
            ax3.plot(0,0,'.') #add the origin (0,0) point
            ax3.set_aspect('equal')
            ax3.set_ylabel('Im')
            ax3.set_xlabel('Re')
            plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.9)

            #plot phase
            ax2.plot(freq, np.angle(data),'o', label = 'Data')
            ax2.plot(freq, np.angle(fit_s21), 'r-', label='best fit')
            ax2.plot(freq, np.angle(guess_s21), '--', label='initial fit')
            ax2.legend()
            ax2.set_ylabel('<S21 (rad)')
            ax2.set_xlabel('Frequency (GHz)')
            for tick in ax2.get_xticklabels():
                tick.set_rotation(45)

            #plot magnitude
            ax1.plot(frequency, 20*np.log(np.abs(data)), 'o', label = 'Data')
            ax1.plot(frequency, 20*np.log(np.abs(fit_s21)), 'r-', label='best fit')
            ax1.plot(frequency, 20*np.log(np.abs(guess_s21)), '--', label='initial fit')
            ax1.legend()
            ax1.set_ylabel('|S21| (dB)')
            ax1.set_xlabel('Frequency (GHz)')
            for tick in ax1.get_xticklabels():
                tick.set_rotation(45)

            plt.show()
            
            print("\nf0: {:.4f} GHz\n"
                "\nBW: {:.4f} MHz\n"
                "\nQ: {:.4f}\n" 
                "\nQi: {:.4f}\n"
                "\nkl: {:.4f}\n"
                "\nkc: {:.4f}\n".format(f0/1e9, BW/1e6, Q, Qi, kl, kc))
            if report == True:
                fig.savefig('circle fit for transmission resonator.pdf',format='pdf')

    if modelinfo == True:
        fig,ax = plt.subplots(1,1,figsize=(8,3))
        print("\nGuess_parameters:\n")
        print(guess_parameters)
        print("\nFit report:\n")
        print(fit_res.fit_report())
        print("\nFit best values:\n")
        print(fit_res.best_values)
        fit_res.plot_fit(show_init=True)
        
    return f0, BW, Q, Qi, kl, kc, offset

def circle_fit_transmission_mag(zdata, freq, guess_BW = 0, guess_kc = 20e6, guess_offset = 0.0025, plot = True, report = True, modelinfo = True):
    """
    Fits a TRANSMISSION resonator into given zdata magnitude.
    Input: 
        zdata (ndarray(np.complex128)), complex numbers with R and I values of the measurement.
        freq (ndarray), x-axis frequency values.
        guess_BW (float), guess of the BW.
        guess_kc (int), guess of the parameter kc/2*pi in MHz (this is the coupling)
        guess_offset (int), guess of the parameter offset
        plot (Boolean), when true plots the data with the fit
        report (Boolean), when true stores the plot in a .pdf file
        modelinfo (Boolean), when true prints the fit report as well as the initial guess
    Output: 
        f0 (float), fitted center frequency of the resonator.
        BW (float), 3dB bandwidth of the resonator.
        Q (float), Loaded Quality factor of the resonator.
        Qi (float), Internal Quality factor of the resonator.
        kl (float), kl =BW expressed as angular frequency (direct result of the fit). This kl is kl/2*pi to express it in MHz.
        kc (float), coupling expressed as angular frequency (direct result of the fit). This kc is kc/2*pi to express it in MHz.
        offset (float), offset 

    """
    frequency = freq
    data = np.abs(zdata)
    voltage_max_i = np.argmax(np.abs(data))

    def linear_resonator(frequency, center, kc, kl,offset):
        w0=center*2*np.pi
        return np.abs((2*kc)/(kl-2*1j*(frequency*2*np.pi-w0))+offset)
    model_complex = lmfit.Model(linear_resonator)

    #Guess parameters
    #to guess center
    guess_center = frequency[voltage_max_i] #Argmax = Returns the indices of the maximum values along an axis.
    model_complex.set_param_hint('center',value=guess_center,vary=True)
    
    #guess kl
    if guess_BW == 0:
        fmin = min(frequency)
        fmax = max(frequency)
        Q_min = 0.1 * (guess_center/(fmax-fmin))  # (minimum Q is a 0.1 than the measured (maximum) BW) assume the user isn't trying to fit just a small part of a resonance curve
        delta_f = np.diff(frequency)  # assume f is sorted
        min_delta_f = delta_f[delta_f > 0].min()
        Q_max = guess_center/min_delta_f  # (max Q is the minimum BW which is the step)assume data actually samples the resonance 
        Q_guess = np.sqrt(Q_min*Q_max)  # geometric mean, why not?
        guess_BW = guess_center/Q_guess

    guess_kl = guess_BW*2*np.pi 
    model_complex.set_param_hint('kl',value=guess_kl,vary=True, min= 0, max= 600e6)

    #guess kc
    model_complex.set_param_hint('kc',value=guess_kc*2*np.pi,vary=True, min= 0)

    #guess offset
    model_complex.set_param_hint('offset',value=guess_offset,vary=True)

    #guessed parameters
    guess_parameters = model_complex.make_params()

    #fit the model with the data and guessed parameters
    fit_res = model_complex.fit(data=data,frequency=frequency,params=guess_parameters)

    #get the values for the output.
    f0 = fit_res.best_values['center']
    kl = fit_res.best_values['kl']
    kc = fit_res.best_values['kc']
    BW = kl/(2*np.pi)
    Q = f0/BW
    Qc = f0/kc
    Qi = 1/(1/Q-1/Qc)
    offset = fit_res.best_values['offset']
    
    fit_s21 = model_complex.eval(params=fit_res.params, frequency=frequency)
    guess_s21 = model_complex.eval(params=guess_parameters, frequency=frequency)

    if plot == True:
            #plot the fitted curves (magnitude)
            fig = plt.figure(figsize=(8, 6))
            ax1 = plt.subplot()

            #plot magnitude
            ax1.plot(frequency, 20*np.log(np.abs(data)), 'o', label = 'Data')
            ax1.plot(frequency, 20*np.log(np.abs(fit_s21)), 'r-', label='best fit')
            ax1.legend()


            plt.show()
            
            print("\nf0: {:.4f} GHz\n"
                "\nBW: {:.4f} MHz\n"
                "\nQ: {:.4f}\n" 
                "\nQi: {:.4f}\n"
                "\nkl: {:.4f}\n"
                "\nkc: {:.4f}\n".format(f0/1e9, BW/1e6, Q, Qi, kl, kc))
            if report == True:
                fig.savefig('circle fit for transmission resonator.pdf',format='pdf')

    if modelinfo == True:
        fig,ax = plt.subplots(1,1,figsize=(8,3))
        print("\nGuess_parameters:\n")
        print(guess_parameters)
        print("\nFit report:\n")
        print(fit_res.fit_report())
        print("\nFit best values:\n")
        print(fit_res.best_values)
        fit_res.plot_fit(show_init=True)
        
    return f0, BW, Q, Qi, kl, kc, offset

#Type C fittings

def circle_fit_hanger(zdata, freq, plot = True, report = True, modelinfo = True, *args, **kwargs):
        """
        Fits a HANGER resonator into given zdata (complex).
        Input: 
            zdata (ndarray(np.complex128)), complex numbers with R and I values of the measurement.
            freq (ndarray), x-axis frequency values.
            guess_Qe (int), guess of the parameter Qe  (external quality factor (coupling))
            plot (Boolean), when true plots the data with the fit
            report (Boolean), when true stores the plot in a .pdf file
            modelinfo (Boolean), when true prints the fit report as well as the initial guess
        Output: 
            f0 (float), fitted center frequency of the resonator.
            BW (float), 3dB bandwidth of the resonator.
            Q (float), Loaded Quality factor of the resonator.
            Qi (float), Internal Quality factor of the resonator.
            Qe (float), External Quality factor of the resonator.
            phi_0 (float), propagation delay from the sample.
            phi_v (float), propagation delay to the sample.
            theta (float), accounts for the cable delay related with the path length of the cables.
            A (float), Amplitude?
            slope (float), also called alpha in literature accounts for any slope in the background transmission surrounding the resonance frequency.


        """
        
        frequency = freq
        data = zdata

        voltage_min_i = np.argmin(np.abs(data))
        voltage_max_i = np.argmax(np.abs(data))
        voltage_min = np.min(np.abs(data))
        voltage_max = np.max(np.abs(data))


        def linear_resonator(frequency, center, Q, Qe, theta, A, slope, phi_v, phi_0, B=0):
            hanger_contribution = 1 - (((Q/Qe) * np.exp(1j*theta)) / (1 + 2j * Q * (frequency - center) / center))
            slope_corr = (1+slope*(frequency-center)/center)
            propagation_delay_corr = np.exp(1j * (phi_v * frequency + phi_0))
            return A * hanger_contribution * slope_corr * propagation_delay_corr + B
        
        model_complex = lmfit.Model(linear_resonator)

        #Guess parameters

        #to guess center
        guess_center = frequency[voltage_max_i] #Argmin = Returns the indices of the minimum values along an axis.
        model_complex.set_param_hint('center',value=guess_center,vary=True)

        #guess bW
        guess_BW = kwargs.get('guess_BW', 0.)

        #guess Q
        fmin = min(frequency)
        fmax = max(frequency)
        Q_min = 0.1 * (guess_center/(fmax-fmin))  # (minimum Q is a 0.1 than the measured (maximum) BW) assume the user isn't trying to fit just a small part of a resonance curve
        delta_f = np.diff(frequency)  # assume f is sorted
        min_delta_f = delta_f[delta_f > 0].min()
        Q_max = guess_center/min_delta_f  # (max Q is the minimum BW which is the step)assume data actually samples the resonance reasonably
        if guess_BW == 0:
            guess_Q = np.sqrt(Q_min*Q_max) #I prefer to input the Q from the Lorentzian #np.sqrt(Q_min*Q_max)  # geometric mean, why not?
        else: 
            guess_Q = guess_center/guess_BW
        model_complex.set_param_hint('Q',value=guess_Q,vary=True, min=Q_min , max=Q_max)

        #guess Qe
        guess_Qe = kwargs.get('guess_Qe', abs(guess_Q/abs(1-(voltage_min))))
        Qe_min = kwargs.get('Qe_min', 1e1)
        Qe_max = kwargs.get('Qe_max', 1e6)
        model_complex.set_param_hint('Qe',value=guess_Qe,vary=True, min=Qe_min , max=Qe_max)

        #guess theta
        theta_guess = 0
        theta_min = -np.pi/2
        theta_max = np.pi/2
        model_complex.set_param_hint('theta',value=theta_guess,vary=True, min=theta_min , max=theta_max)

        #guess A
        measured_powers = (np.abs(data))
        A_min = min((np.abs(data)))#np.abs(data[0])
        A_max = max((np.abs(data)))
        A_guess = A_max-A_min#max(measured_powers) - min(measured_powers)
        model_complex.set_param_hint('A',value=A_guess,vary=True)
        
        #guess slope [alpha] (this is the slope in the background transmission surrounding the resonance)
        x_values = [frequency[0], frequency[len(frequency)-1]] #the same = f[-1]]
        y_values = [abs(data[0]), abs(data[len(data)-1])]
        slopef = linregress(x_values, y_values)
        slope_guess = slopef.slope
        model_complex.set_param_hint('slope',value=slope_guess,vary=True)
        
        #guess phi_0
        model_complex.set_param_hint('phi_0',value=288.,vary=True)
        
        #guess phi_v
        model_complex.set_param_hint('phi_v',value=0.,vary=True)

        #guess phi_v
        model_complex.set_param_hint('B',value=0.,vary=True)

        #guessed parameters
        guess_parameters = model_complex.make_params()

        #fit the model with the data and guessed parameters
        fit_res = model_complex.fit(data=data,frequency=frequency,params=guess_parameters)

        #get the values for postprocessing and for legend.
        f0 = fit_res.best_values['center']
        Q = fit_res.best_values['Q']
        Qe = fit_res.best_values['Qe']
        theta = fit_res.best_values['theta']
        # Qi =1./((1./Q) - np.real(1./Qe))
        Qi =1./((1./Q) - abs(np.cos(-theta)/Qe))
        A = fit_res.best_values['A']
        slope = fit_res.best_values['slope']
        phi_v = fit_res.best_values['phi_v']
        phi_0 = fit_res.best_values['phi_0']
        BW = f0/Q
        kc = f0/Qe*1e-6
        #fit_res.best_values.update({'Q_i':Qi})
        
        fit_s21 = model_complex.eval(params=fit_res.params, frequency=frequency)
        guess_s21 = model_complex.eval(params=guess_parameters, frequency=frequency)

        if plot == True:
                #plot the fitted curves (magnitude, phase and RI)
                fig = plt.figure(figsize=(14, 4))
                gs = gridspec.GridSpec(1, 3)
                gs.update(wspace=0.2)
                ax1 = plt.subplot(gs[0, 0])
                ax2 = plt.subplot(gs[0, 1])
                ax3 = plt.subplot(gs[0, 2])

                #plot RI
                ax3.plot(data.real,data.imag,'o', label='Data') #add the data
                ax3.plot(fit_s21.real,fit_s21.imag, 'r-', label='best fit')
                #ax3.plot(guess_s21.real, guess_s21.imag,'--', label='initial fit')
                ax3.legend()
                ax3.plot(0,0,'.') #add the origin (0,0) point
                ax3.set_aspect('equal')
                ax3.set_ylabel('Im')
                ax3.set_xlabel('Re')
                plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.9)

                #plot phase
                ax2.plot(freq, np.angle(data),'o', label = 'Data')
                ax2.plot(freq, np.angle(fit_s21), 'r-', label='best fit')
                #ax2.plot(freq, np.angle(guess_s21), '--', label='initial fit')
                ax2.legend()
                ax2.set_ylabel('<S21 (rad)')
                ax2.set_xlabel('Frequency (GHz)')
                for tick in ax2.get_xticklabels():
                    tick.set_rotation(45)

                #plot magnitude
                ax1.plot(frequency, 20*np.log(np.abs(data)), 'o', label = 'Data')
                ax1.plot(frequency, 20*np.log(np.abs(fit_s21)), 'r-', label='best fit')
                #ax1.plot(frequency, 20*np.log(np.abs(guess_s21)), '--', label='initial fit')
                ax1.legend()
                ax1.set_ylabel('|S21| (dB)')
                ax1.set_xlabel('Frequency (GHz)')
                for tick in ax1.get_xticklabels():
                    tick.set_rotation(45)

                plt.show()
                
                print("\nf0: {:.4f} GHz\n"
                    "\nBW: {:.4f} MHz\n"
                    "\nQ: {:.4f}\n" 
                    "\nQi: {:.4f}\n"
                    "\nQe: {:.4f}\n"
                    "\nkc: {:.4f} MHz\n"
                    "\ntheta: {:.4f}\n"
                    "\nA: {:.4f}\n"
                    "\nslope: {:.4f}\n"
                    "\nphi_v: {:.4f}\n"
                    "\nphi_0: {:.4f}\n".format(f0/1e9, BW/1e6, Q, Qi, Qe, kc,theta, A, slope, phi_v, phi_0))
                if report == True:
                    fig.savefig('circle fit for hanger resonator.pdf',format='pdf')

        if modelinfo == True:
            #fig,ax = plt.subplots(1,1,figsize=(8,3))
            print("\nGuess_parameters:\n")
            print(guess_parameters)
            print("\nFit report:\n")
            print(fit_res.fit_report())
            print("\nFit best values:\n")
            print(fit_res.best_values)
            #fit_res.plot_fit(show_init=True)
            
        return f0, BW, Q, Qi, Qe, kc, theta, A, slope, phi_v, phi_0, model_complex


def circle_fit_hanger_Q(zdata, freq, guess_BW = 0, Qe = 20e6, plot = True, report = True, modelinfo = True, Q = 1):
        """
        Fits a HANGER resonator into given zdata (complex).
        Input: 
            zdata (ndarray(np.complex128)), complex numbers with R and I values of the measurement.
            freq (ndarray), x-axis frequency values.
            guess_Qe (int), guess of the parameter Qe  (external quality factor (coupling))
            plot (Boolean), when true plots the data with the fit
            report (Boolean), when true stores the plot in a .pdf file
            modelinfo (Boolean), when true prints the fit report as well as the initial guess
        Output: 
            f0 (float), fitted center frequency of the resonator.
            BW (float), 3dB bandwidth of the resonator.
            Q (float), Loaded Quality factor of the resonator.
            Qi (float), Internal Quality factor of the resonator.
            Qe (float), External Quality factor of the resonator.
            phi_0 (float), propagation delay from the sample.
            phi_v (float), propagation delay to the sample.
            theta (float), accounts for the cable delay related with the path length of the cables.
            A (float), Amplitude?
            slope (float), also called alpha in literature accounts for any slope in the background transmission surrounding the resonance frequency.


        """
        
        frequency = freq
        data = zdata

        voltage_min_i = np.argmin(np.abs(data))
        voltage_max_i = np.argmax(np.abs(data))
        voltage_min = np.min(np.abs(data))
        voltage_max = np.max(np.abs(data))


        def linear_resonator(frequency, center, Qi, Qe, theta, A, slope, phi_v, phi_0):


            hanger_contribution = 1 - (((Qi/(Qe+Qi)) * np.exp(1j*theta)) / (1 + 2j * (Qi/(Qe+Qi)) * (frequency - center) / center))
            slope_corr = (1+slope*(frequency-center)/center)
            propagation_delay_corr = np.exp(1j * (phi_v * frequency + phi_0))
            return A * hanger_contribution * slope_corr * propagation_delay_corr
        model_complex = lmfit.Model(linear_resonator)

        #Guess parameters
        #to guess center
        guess_center = frequency[voltage_max_i] #Argmin = Returns the indices of the minimum values along an axis.
        model_complex.set_param_hint('center',value=guess_center,vary=True)

        #guess Q
        fmin = min(frequency)
        fmax = max(frequency)
        Q_min = 0.1 * (guess_center/(fmax-fmin))  # (minimum Q is a 0.1 than the measured (maximum) BW) assume the user isn't trying to fit just a small part of a resonance curve
        delta_f = np.diff(frequency)  # assume f is sorted
        min_delta_f = delta_f[delta_f > 0].min()
        Q_max = guess_center/min_delta_f  # (max Q is the minimum BW which is the step)assume data actually samples the resonance reasonably
        if guess_BW == 0:
            guess_Q = np.sqrt(Q_min*Q_max) #I prefer to input the Q from the Lorentzian #np.sqrt(Q_min*Q_max)  # geometric mean, why not?
        else: 
            guess_Q = guess_center/guess_BW
        model_complex.set_param_hint('Qi',value=guess_Q,vary=True, min=Q_min , max=Q_max)

        #guess Qe
        guess_Qe = abs(guess_Q/abs(1-(voltage_min))) #Q_guess/(1-np.abs(data[argmin_s21]))  #From mazin p29 or S.Asaad p15 eq 2.3. guess_Qe = 5e3
        Qe_min = 1
        Qe_max = 8e5
        model_complex.set_param_hint('Qe',value=guess_Qe,vary=True, min=Qe_min , max=Qe_max)

        #guess theta
        theta_guess = 0
        theta_min = -np.pi/2
        theta_max = np.pi/2
        model_complex.set_param_hint('theta',value=theta_guess,vary=True, min=theta_min , max=theta_max)

        #guess A
        measured_powers = (np.abs(data))
        A_min = min((np.abs(data)))#np.abs(data[0])
        A_max = max((np.abs(data)))
        A_guess = A_max-A_min#max(measured_powers) - min(measured_powers)
        model_complex.set_param_hint('A',value=A_guess,vary=True)
        
        #guess slope [alpha] (this is the slope in the background transmission surrounding the resonance)
        x_values = [frequency[0], frequency[len(frequency)-1]] #the same = f[-1]]
        y_values = [abs(data[0]), abs(data[len(data)-1])]
        slopef = linregress(x_values, y_values)
        slope_guess = slopef.slope
        model_complex.set_param_hint('slope',value=slope_guess,vary=True)
        
        #guess phi_0
        model_complex.set_param_hint('phi_0',value=288,vary=True)
        
        #guess phi_v
        model_complex.set_param_hint('phi_v',value=0,vary=True)

        #guessed parameters
        guess_parameters = model_complex.make_params()

        #fit the model with the data and guessed parameters
        fit_res = model_complex.fit(data=data,frequency=frequency,params=guess_parameters)

        #get the values for postprocessing and for legend.
        f0 = fit_res.best_values['center']
        Qi = fit_res.best_values['Qi']
        Qe = fit_res.best_values['Qe']
        theta = fit_res.best_values['theta']
        # Qi =1./((1./Q) - np.real(1./Qe))
        Q =1./(1./Qi + 1./Qe)
        A = fit_res.best_values['A']
        slope = fit_res.best_values['slope']
        phi_v = fit_res.best_values['phi_v']
        phi_0 = fit_res.best_values['phi_0']
        BW = f0/Q
        kc = f0/Qe*1e-6
        #fit_res.best_values.update({'Q_i':Qi})
        
        fit_s21 = model_complex.eval(params=fit_res.params, frequency=frequency)
        guess_s21 = model_complex.eval(params=guess_parameters, frequency=frequency)

        if plot == True:
                #plot the fitted curves (magnitude, phase and RI)
                fig = plt.figure(figsize=(14, 4))
                gs = gridspec.GridSpec(1, 3)
                gs.update(wspace=0.2)
                ax1 = plt.subplot(gs[0, 0])
                ax2 = plt.subplot(gs[0, 1])
                ax3 = plt.subplot(gs[0, 2])

                #plot RI
                ax3.plot(data.real,data.imag,'o', label='Data') #add the data
                ax3.plot(fit_s21.real,fit_s21.imag, 'r-', label='best fit')
                #ax3.plot(guess_s21.real, guess_s21.imag,'--', label='initial fit')
                ax3.legend()
                ax3.plot(0,0,'.') #add the origin (0,0) point
                ax3.set_aspect('equal')
                ax3.set_ylabel('Im')
                ax3.set_xlabel('Re')
                plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.9)

                #plot phase
                ax2.plot(freq, np.angle(data),'o', label = 'Data')
                ax2.plot(freq, np.angle(fit_s21), 'r-', label='best fit')
                #ax2.plot(freq, np.angle(guess_s21), '--', label='initial fit')
                ax2.legend()
                ax2.set_ylabel('<S21 (rad)')
                ax2.set_xlabel('Frequency (GHz)')
                for tick in ax2.get_xticklabels():
                    tick.set_rotation(45)

                #plot magnitude
                ax1.plot(frequency, 20*np.log(np.abs(data)), 'o', label = 'Data')
                ax1.plot(frequency, 20*np.log(np.abs(fit_s21)), 'r-', label='best fit')
                #ax1.plot(frequency, 20*np.log(np.abs(guess_s21)), '--', label='initial fit')
                ax1.legend()
                ax1.set_ylabel('|S21| (dB)')
                ax1.set_xlabel('Frequency (GHz)')
                for tick in ax1.get_xticklabels():
                    tick.set_rotation(45)

                plt.show()
                
                print("\nf0: {:.4f} GHz\n"
                    "\nBW: {:.4f} MHz\n"
                    "\nQ: {:.4f}\n" 
                    "\nQi: {:.4f}\n"
                    "\nQe: {:.4f}\n"
                    "\nkc: {:.4f} MHz\n"
                    "\ntheta: {:.4f}\n"
                    "\nA: {:.4f}\n"
                    "\nslope: {:.4f}\n"
                    "\nphi_v: {:.4f}\n"
                    "\nphi_0: {:.4f}\n".format(f0/1e9, BW/1e6, Q, Qi, Qe, kc,theta, A, slope, phi_v, phi_0))
                if report == True:
                    fig.savefig('circle fit for hanger resonator.pdf',format='pdf')

        if modelinfo == True:
            fig,ax = plt.subplots(1,1,figsize=(8,3))
            print("\nGuess_parameters:\n")
            print(guess_parameters)
            print("\nFit report:\n")
            print(fit_res.fit_report())
            print("\nFit best values:\n")
            print(fit_res.best_values)
            fit_res.plot_fit(show_init=True)
            
        return f0, BW, Q, Qi, Qe, kc, theta, A, slope, phi_v, phi_0


def circle_fit_hanger_mag(zdata, freq, plot = True, report = True, modelinfo = True, *args, **kwargs):
        """
        Fits a HANGER resonator into given zdata magnitude.
        Input: 
            zdata (ndarray(np.complex128)), complex numbers with R and I values of the measurement.
            freq (ndarray), x-axis frequency values.
            guess_Qe (int), guess of the parameter Qe  (external quality factor (coupling))
            plot (Boolean), when true plots the data with the fit
            report (Boolean), when true stores the plot in a .pdf file
            modelinfo (Boolean), when true prints the fit report as well as the initial guess
        Output: 
            f0 (float), fitted center frequency of the resonator.
            BW (float), 3dB bandwidth of the resonator.
            Q (float), Loaded Quality factor of the resonator.
            Qi (float), Internal Quality factor of the resonator.
            Qe (float), External Quality factor of the resonator.
            phi_0 (float), propagation delay from the sample.
            phi_v (float), propagation delay to the sample.
            theta (float), accounts for the cable delay related with the path length of the cables.
            A (float), Amplitude?
            slope (float), also called alpha in literature accounts for any slope in the background transmission surrounding the resonance frequency.


        """
        
        frequency = freq
        data = np.abs(zdata)

        voltage_min_i = np.argmin(np.abs(data))
        voltage_max_i = np.argmax(np.abs(data))
        voltage_min = np.min(np.abs(data))
        voltage_max = np.max(np.abs(data))


        def linear_resonator(frequency, center, Q, Qe, theta, A, slope, phi_v, phi_0):
            hanger_contribution = 1 - (((Q/Qe) * np.exp(1j*theta)) / (1 + 2j * Q * (frequency - center) / center))
            slope_corr = (1+slope*(frequency-center)/center)
            propagation_delay_corr = np.exp(1j * (phi_v * frequency + phi_0))
            return np.abs(A * hanger_contribution * slope_corr * propagation_delay_corr)
        model_complex = lmfit.Model(linear_resonator)

        #Guess parameters
        #to guess center
        guess_center = frequency[voltage_min_i] #Argmin = Returns the indices of the minimum values along an axis.
        model_complex.set_param_hint('center',value=guess_center,vary=True)

        # Guess BW
        guess_BW = kwargs.get('guess_BW', 0.)

        #guess Q
        fmin = min(frequency)
        fmax = max(frequency)
        Q_min = 0.1 * (guess_center/(fmax-fmin))  # (minimum Q is a 0.1 than the measured (maximum) BW) assume the user isn't trying to fit just a small part of a resonance curve
        delta_f = np.diff(frequency)  # assume f is sorted
        min_delta_f = delta_f[delta_f > 0].min()
        Q_max = guess_center/min_delta_f  # (max Q is the minimum BW which is the step)assume data actually samples the resonance reasonably
        if guess_BW == 0:
            guess_Q = np.sqrt(Q_min*Q_max) #I prefer to input the Q from the Lorentzian #np.sqrt(Q_min*Q_max)  # geometric mean, why not?
        else: 
            guess_Q = guess_center/guess_BW
        model_complex.set_param_hint('Q',value=guess_Q,vary=True, min=Q_min , max=Q_max)

        #guess Qe
        guess_Qe = kwargs.get('guess_Qe', abs(guess_Q/abs(1-(voltage_min))))
        Qe_min = kwargs.get('Qe_min', 1e3)
        Qe_max = kwargs.get('Qe_max', 1e6)
        model_complex.set_param_hint('Qe',value=guess_Qe,vary=True, min=Qe_min , max=Qe_max)

        #guess theta
        theta_guess = 0
        theta_min = -np.pi/2
        theta_max = np.pi/2
        model_complex.set_param_hint('theta',value=theta_guess,vary=True, min=theta_min , max=theta_max)

        #guess A
        measured_powers = (np.abs(data))
        A_min = min((np.abs(data)))#np.abs(data[0])
        A_max = max((np.abs(data)))
        A_guess = A_max-A_min#max(measured_powers) - min(measured_powers)
        model_complex.set_param_hint('A',value=A_guess,vary=True)
        
        #guess slope [alpha] (this is the slope in the background transmission surrounding the resonance)
        x_values = [frequency[0], frequency[len(frequency)-1]] #the same = f[-1]]
        y_values = [abs(data[0]), abs(data[len(data)-1])]
        slopef = linregress(x_values, y_values)
        slope_guess = slopef.slope
        model_complex.set_param_hint('slope',value=slope_guess,vary=True)
        
        #guess phi_0
        model_complex.set_param_hint('phi_0',value=288,vary=True)
        
        #guess phi_v
        model_complex.set_param_hint('phi_v',value=0,vary=True)

        #guessed parameters
        guess_parameters = model_complex.make_params()

        #fit the model with the data and guessed parameters
        fit_res = model_complex.fit(data=data,frequency=frequency,params=guess_parameters)

        #get the values for postprocessing and for legend.
        f0 = fit_res.best_values['center']
        Q = fit_res.best_values['Q']
        Qe = fit_res.best_values['Qe']
        Qi =1./((1./Q) - np.real(1./Qe))
        theta = fit_res.best_values['theta']
        A = fit_res.best_values['A']
        slope = fit_res.best_values['slope']
        phi_v = fit_res.best_values['phi_v']
        phi_0 = fit_res.best_values['phi_0']
        BW = f0/Q
        kc = f0/(np.abs(Qe))
        #fit_res.best_values.update({'Q_i':Qi})
        
        fit_s21 = model_complex.eval(params=fit_res.params, frequency=frequency)
        guess_s21 = model_complex.eval(params=guess_parameters, frequency=frequency)

        if plot == True:
                #plot the fitted curves (magnitude)
                fig = plt.figure(figsize=(8, 6))
                ax1 = plt.subplot()

                #plot magnitude
                ax1.plot(frequency, 20*np.log(np.abs(data)), 'o', label = 'Data')
                ax1.plot(frequency, 20*np.log(np.abs(fit_s21)), 'r-', label='best fit')
                ax1.legend()


                plt.show()
                
                print("\nf0: {:.4f} GHz\n"
                    "\nBW: {:.4f} MHz\n"
                    "\nQ: {:.4f}\n" 
                    "\nQi: {:.4f}\n"
                    "\nQe: {:.4f}\n"
                    "\nkc: {:.4f}\n"
                    "\ntheta: {:.4f}\n"
                    "\nA: {:.4f}\n"
                    "\nslope: {:.4f}\n"
                    "\nphi_v: {:.4f}\n"
                    "\nphi_0: {:.4f}\n".format(f0/1e9, BW/1e6, Q, Qi, Qe, kc,theta, A, slope, phi_v, phi_0))
                if report == True:
                    fig.savefig('circle fit magnitude for hanger resonator.pdf',format='pdf')
            
        return f0, BW, Q, Qi, Qe, kc, theta, A, slope, phi_v, phi_0


#Batch fit

def complex_fit_amplitude_batch_fit(label, tau, theta, c, report = False, modelinfo = False, plot = False, datatype = "zdata"):
    material = label[:2]
    attenuation = label[-2:]
    filename = f'data/*{material}*-{attenuation}.S2P'
    res_files = glob.glob(filename)
    if not res_files:
        filename = f'data/*{material}*-{attenuation}.txt'
        res_files = glob.glob(filename)
    
    f0_values = []
    working_files = []
    Q_values = []
    Qi_values = []
    res_files
    for file in res_files:
        try:
            s21_complex, zdata, zdata_r, freq = file2data(file, tau, c, theta)
            if datatype == "zdata":
                f0, BW, Q, offset = lorentzian_fit(zdata,freq,report = report, modelinfo = modelinfo, plot = plot ) 
                f0, BW, Q, Qi, Qe, kc, theta, A, slope, phi_v, phi_0 = circle_fit_hanger(zdata, freq, report = report, modelinfo = modelinfo, plot = plot, guess_BW=BW) 
            elif datatype == "s21_complex":
                f0, BW, Q, offset = lorentzian_fit(s21_complex,freq,report = report, modelinfo = modelinfo, plot = plot ) 
                f0, BW, Q, Qi, Qe, kc, theta, A, slope, phi_v, phi_0 = circle_fit_hanger(s21_complex, freq, report = report, modelinfo = modelinfo, plot = plot, guess_BW=BW)
            elif datatype == "zdata_r":
                f0, BW, Q, offset = lorentzian_fit(zdata_r,freq,report = report, modelinfo = modelinfo, plot = plot ) 
                f0, BW, Q, Qi, Qe, kc, theta, A, slope, phi_v, phi_0 = circle_fit_hanger(zdata_r, freq, report = report, modelinfo = modelinfo, plot = plot, guess_BW=BW) 

            f0 = np.float(f0/1e9)
            f0_values.append(f0)
            Q_values.append(Q)
            Qi_values.append(Qi)
            working_files.append(file)
        except Exception as e:
            print(e)
            print(f"Problem with file - {file}")
            continue
    # create dataframe and write to csv
    res_files = pd.DataFrame(
        {
            "TUID": working_files,
            "f0": f0_values,
            "Q": Q_values,
            "Qi": Qi_values,
        }
    )
    res_files = res_files[res_files.f0 > 1]
    res_files = res_files[res_files.f0 < 8]
    filename2 = f'{label}_resonator_data.csv'
    res_files.to_csv(filename2)
    return res_files

def complex_fit_amplitude_batch_fit_mag(label, tau, theta, c, report = False, modelinfo = False, plot = False):
    material = label[:2]
    attenuation = label[-2:]
    filename = f'data/*{material}*-{attenuation}.S2P'
    res_files = glob.glob(filename)
    if not res_files:
        filename = f'data/*{material}*-{attenuation}.txt'
        res_files = glob.glob(filename)
    
    f0_values = []
    working_files = []
    Q_values = []
    Qi_values = []
    res_files
    for file in res_files:
        try:
            s21_complex, zdata, zdata_r, freq = file2data(file, tau, c, theta)
            f0, BW, Q, offset = lorentzian_fit(zdata,freq,report = report, modelinfo = modelinfo, plot = plot ) 
            f0, BW, Q, Qi, Qe, kc, theta, A, slope, phi_v, phi_0 = circle_fit_hanger_mag(zdata, freq, report = report, modelinfo = modelinfo, plot = plot, guess_BW=BW) 
            f0 = np.float(f0/1e9)
            f0_values.append(f0)
            Q_values.append(Q)
            Qi_values.append(Qi)
            working_files.append(file)
        except Exception as e:
            print(e)
            print(f"Problem with file - {file}")
            continue
    # create dataframe and write to csv
    res_files = pd.DataFrame(
        {
            "TUID": working_files,
            "f0": f0_values,
            "Q": Q_values,
            "Qi": Qi_values,
        }
    )
    res_files = res_files[res_files.f0 > 1]
    res_files = res_files[res_files.f0 < 8]
    filename2 = f'{label}_resonator_data.csv'
    res_files.to_csv(filename2)
    return res_files

def complex_fit_batch_hdf5(label, tau, theta, c, report = False, modelinfo = False, plot = False, datatype = "zdata"):
    f0_values = []
    Q_values = []
    Qi_values = []
    for freq in range(8):
        df = read_hdf5(label,freq,5) #for the frist frequency center (0) index, the power (5) index (corresponding to 5dBm).
        voltage = 10**(df.S21_dB/20)
        s21_complex = np.array(voltage*np.exp(1j*df.S21_ang*np.pi/180.),dtype=np.complex128) #dtype = # 128-bit complex floating-point number
        freq = df.frequency.values
        zdata = np.array(s21_complex*np.exp(1j*freq*tau*1e-9*2*np.pi),dtype=np.complex128) # This is the unwinded data.Add the electrical delay as exp^(j·w·tau)
        zdata_r = c*zdata*np.exp(1j*theta) #rotated zdata

        try:
            if datatype == "zdata":
                f0, BW, Q, offset = lorentzian_fit(zdata,freq,report = report, modelinfo = modelinfo, plot = plot ) 
                f0, BW, Q, Qi, Qe, kc, theta, A, slope, phi_v, phi_0 = circle_fit_hanger(zdata, freq, report = report, modelinfo = modelinfo, plot = plot, guess_BW=BW) 
            elif datatype == "s21_complex":
                f0, BW, Q, offset = lorentzian_fit(s21_complex,freq,report = report, modelinfo = modelinfo, plot = plot ) 
                f0, BW, Q, Qi, Qe, kc, theta, A, slope, phi_v, phi_0 = circle_fit_hanger(s21_complex, freq, report = report, modelinfo = modelinfo, plot = plot, guess_BW=BW)
            elif datatype == "zdata_r":
                f0, BW, Q, offset = lorentzian_fit(zdata_r,freq,report = report, modelinfo = modelinfo, plot = plot ) 
                f0, BW, Q, Qi, Qe, kc, theta, A, slope, phi_v, phi_0 = circle_fit_hanger(zdata_r, freq, report = report, modelinfo = modelinfo, plot = plot, guess_BW=BW) 

            f0 = float(f0/1e9)
            f0_values.append(f0)
            Q_values.append(Q)
            Qi_values.append(Qi)
        except Exception as e:
            print(e)
            print(f"Problem with freq index - {freq}")
            continue
    # create dataframe and write to csv
    res_files = pd.DataFrame(
        {
            "f0": f0_values,
            "Q": Q_values,
            "Qi": Qi_values,
        }
    )
    res_files = res_files[res_files.f0 > 1]
    res_files = res_files[res_files.f0 < 8]
    filename2 = f'resonator_data.csv'
    res_files.to_csv(filename2)
    return res_files

