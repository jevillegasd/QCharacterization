import numpy as np
import matplotlib.pyplot as plt

import lmfit
from lmfit import Model
import matplotlib.gridspec as gridspec


#Models
def hanger_resonator(x_data, f_c, Q, Qe, A, A0, slope, theta, phi_v, phi_v0):
            transmission = 1 - (  ((Q/Qe) * np.exp(1j*theta)) / (1 + 2j * Q * (x_data - f_c) / f_c))
            slope_corr = (1+slope*(x_data-f_c)/f_c)
            propagation_delay_corr = np.exp(1j * (phi_v * x_data + phi_v0))
            return A*transmission * slope_corr * propagation_delay_corr 

def inline_resonator(x_data, f_c, Q, Qe, A, A0, slope, theta, phi_v, phi_v0):
            transmission = (  ((Q/Qe) * np.exp(1j*theta)) / (1 + 2j * Q * (x_data - f_c) / f_c))
            slope_corr = (1+slope*(x_data-f_c)/f_c)
            propagation_delay_corr = np.exp(1j * (phi_v * x_data + phi_v0))
            return A*transmission * slope_corr * propagation_delay_corr 


def circle_fit(x_data,z, plot=False, *args, **kwargs):
    '''
    Circle fit for resonators in hanger mode
    '''
    
    model_complex = Model(hanger_resonator)

    # put z data in linear format
    format = kwargs.get('format', 'MA')
    if format == 'DB':
        norm_z = 10**(z/20)
    else:
        norm_z = z

    df = np.diff(x_data)[0]

    # Set parameters of the fit

    #Center frequemcy
    if ('f_c' in kwargs):  vary = True
    else: vary = True
    f_c =  kwargs.get('f_c', x_data[z.argmin()])  #Argmin = Returns the indices of the minimum values along an axis.
    model_complex.set_param_hint('f_c',
                                 value=f_c,
                                 min= kwargs.get('f_c_min',f_c-10e6),
                                 max= kwargs.get('f_c_max',f_c+10e6),
                                 vary= kwargs.get('f_c_vary', True))

    #bandwidth
    fwhm = kwargs.get('fwhm', (x_data.max()-x_data.min())/10)

    #Quality factor
    Q = kwargs.get('Q', f_c/fwhm)
    model_complex.set_param_hint('Q',
                                 value= Q,
                                 min=kwargs.get('Q_min', Q/2),
                                 max=kwargs.get('Q_max', f_c/df/4),
                                 vary=kwargs.get('Q_vary', True))

    #External quality factor
    model_complex.set_param_hint('Qe',
                                 value=kwargs.get('Qe', Q),
                                 min=kwargs.get('Qe_min', Q/2), 
                                 max=kwargs.get('Qe_max', 1e7),
                                 vary=kwargs.get('Qe_vary', True))

    #Constant amplitude, should be zero becasue we normalize
    model_complex.set_param_hint('A',
                                 value= kwargs.get('A',1),
                                 vary=  kwargs.get('A_vary',True),
                                 min=  kwargs.get('A_min',0),
                                 max=  kwargs.get('A_max',10))

    #Constant amplitude, should be zero becasue we normalize
    model_complex.set_param_hint('A0',
                                 min= kwargs.get('A0_min',0),
                                 max = kwargs.get('A0_max',1),
                                 value=kwargs.get('A0',0),
                                 vary=kwargs.get('A0_vary',False))  # By default a background amplification of 0 is taken

    #Slope correction
    model_complex.set_param_hint('slope',
                                 value=kwargs.get('slope',0.),
                                 min= kwargs.get('slope_min',-10),
                                 max = kwargs.get('slope_max',10),
                                 vary=kwargs.get('slope_vary',True))

    #tetha
    model_complex.set_param_hint('theta',
                                 value=kwargs.get('theta',0.),
                                 vary=kwargs.get('theta_vary', True), 
                                 min=-np.pi, 
                                 max=+np.pi)

    #propagation delay dispersion
    model_complex.set_param_hint('phi_v',
                                 value=kwargs.get('phi_v',0.0),
                                 min= kwargs.get('phi_v_min',-np.pi),
                                 max= kwargs.get('phi_v_max',np.pi),
                                 vary=kwargs.get('phi_v_vary', True))

    #porpagation delay constant
    model_complex.set_param_hint('phi_v0',
                                 value=kwargs.get('phi_v0',288.0),
                                 min=kwargs.get('phi_v0_min',-1e3),
                                 max=kwargs.get('phi_v0_max',1e3),
                                 vary=kwargs.get('phi_v0_vary', True))

    # run the fit
    guess_parameters = model_complex.make_params()
    fit_res = model_complex.fit(data=norm_z,x_data=x_data, params=guess_parameters)
    res_params = fit_res.best_values
    r_sqrd = fit_res.rsquared
   
    # fit_params['Qi'] =1./((1./fit_params['Q']) - abs(np.cos(-fit_params['theta'])/fit_params['Qe']))
    res_params['Qi'] =1./((1./res_params['Q']) - (np.cos(-res_params['theta'])/res_params['Qe']))
    res_params['fwhm'] = res_params['f_c']/res_params['Q']
    res_params['kc'] = res_params['f_c']/res_params['Qe']*1e-6

    fit_s21 = model_complex.eval(params=fit_res.params, x_data=x_data)
    initial_s21 = model_complex.eval(params=guess_parameters, x_data=x_data)
    
    if plot:
        plot_all(x_data, norm_z,fit_s21,initial_s21)

    return res_params, fit_res, r_sqrd

def plot_all(freq, data,fit_s21,initial_s21):
    #plot the fitted curves (magnitude, phase and RI)
    fig = plt.figure(figsize=(14, 4))
    gs = gridspec.GridSpec(1, 3)
    gs.update(wspace=0.2)
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])
    ax3 = plt.subplot(gs[0, 2])

    #plot RI
    ax3.plot(data.real,data.imag,'x', label='Data') #add the data
    ax3.plot(fit_s21.real,fit_s21.imag, 'r-', label='best fit')
    ax3.plot(initial_s21.real, initial_s21.imag, '--', label='initial fit')

    #ax3.plot(guess_s21.real, guess_s21.imag,'--', label='initial fit')
    ax3.legend()
    ax3.plot(0,0,'gx') #add the origin (0,0) point
    ax3.set_aspect('equal')
    ax3.set_ylabel('Im')
    ax3.set_xlabel('Re')
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.9)

    #plot phase
    ax2.plot(freq, np.angle(data),'o', label = 'Data')
    ax2.plot(freq, np.angle(fit_s21), 'r-', label='best fit')
    ax2.plot(freq, np.angle(initial_s21), '--', label='initial fit')

    #ax2.plot(freq, np.angle(guess_s21), '--', label='initial fit')
    ax2.legend()
    ax2.set_ylabel('<S21 (rad)')
    ax2.set_xlabel('Frequency (GHz)')
    for tick in ax2.get_xticklabels():
        tick.set_rotation(45)

    #plot magnitude
    ax1.plot(freq, 10*np.log(np.abs(data)), 'o', label = 'Data')
    ax1.plot(freq, 10*np.log(np.abs(fit_s21)), 'r-', label='best fit')
    ax1.plot(freq, 10*np.log(np.abs(initial_s21)), '--', label='initial fit')

    #ax1.plot(frequency, 20*np.log(np.abs(guess_s21)), '--', label='initial fit')
    ax1.legend()
    ax1.set_ylabel('|S21| (dB)')
    ax1.set_xlabel('Frequency (GHz)')
    for tick in ax1.get_xticklabels():
        tick.set_rotation(45)

    plt.show()

