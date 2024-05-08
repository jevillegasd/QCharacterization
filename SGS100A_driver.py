import pyvisa
import numpy as np

rm = pyvisa.ResourceManager()

class SGS100A_driver():
    def __init__(self, address):
        self._address = address
        self.connected = False
        try:
            self._inst = rm.open_resource(self._address)
            self.idn = self._inst.query("*IDN?")
            self.stop_rf()
            self._inst.close()
            self.connected = True

        except:
            self._inst = False
            print('Could not connect to %s'%address)


    def __enter__(self):
        if self.connected:
            self.open()
            self.start_rf()

    def __exit__(self, ex_type, ex_value, ex_traceback):
        if self.connected:
            self.stop_rf()
            self.close()

    def open(self):
        '''Open conenction to instrument using address given in the constructor'''
        if self.connected:
            self._inst = rm.open_resource(self._address)

    def close(self):
        if self.connected:
            self.stop_rf()
            self._inst.close()

    def start_rf(self):
        self._inst.write(':OUTPut 1')
    
    def stop_rf(self):
        self._inst.write(':OUTPut 0')
    
    def set_RF(self, freq_mhz):
        self._inst.write(':FREQuency:CW {} MHz'.format(freq_mhz))
    
    def set_level(self,level_dbm):
        '''Set the power output in dBm'''
        self._inst.write(':POWer:POWer {}'.format(level_dbm))
    
    def setup(self, frequency_MHz:float = 1e3, power_dbm:float = 0.):
        if self.connected:
            self.open()
            self.set_level(power_dbm)
            self.set_RF(frequency_MHz)
            self.close()


# Extended functions

    def set_PM_source(self, source='INT'):
        '''Set the the pulse modulation source.'''
        '''<source> = INT | EXT'''
        self._inst.write(':PULM:SOURce {}'.format(source))

    def set_PM_trigger_mode(self, mode='AUTO'):
        '''Set the pulse modulation trigger mode'''
        self._inst.write(':PULM:TRIGger:MODE AUTO'.format(mode))

    def set_PM_mode(self, mode='SING'):
        '''Set the pulse modulation trigger mode'''
        '''<mode> = SINGle | DOUBle'''
        self._inst.write(':PULM:MODE {}'.format(mode))

    def set_PM_pol(self, pol='NORM'):
        '''Set the polarity of the external pulse modulation source'''
        self._inst.write(':PULM:POLarity {}'.format(pol))

    def set_PM_trigger_imp(self, imp = 'G10K'):
        '''Select the impedance for the external pulse modulation trigger input'''
        self._inst.write(':PULM:TRIGger:EXTernal:IMPedance {}'.format(imp))

    def set_PG_period(self,period=0.0):
        '''Set pulse period in us'''
        self._inst.write(':PULM:PERiod {} us'.format(period))

    def set_PG_width(self,width):
        '''Set pulse width in us'''
        self._inst.write(':PULM:WIDth {} us'.format(width))

    def set_PG_double_width(self,width):
        '''Set double pulse width in case of a double pulse'''
        self._inst.write(':PULM:DOUBle:WIDTh {}'.format(width))

    def set_PG_double_delay(self,delay):
        '''Set double pulse delay'''
        self._inst.write(':PULM:DOUBle:DELay {}'.format(delay))

    def set_PG_state(self,state=0):
        '''Set pulse generator output state'''
        self._inst.write(':PGENerator:OUTPut:STATe {}'.format(state))

    def set_PM_state(self,state=0):
        '''Set pulse modulator state'''
        self._inst.write(':PULM:STATe {}'.format(state))

    def start_PG(self):
        self.set_PG_state(1)
    
    def stop_PG(self):
        self.set_PG_state(0)

    def start_PM(self):
        self.set_PM_state(1)

    def stop_PM(self):
        self.set_PM_state(0)


    def setup_PulseMod(self,freq=4000, amp=-25, period_us = 10, width_us = 0.2):
        '''Configure Pulse Width Modulation'''

        # Configure the pulse modulation settings 
        self.set_RF(freq)
        self.set_level(amp)
        self.set_PM_source('INT')
        self.set_PM_trigger_mode('AUTO')
        self.set_PM_mode('SING')

        # Configure the pulse generation settings
        self.set_PG_period(period_us)
        self.set_PG_width(width_us)

    def start_PulseMod(self):
        ''' '''
        self.start_PM()
        self.start_PG()
        self.start_rf()

    def stop_PulseMod(self):
        ''' '''
        self.stop_rf()
        self.stop_PG()
        self.stop_PM()
        