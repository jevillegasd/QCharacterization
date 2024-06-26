import pyvisa
from time import sleep

rm = pyvisa.ResourceManager()

import numpy as np

class E5080B_driver():
    '''Driver for the VNA E5080B'''
    
    def __init__(self,address):
        self._address = address

    def send(self, cmd):
        self._inst.write(cmd)

    def query(self, cmd: str):
        rep = None
        if cmd.find('?') != -1:
            rep = self._inst.query(cmd)
        return rep

    def open(self):
        '''Open conenction to instrument using address given in the constructor'''
        self._inst = rm.open_resource(self._address)
        self.idn = self._inst.query("*IDN?")
        self.stop_rf()
        self._npoints = self.get_sweep_npoints()

    def reset_average(self):
        '''Reset average'''
        self._inst.write('SENS:AVER:CLE')

    def autoscale(self, trace:int = 1):
        '''Auto Scale'''
        self._inst.write(':DISP:WIND1:TRAC%1i:Y:AUTO'%trace)

    def start_averaging(self):
        self._inst.write('SENSe:AVERage:STATe 1')

    def stop_averaging(self):
        self._inst.write('SENSe:AVERage:STATe 0')


    def set_data_format(self,data_format):
        ''' 
        The format available are:
        MA - Linear Magnitude / degrees
        DB - Log Magnitude / degrees
        RI - Real / Imaginary
        '''

        if data_format != 'MA' and data_format != 'DB' and data_format != 'RI':
            raise ValueError("Format not available.")
        
        self._inst.write("MMEM:STOR:TRAC:FORM:SNP {}".format(data_format))

    def set_average_count(self,count):
        if count < 0:
            raise ValueError("The number of average count should be positive and integer.")
            
        self._inst.write('SENSe:AVERage:COUNT {}'.format(int(count)))

    def get_average_count(self):
        result = self._inst.query("SENSe:AVERage:COUNT?")

        return int(result[:-1])

    def stop_rf(self):
        self._inst.write("OUTPut:STATe 0") 

    def start_rf(self):
        self._inst.write("OUTPut:STATe 1") 

    def set_power(self,power_dBm, overule_power = False):
        '''Set the power in dBm'''

        if power_dBm > 15:
            raise ValueError("Power over instrument limit. Should  be between -100 and 15.")
        if power_dBm < -100:
            raise ValueError("Power below instrument limit. Should  be between -100 and 15.")

        if not overule_power and power_dBm > -5:
            raise ValueError("Power may be too high and may damage qubits. If you are sure, set overule_power to True.")

        self._inst.write("SOUR:POW1 {}".format(power_dBm))

    def set_sweep_npoints(self,npoints):
        '''Set the number of points for sweep.'''
        
        if npoints < 0:
            raise ValueError("Negative numbers are not allowed.")
        
        self._inst.write("SENSe:SWEep:POINts {}".format(int(npoints)))

        self._npoints = int(npoints)

    def get_sweep_npoints(self):
        '''Get the number of points for the sweep.'''

        result = self._inst.query("SENSe:SWEep:POINts?")

        return int(result[:-1])

    def set_span_frequency(self, span_MHz):
        '''Set the span frequency in MHz'''
        self._inst.write("SENSe:FREQuency:SPAN {}".format(int(span_MHz*1e6)))

    def get_span_frequency(self):
        result = self._inst.query("SENSe:FREQuency:SPAN?")
        return float(result[:-1])

    def set_center_frequency(self, center):
        self._inst.write("SENSe:FREQuency:CENTer {}".format(center))

    def get_center_frequency(self):
        result = self._inst.query("SENSe:FREQuency:CENTer?")
        return float(result[:-1])

    def set_stop_frequency(self, stop):
        self._inst.write("SENSe:FREQuency:STOP {}".format(stop))

    def get_stop_frequency(self):
        result = self._inst.query("SENSe:FREQuency:STOP?")
        return float(result[:-1])

    def set_start_frequency(self, start):
        self._inst.write("SENSe:FREQuency:STARt {}".format(start))

    def get_start_frequency(self):
        result = self._inst.query("SENSe:FREQuency:STARt?")
        return float(result[:-1])

    def get_freq_array(self):
        '''Return the data from the VNA. The result is big array with frequency, S11,S21,S12,S22. The first npoints is the frequency.'''
        data = self._inst.query("CALCulate:MEASure:DATA:SNP?")
        f = np.array(data.split(','),dtype=float)
        return f[:self._npoints]

    def get_data(self, meas = 'S21'):
        '''Return the data from the VNA. The result is big array with frequency, S11,S21,S12,S22. the third array of npoints is the S21 mag.
           and the fourth array of npoints is the S21 phase. '''
        data = self._inst.query("CALCulate:MEASure:DATA:SNP?")
        f = np.array(data.split(','),dtype=float)

        npoints = self._npoints
        
        if meas == 'S11':
            j = 1
            i = 2
        elif meas == 'S12':
            j = 5
            i = 6
        elif meas == 'S22':
            j = 7
            i = 8
        else:
            j = 3
            i = 4

        
        mag = f[npoints*j:npoints*(j+1)]    
        phase = f[npoints*i:npoints*(i+1)]

        return mag,phase


    
    def close(self):
        self.stop_rf()
        self._inst.close()


#Context manager handle to start the RF output
class run():
    def __init__(self,driver: E5080B_driver):
        self.driver = driver

    def __enter__(self):
        self.driver.start_rf()
        self.driver.reset_average()
        sleep(0.1)
        return self.driver

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.driver.stop_rf()