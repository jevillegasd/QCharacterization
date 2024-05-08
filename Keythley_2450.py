# ----------------------------------------------------------------------------
# Description    : Keythley 2450 Driver
# Git repository : 
# Copyright (C) TII (2023)
# ----------------------------------------------------------------------------

# -- include -----------------------------------------------------------------
import pyvisa
from pyvisa.errors import *
import numpy as np

from typing import Any, Callable, Dict, List, Optional, Union
from qcodes import validators as vals
from qcodes import VisaInstrument, InstrumentChannel, Parameter

# -- class -------------------------------------------------------------------
class Keithley2450(VisaInstrument):
    '''Driver for the Keythley 2450 Source Meter Unit'''
    
    def __init__(self,
                 name: str='',
                 address: str = '192.168.0.4',
                 terminator="\r",
                 identifier: Optional[str] = None,
                 **kwargs
                 ):
        super().__init__(name, address, terminator="\r", **kwargs)

        #self._name = name
        #self._address = address
        #if identifier is None:
        #    identifier = name
        #self._identifier = identifier
        #self._rm = pyvisa.ResourceManager()
        self._inst = pyvisa.Resource()

        # --- QCodes Parameters ----------------------------------------------
        self.add_parameter("Voltage",
            label="Source Voltage",
            docstring="Sets/gets SMU voltage output.",
            unit="V",
            vals=vals.Numbers(),
            set_parser=float,
            get_parser=float,
            set_cmd=self._set_voltage,
            get_cmd=self._get_voltage,
        )

        self.add_parameter(            "Current",
            label="Source Current",
            docstring="Sets/gets SMU Current output.",
            unit="A",
            vals=vals.Numbers(),
            set_parser=float,
            get_parser=float,
            set_cmd=self._set_voltage,
            get_cmd=self._get_voltage,
        )

        self.add_parameter(            "Compliance",
            label="Source Compliance",
            docstring="Sets/gets SMU Compliance output.",
            unit="",
            vals=vals.Numbers(0, 1e-3),
            set_parser=float,
            get_parser=float,
            set_cmd=self._set_voltage,  
            get_cmd=self._get_voltage,
        )

        self.connect_message()


    # ------------------------------------------------------------------------
    # --- Function Definitios ------------------------------------------------
    # ------------------------------------------------------------------------

    def reset(self) -> None:
        """
        Resets device, invalidates QCoDeS parameter cache and clears all
        status and event registers (see
        `SCPI <https://www.ivifoundation.org/docs/scpi-99.pdf>`_).

        Parameters
        ----------

        Returns
        ----------

        Raises
        ----------
        """

        # Reset
        self._reset()


    # ------------------------------------------------------------------------
    def __query__(self,cmd):
        """
        Sends a query to an instrument.

        Parameters
        ----------
        slot: int
            Slot index

        Returns
        ----------

        Raises
        ----------
        """
        rep = ''
        if cmd.find('?')==-1:
            self._inst.write(cmd)
        else:
            try:
                rep = self._inst.query(cmd)
            except VisaIOError as e:
                if(e.error_code != StatusCode.error_timeout):
                    try:
                        self._inst.write("*IDN?")
                    except VisaIOError as e2:
                        raise e2
                    finally:
                        self._inst.close()
                print("No reply from device. Command \'%s\' may be wrongly formatted or not supported."%cmd)
        return rep
        
    def open(self):
        '''Open conenction to instrument using address given in the constructor'''
        self._inst = self._rm.open_resource(self._address, timeout=1)
        self.idn = self.__query__("*IDN?")
        self.version = self.__query__(":SYST:VERSION?")
    
    def close(self):
        self._inst.close()


    def _set_voltage(self,voltage):
        return
    
    def _get_voltage(self):
        return
    

    def _set_current(self,voltage):
        return
    
    def _get_current(self):
        return
    

    def _set_compliance(self, comp):
        return
    
    def _get_compliance(self):
        return


    def _set_source_type(self,type):
        return
        
    def _get_source_type(self):
        return
    
     
    