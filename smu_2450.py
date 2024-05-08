from qibolab.instruments.abstract import Instrument, InstrumentException
from qibo.config import log

from K2450_driver import K2450_driver as K2450

class smu_2450(Instrument):
    """A class to control a Keithley SMU 2450.

    Attributes:
        name (str): A unique name given to the instrument.
        address (str): IP_address:module_number (the IP address of the cluster and the 
            module number)
        device (K2450_driver): A reference to the underlying
            `K2450_driver.K2450_driver` object. It can be used to access other
            features not directly exposed by this wrapper.

    The class inherits from :class:`qibolab.instruments.abstract.Instrument` and implements its interface methods:
        __init__()
        connect()
        setup()
        start()
        stop()
        disconnect()
    """

    property_wrapper = lambda parent, *parameter: property(
        lambda self: parent.device.get(parameter[0]),
        lambda self, x: parent._set_device_parameter(parent.device, *parameter, value=x),
    )
    property_wrapper.__doc__ = """A lambda function used to create properties that wrap around the device parameters and
    caches their value using `_set_device_parameter()`.
    """

    def __init__(self, name: str, address: str):
        """Initialises the instrument storing its name, address and settings."""
        super().__init__(name, address)
        self.device: K2450 = None
        """Reference to the underlying `K2450` object."""

    def _set_device_parameter(self, target, *parameters, value):
        """Sets a parameter of the instrument, if it changed from the last stored in the cache.

        Args:
            target = an instance of K2450.K2450 or
            *parameters (list): A list of parameters to be cached and set.
            value = The value to set the paramters.
        Raises:
            Exception = If attempting to set a parameter without a connection to the instrument.
        """
        if self.is_connected:
            key = target.name + "." + parameters[0]
            if not key in self._device_parameters:
                for parameter in parameters:
                    if not hasattr(target, parameter):
                        raise Exception(f"The instrument {self.name} does not have parameters {parameter}")
                    target.set(parameter, value)
                self._device_parameters[key] = value
            elif self._device_parameters[key] != value:
                for parameter in parameters:
                    target.set(parameter, value)
                self._device_parameters[key] = value
        else:
            raise Exception("There is no connection to the instrument {self.name}")

    def _erase_device_parameters_cache(self):
        """Erases the cache of instrument parameters."""
        self._device_parameters = {}



    def connect(self):
        """Connects to the SMU.

        If the connection is successful, it resets the instrument and configures it with the stored settings.
        A reference to the underlying object is saved in the attribute `device`.
        """
        if not self.is_connected:
            for attempt in range(3):
                try:
                    self.device = K2450(self.name, self.address)
                    self.device.reset()
                    self.is_connected = True
                    break
                except Exception as exc:
                    log.info(f"Unable to connect:\n{str(exc)}\nRetrying...")
                # TODO: if not able to connect after 3 attempts, check for ping response and reboot
            if not self.is_connected:
                raise InstrumentException(self, f"Unable to connect to {self.name}")

        # apply stored settings
        self._setup()

    def _setup(self):
        pass

    def setup(self):
        """Configures the instrument with the stored settings."""
        self._setup()

    def start(self):
        """Empty method to comply with Instrument interface."""

    def stop(self):
        """Empty method to comply with Instrument interface."""

    def disconnect(self):
        """Closes the connection to the instrument."""
        self.is_connected = False
        self.device.close()