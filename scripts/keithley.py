# software for sourcing and measuring with the keithley 2400 over serial
# interface with baud 9600,
import serial
import time
import datetime as dt


test_string = 'KEITHLEY INSTRUMENTS INC.,MODEL 2400,1232220,C32'
delay = 0.2


def delayed(delay):
    def decorator(func):
        def inner(*args):
            ret = func(*args)
            time.sleep(delay)
            return ret
        return inner
    return decorator


class Keithley(object):
    """docstring for Keithley"""
    def __init__(self):
        super(Keithley, self).__init__()
        self._voltage = None
        self.time_offset = None
        # standard setup
        self.ser = serial.Serial()
        self.ser.baudrate = 9600
        self.ser.port = '/dev/ttyUSB0'
        self.ser.timeout = .1
        self.ser.open()
        # check connection
        self.to_keithley('*idn?')
        if test_string not in self.from_keithley():
            raise Exception('something went wrong')
        # reset
        self.to_keithley('*rst')
        self.start_config()
        print 'Keithley 2400 is ready'

    @property
    def voltage(self):
        return self._voltage

    @voltage.setter
    def voltage(self, volt):
        self.change_voltage(volt)

    @delayed(delay)
    def to_keithley(self, string):
        self.ser.write(''.join((string, '\r')).encode())

    def from_keithley(self):
        return self.ser.readall()

    def start_config(self):
        # disable beep
        self.to_keithley(':syst:beep:stat 0')
        self.to_keithley(':sour:func volt')
        self.to_keithley(':sour:volt:mode fixed')
        self.to_keithley(':sour:volt:rang 30')
        # set compliace
        self.to_keithley(':sens:curr:prot .1')
        self.to_keithley(':sens:func "curr:dc"')
        self.to_keithley(':sens:curr:rang .1')
        self.to_keithley(':form:elem:sens curr, volt, time')
        self.to_keithley(':syst:time:res')
        self.time_offset = dt.datetime.now() - dt.timedelta(seconds=delay)
        self.to_keithley(':form:elem:sens?')
        print 'configured with data format {}'.format(self.from_keithley())

    def source_on(self):
        self.to_keithley(':outp:stat 1')
        'turned on with voltage {}'.format(self.voltage)

    def source_off(self):
        self.to_keithley(':outp:stat 0')

    def change_voltage(self, volt):
        self.to_keithley(':sour:volt:lev {}'.format(float(volt)))
        self._voltage = float(volt)

    def measure(self):
        self.to_keithley(':meas?')
        volt, curr, time = map(float, self.from_keithley().strip().split(','))
        return volt, curr, time
