import serial
import time
from datetime import datetime
from gevent import socket, sleep


test_control = 'IAAT_ClimateChamber'
test_readout = 'IAAT_PT100_Sensors'


class SerialStandard(object):
    """docstring for SerialStandard"""
    def __init__(self, tty_id=0):
        super(SerialStandard, self).__init__()
        assert isinstance(tty_id, int)
        # setup serial port
        self.ser = serial.Serial()
        self.ser.port = '/dev/ttyUSB{}'.format(tty_id)
        self.ser.baudrate = 115200
        self.ser.timeout = .1
        self.ser.open()

    def write(self, st):
        self.ser.write(''.join((st, '\r\n')).encode())

    def read(self):
        return self.ser.readall()


class TempReadout(SerialStandard):
    """docstring for temp_readout"""
    def __init__(self, tty_id=0):
        super(TempReadout, self).__init__(tty_id)
        # check connection
        self.write('*idn?')
        answer = self.read()
        if test_readout not in answer:
            raise Exception('could not verify readout')
        print answer.strip()

    def read_temperature(self, channel=7):
        # flush buffer
        try:
            self.ser.readall()
            self.write('tempf? 4')
            offset = float(self.read())
            self.write('tempf? {}'.format(channel))
            temp = float(self.read()) - offset
            return temp
        except serial.SerialException as err:
            time.sleep(.2)
            return self.read_temperature(channel)
        except ValueError as err:
            time.sleep(1)
            return self.read_temperature(channel)
        except OSError as err:
            time.sleep(2)
            return self.read_temperature(channel)


def logged(func):
    def inner(*args):
        ret = func(*args)
        temp = args[-1]
        with open('log.txt', 'a') as f:
            f.write(','.join(map(str, (datetime.now(), temp))) + '\n')
        return ret
    return inner


class TempControl(SerialStandard):
    """docstring for TempControl"""
    def __init__(self, tty_id=0):
        super(TempControl, self).__init__(tty_id)
        # check connection
        self.write('*idn?')
        answer = self.read()
        if test_control not in answer:
            raise Exception('could not verify control')
        print answer.strip()
        self.temperature = None

    @logged
    def set_temperature(self, temp):
        assert isinstance(temp, (int, float))
        assert -30 < temp <= 50
        self.write('set {}'.format(temp))
        self.temperature = temp
        self.read()


def read_temp(channel):
    assert isinstance(channel, int)
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.connect(('localhost', 19995))
    f = sock.makefile()
    f.write('get_temp {}\n'.format(channel))
    f.flush()
    temp = f.readline().strip()
    sock.close()
    return float(temp)

