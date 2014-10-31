#!/usr/bin/python

import time
import readline as rl
import threading
import serial
import cmd


def bias_to_12bit(value):
    try:
        value = float(value)
        val = int(4095. / 19 * (value - 15))
    except ValueError:
        raise ValueError('voltage must be a number')
    if 0 <= val < 4096:
        return val
    else:
        raise ValueError('invalid voltage\nmust be between 15 V and 34 V')


def offset_to_12bit(value):
    try:
        val = int(value)
    except ValueError:
        raise ValueError('voltage must be a number')
    if 0 <= val < 4096:
        return val
    else:
        raise ValueError('invalid voltage\nmust be between 15 V and 34 V')


def channel(value):
    try:
        val = int(value)
    except ValueError:
        raise ValueError('channel must be between a number between 0 and 15')
    if 0 <= val < 16:
        return val
    else:
        raise ValueError('invalid channel number\nmust be between 0 and 15')


class CameraCommunication(cmd.Cmd):
    """serial communication between adapter board and pc"""
    def __init__(self):
        cmd.Cmd.__init__(self)
        self.ser = serial.Serial()
        self.ser.port = '/dev/ttyUSB0'
        self.ser.baudrate = 115200
        self.ser.timeout = .05
        self.prompt = 'cam> '
        self.running = False

    def emptyline(self):
        pass

    def run(self):
            self.t1 = threading.Thread(target=self._run)
            self.t1.start()
            if self.running:
                self.cmdloop()

    def _run(self):
        try:
            self.ser.open()
            self.running = True
        except OSError as er:
            print er
        while self.running:
            time.sleep(.1)
            if self.ser.inWaiting() > 1:
                self.stdout.write('\r'+' '*(len(rl.get_line_buffer())+2)+'\r')
                self.handle_output(self.ser.readall())
                self.stdout.write(self.prompt + rl.get_line_buffer())
                self.stdout.flush()

    def shutdown(self):
        self.running = False
        try:
            self.t1.join()
            self.ser.close()
        except AttributeError:
            pass

    def send(self, s):
        self.ser.write(''.join((s, '\r')).encode())

    def handle_output(self, readall):
        print readall.strip()

    #command line methods
    def do_exit(self, line):
        self.shutdown()
        return True

    def help_exit(self):
        print 'exits the command line program.'

    def do_id(self, line):
        self.send('*idn?')

    def help_id(self, line):
        print 'displays the id-string of the camera.'

    def do_set_bias(self, args):
        try:
            chnl, num = args.strip().split()
            try:
                self.send(' '.join(('sb', str(channel(chnl)),
                          str(bias_to_12bit(num)))))
            except ValueError as er:
                print er
        except ValueError as er:
            print er

    def help_set_bias(self):
        print '\n'.join(('set bias voltage in one channel.',
                        'syntax:', 'set_bias #channel #voltage',
                        '0 <= #channel <= 15', '15 <= #voltage(volts) <= 34'))

    def do_set_bias_all(self, args):
        """set bias voltage for all channels.\nsyntax:\nset_bias_all #voltage
15 <= #voltage(volts) <= 34"""
        try:
            num = args.strip()
            try:
                self.send(' '.join(('sba', str(bias_to_12bit(num)))))
            except ValueError as er:
                print er
        except ValueError as er:
            print er

    def do_set_offset(self, args):
        try:
            chnl, num = args.strip().split()
            try:
                self.send(' '.join(('so', str(channel(chnl)),
                          str(offset_to_12bit(num)))))
            except ValueError as er:
                print er
        except ValueError as er:
            print er

    def help_set_offset(self):
        print '\n'.join(('set offset voltage in one channel.',
                        'syntax:', 'set_offset #channel #voltage',
                        '0 <= #channel <= 15', '? <= #voltage(volts) <= ?'))

    def do_set_offset_all(self, args):
        """set bias voltage for all channels.\nsyntax:\nset_offset_all #voltage
15 <= #voltage(volts) <= 34"""
        try:
            num = args.strip()
            try:
                self.send(' '.join(('soa', str(offset_to_12bit(num)))))
            except ValueError as er:
                print er
        except ValueError as er:
            print er

cam = CameraCommunication()
try:
    cam.run()
except KeyboardInterrupt:
    cam.shutdown()
