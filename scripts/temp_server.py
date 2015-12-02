from gevent.server import StreamServer
from gevent.event import Event
import gevent
from time import time
from temperature_readout import TempReadout


def get_temp(channel):
    pass


def handle(sock, address):
    fp = sock.makefile()
    while True:
        line = fp.readline().strip()
        if line.startswith('get_temp'):
            channel = int(line.split()[1])
            if channel in (6, 7):
                fp.write(str(get_temp.ts[channel]) + '\n')
                fp.flush()
        else:
            break
    fp.close()
            

def measure_temperature(ev):
    tr = TempReadout(1)
    t = time()
    while not ev.is_set():
        if time() - t > 2:
            get_temp.ts[6] = tr.read_temperature(6)
            get_temp.ts[7] = tr.read_temperature(7)
            t = time()
        gevent.sleep(.1)
            
            


if __name__ == '__main__':
    get_temp.ts = {6: 25, 7: 25}
    
    ev = Event()
    gevent.spawn(measure_temperature, ev)
    server = StreamServer(('localhost', 19995), handle)
    server.serve_forever() 
