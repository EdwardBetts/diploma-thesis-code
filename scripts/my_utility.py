import datetime as dt
from sys import stdout


def build_time_string():
    now = dt.datetime.now()
    return '{:04}-{:02}-{:02}_{:02}:{:02}'.format(now.year,
                                                  now.month, now.day, now.hour, now.minute)


def print_progress(frac):
    length = 30
    block = int(round(length*frac))
    text = '\r[{}] {:.3g}%'.format('#'*block + '-'*(length - block), frac*100)
    stdout.write(text)
    stdout.flush()
