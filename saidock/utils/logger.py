import time
from pathlib import Path

LEVELS = {'DEBUG':0, 'INFO':1, 'OK':2, 'WARN':3, 'ERROR':4, 'STEP':5}

class SAILogger:
    def __init__(self, log_path=None, verbose=True):
        self.log_path = Path(log_path) if log_path else None
        self.verbose  = verbose
        self.entries  = []

    def log(self, msg, level='INFO'):
        ts    = time.strftime('%H:%M:%S')
        entry = f'[{ts}][{level}] {msg}'
        self.entries.append({'time': ts, 'level': level, 'msg': msg})
        if self.verbose:
            print(f'  {entry}')
        if self.log_path:
            with open(self.log_path, 'a') as f:
                f.write(entry + '\n')

    def step(self, msg):   self.log(msg, 'STEP')
    def ok(self, msg):     self.log(msg, 'OK')
    def warn(self, msg):   self.log(msg, 'WARN')
    def error(self, msg):  self.log(msg, 'ERROR')
    def info(self, msg):   self.log(msg, 'INFO')
