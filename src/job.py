
import os
import subprocess

from pathlib import Path
from typing import List, Union
from datetime import datetime


class Job:
    def __init__(self,
                 command: Union[str, List[str]],
                 working_directory: Union[str, Path, None],
                 sample=None):

        self.command = command
        self.working_directory = working_directory
        self.returncode = None
        self.time_start = str(datetime.now()).replace(' ', '')
        self.time_end = ''
        self.stdout: bytes = b''
        self.stderr: bytes = b''
        self.sample = sample

    def run(self):
        if self.working_directory is not None:
            os.chdir(self.working_directory)

        command = self.command.split(' ') if isinstance(self.command, str) else self.command
        p = subprocess.Popen(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        self.stdout, self.stderr = p.communicate()
        self.time_end = str(datetime.now()).replace(' ', '')
        self.returncode = p.returncode
