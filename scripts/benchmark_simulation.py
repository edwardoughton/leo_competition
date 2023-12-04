"""
Script to solve for equilibrium and optima of the
orbital location game under benchmark calibration.

"""
import os
import configparser
import numpy as np
import pandas as pd 

from misc import shell_volume, velocity
from calibration import load_params_element

CONFIG = configparser.ConfigParser()
CONFIG.read(os.path.join(os.path.dirname(__file__), 'script_config.ini'))
BASE_PATH = CONFIG['file_locations']['base_path']

DATA_RAW = os.path.join(BASE_PATH, 'raw')
DATA_PROCESSED = os.path.join(BASE_PATH, 'processed')




if __name__ == '__main__':

    params = load_params_element()

    params['atmo_damages'] = 0

