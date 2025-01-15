import logging
from timsconvert.timestamp import *


logging.info(get_iso8601_timestamp() + ':' + 'Initialize constants...')

VERSION = '2.0.1'

# Only type "1" has been tested with the current test dataset.
INSTRUMENT_SOURCE_TYPE = {'0': 'unspecified',
                          '1': 'electrospray ionization',
                          '2': 'atmospheric pressure chemical ionization',
                          '3': 'nanoelectrospray',
                          '4': 'nanoelectrospray',
                          '5': 'atmospheric pressure photoionization',
                          '6': 'multimode ionization',
                          '9': 'nanoflow electrospray ionization',
                          '10': 'ionBooster',
                          '11': 'CaptiveSpray',
                          '12': 'GC-APCI',
                          '13': 'VIP-HESI-APCI',
                          '18': 'VIP-HESI'}

MSMS_TYPE = {'0': 'MS',
             '2': 'MS/MS',
             '8': 'dda-PASEF',
             '9': 'dia-PASEF'}

MSMS_TYPE_CATEGORY = {'ms1': [0],
                      'ms2': [2, 8, 9]}
