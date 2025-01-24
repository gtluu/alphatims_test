import os
import logging
from timsconvert.timestamp import get_iso8601_timestamp


def dot_d_detection(input_directory):
    """
    Search the input directory and any subdirectories for .d directory paths.

    :param input_directory: Path to the directory to be searched.
    :type: str
    :return: List of absolute paths for each .d directory found.
    :rtype: list
    """
    if not os.path.exists(input_directory):
        logging.warning(get_iso8601_timestamp() + ':' + f'{input_directory} does not exist...')
        logging.warning(get_iso8601_timestamp() + ':' + 'Skipping current file conversion...')
    else:
        return [os.path.join(dirpath, directory) for dirpath, dirnames, filenames in os.walk(input_directory)
                for directory in dirnames if directory.endswith('.d')]


def schema_detection(bruker_dot_d_file):
    """
    Detect the schema used by the raw data in the Bruker .d directory.

    :param bruker_dot_d_file: Path to the .d directory of interest.
    :type: str
    :return: Capitalized schema extension (TDF, TSF, or BAF).
    :rtype: str
    """
    exts = [os.path.splitext(fname)[1] for dirpath, dirnames, filenames in os.walk(bruker_dot_d_file)
            for fname in filenames]
    if '.tdf' in exts and '.tsf' not in exts and '.baf' not in exts:
        return 'TDF'
    elif '.tsf' in exts and '.tdf' not in exts and '.baf' not in exts:
        return 'TSF'
    elif '.baf' in exts and '.tdf' not in exts and '.tsf' not in exts:
        return 'BAF'


def check_for_multiple_analysis(bruker_dot_d_file):
    """
    Check to ensure that only a single .baf/.tsf/.tdf and associated .tsf_bin/.tdf_bin exists within the .d directory.

    :param bruker_dot_d_file: Path to the .d directory of interest.
    :type: str
    """
    fnames = [fname for dirpath, dirnames, filenames in os.walk(bruker_dot_d_file) for fname in filenames
              if fname in ['analysis.baf', 'analysis.tsf', 'analysis.tdf', 'analysis.tsf_bin', 'analysis.tdf_bin']]
    if len(fnames) != len(set(fnames)):
        logging.error(get_iso8601_timestamp() + ':' + 'Duplicate analysis file detected within .d directory...')
        logging.error(get_iso8601_timestamp() + ':' + 'Check the *.d directory for duplicate analysis.sqlite, '
                                                      'analysis.baf, analysis.tsf, analysis.tsf_bin, analysis.tdf, '
                                                      'or analysis.tdf_bin files...')
        logging.error(get_iso8601_timestamp() + ':' + 'Skipping conversion of ' + bruker_dot_d_file + '...')
        return True
    if 'analysis.tdf' in fnames and 'analysis.tsf' in fnames:
        logging.error(get_iso8601_timestamp() + ':' + 'Duplicate analysis file detected within .d directory...')
        logging.error(get_iso8601_timestamp() + ':' + 'Check the *.d directory for duplicate analysis.sqlite, '
                                                      'analysis.baf, analysis.tsf, analysis.tsf_bin, analysis.tdf, '
                                                      'or analysis.tdf_bin files...')
        logging.error(get_iso8601_timestamp() + ':' + 'Skipping conversion of ' + bruker_dot_d_file + '...')
        return True
    if 'analysis.tdf' in fnames and 'analysis.baf' in fnames:
        logging.error(get_iso8601_timestamp() + ':' + 'Duplicate analysis file detected within .d directory...')
        logging.error(get_iso8601_timestamp() + ':' + 'Check the *.d directory for duplicate analysis.sqlite, '
                                                      'analysis.baf, analysis.tsf, analysis.tsf_bin, analysis.tdf, '
                                                      'or analysis.tdf_bin files...')
        logging.error(get_iso8601_timestamp() + ':' + 'Skipping conversion of ' + bruker_dot_d_file + '...')
        return True
    if 'analysis.tsf' in fnames and 'analysis.baf' in fnames:
        logging.error(get_iso8601_timestamp() + ':' + 'Duplicate analysis file detected within .d directory...')
        logging.error(get_iso8601_timestamp() + ':' + 'Check the *.d directory for duplicate analysis.sqlite, '
                                                      'analysis.baf, analysis.tsf, analysis.tsf_bin, analysis.tdf, '
                                                      'or analysis.tdf_bin files...')
        logging.error(get_iso8601_timestamp() + ':' + 'Skipping conversion of ' + bruker_dot_d_file + '...')
        return True
    return False
