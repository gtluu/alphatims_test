import os
import sys
import logging
from timsconvert.timestamp import get_iso8601_timestamp, get_timestamp
from timsconvert.data_input import check_for_multiple_analysis, schema_detection
from timsconvert.classes import TimsconvertBafData, TimsconvertTsfData, TimsconvertTdfData
from timsconvert_nanodesi.write import write_nanodesi_imzml
from pyTDFSDK.init_tdf_sdk import init_tdf_sdk_api
from pyTDFSDK.ctypes_data_structures import PressureCompensationStrategy
from pyBaf2Sql.init_baf2sql import init_baf2sql_api


def convert_raw_file(tuple_args):
    run_args = tuple_args[0]
    infile = tuple_args[1]
    # Set output directory to default if not specified.
    if run_args['outdir'] == '':
        run_args['outdir'] = os.path.split(infile)[0]

    # Initialize logger if not running on server.
    logname = 'tmp_log_' + os.path.splitext(os.path.split(infile)[-1])[0] + '.log'
    if run_args['outdir'] == '' and os.path.isdir(run_args['input']) and os.path.splitext(run_args['input'])[
        -1] != '.d':
        logfile = os.path.join(run_args['input'], logname)
    elif run_args['outdir'] == '' and os.path.isdir(run_args['input']) and os.path.splitext(run_args['input'])[
        -1] == '.d':
        logfile = os.path.split(run_args['input'])[0]
        logfile = os.path.join(logfile, logname)
    else:
        logfile = os.path.join(run_args['outdir'], logname)
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=logfile, level=logging.INFO)
    if run_args['verbose']:
        logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    # Initialize Bruker DLL.
    logging.info(get_iso8601_timestamp() + ':' + 'Initialize Bruker .dll file...')
    tdf_sdk_dll = init_tdf_sdk_api()
    baf2sql_dll = init_baf2sql_api()

    # Read in input file (infile).
    logging.info(get_iso8601_timestamp() + ':' + 'Reading file: ' + infile)
    if not check_for_multiple_analysis(infile):
        schema = schema_detection(infile)
        if schema == 'TSF':
            if run_args['use_raw_calibration']:
                use_recalibrated_state = False
            elif not run_args['use_raw_calibration']:
                use_recalibrated_state = True
            data = TimsconvertTsfData(infile, tdf_sdk_dll, use_recalibrated_state=use_recalibrated_state)
        elif schema == 'TDF':
            if run_args['use_raw_calibration']:
                use_recalibrated_state = False
            elif not run_args['use_raw_calibration']:
                use_recalibrated_state = True
            if run_args['pressure_compensation_strategy'] == 'none':
                pressure_compensation_strategy = PressureCompensationStrategy.NoPressureCompensation
            elif run_args['pressure_compensation_strategy'] == 'global':
                pressure_compensation_strategy = PressureCompensationStrategy.AnalyisGlobalPressureCompensation
            elif run_args['pressure_compensation_strategy'] == 'frame':
                pressure_compensation_strategy = PressureCompensationStrategy.PerFramePressureCompensation
            data = TimsconvertTdfData(infile,
                                      tdf_sdk_dll,
                                      use_recalibrated_state=use_recalibrated_state,
                                      pressure_compensation_strategy=pressure_compensation_strategy)
        elif schema == 'BAF':
            if run_args['use_raw_calibration']:
                raw_calibration = True
            elif not run_args['use_raw_calibration']:
                raw_calibration = False
            data = TimsconvertBafData(infile, baf2sql_dll, raw_calibration=raw_calibration)
    else:
        return

    # Log arguments.
    for key, value in run_args.items():
        logging.info(get_iso8601_timestamp() + ':' + str(key) + ': ' + str(value))

    if schema == 'BAF':
        logging.info(get_iso8601_timestamp() + ':' + '.baf file detected...')
    elif schema == 'TSF':
        logging.info(get_iso8601_timestamp() + ':' + '.baf file detected...')
    elif schema == 'TDF':
        logging.info(get_iso8601_timestamp() + ':' + '.tdf file detected...')
    else:
        logging.warning(get_iso8601_timestamp() + ':' + 'Unable to determine acquisition mode using metadata for' +
                        infile + '...')
        logging.warning(get_iso8601_timestamp() + ':' + 'Exiting...')
    outfile = os.path.splitext(os.path.split(infile)[-1])[0] + '.imzML'
    logging.info(get_iso8601_timestamp() + ':' + 'Processing nano-DESI data...')
    write_nanodesi_imzml(data,
                         outdir=run_args['outdir'],
                         outfile=outfile,
                         mode=run_args['mode'],
                         exclude_mobility=run_args['exclude_mobility'],
                         profile_bins=run_args['profile_bins'],
                         imzml_mode=run_args['imzml_mode'],
                         mz_encoding=run_args['mz_encoding'],
                         intensity_encoding=run_args['intensity_encoding'],
                         mobility_encoding=run_args['mobility_encoding'],
                         compression=run_args['compression'],
                         chunk_size=10)

    logging.info('\n')

    for hand in logging.getLogger().handlers:
        logging.getLogger().removeHandler(hand)

    return logfile


def clean_up_logfiles(args, list_of_logfiles):
    # Concatenate log files.
    concat_logfile = ''
    for logfile in list_of_logfiles:
        with open(logfile, 'r') as logfile_obj:
            concat_logfile += logfile_obj.read()
    # Determine final log filename.
    logname = 'log_' + get_timestamp() + '.log'
    if args['outdir'] == '' and os.path.isdir(args['input']) and os.path.splitext(args['input'])[-1] != '.d':
        final_logfile = os.path.join(args['input'], logname)
    elif args['outdir'] == '' and os.path.isdir(args['input']) and os.path.splitext(args['input'])[-1] == '.d':
        final_logfile = os.path.split(args['input'])[0]
        final_logfile = os.path.join(final_logfile, logname)
    else:
        final_logfile = os.path.join(args['outdir'], logname)
    with open(final_logfile, 'w') as final_logfile_obj:
        final_logfile_obj.write(concat_logfile)
        print(get_iso8601_timestamp() + ':' + 'Final log file written to ' + final_logfile + '...')
    # Delete temporary log files.
    for logfile in list_of_logfiles:
        try:
            os.remove(logfile)
            print(get_iso8601_timestamp() + ':' + 'Removed temporary log file ' + logfile + '...')
        except OSError:
            print(get_iso8601_timestamp() + ':' + 'Unable to remove temporary log file ' + logfile + '...')
