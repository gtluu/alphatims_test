import os
import sys
import logging
from timsconvert.timestamp import get_iso8601_timestamp, get_timestamp
from timsconvert.data_input import check_for_multiple_analysis, schema_detection
from timsconvert.classes import TimsconvertBafData, TimsconvertTsfData, TimsconvertTdfData
from timsconvert.write import (write_lcms_mzml, write_maldi_dd_mzml, write_maldi_ims_imzml, write_maldi_ims_iprm_imzml,
                               write_maldi_ims_iprm_mzml, write_maldi_ims_iprm_mgf)
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

    # BAF ESI-MS Dataset
    if schema == 'BAF':
        logging.info(get_iso8601_timestamp() + ':' + '.baf file detected...')
        outfile = os.path.splitext(os.path.split(infile)[-1])[0] + '.mzML'
        logging.info(get_iso8601_timestamp() + ':' + 'Processing LC-MS data...')
        write_lcms_mzml(data,
                        infile,
                        run_args['outdir'],
                        outfile,
                        run_args['mode'],
                        run_args['ms2_only'],
                        run_args['exclude_mobility'],
                        run_args['profile_bins'],
                        run_args['mz_encoding'],
                        run_args['intensity_encoding'],
                        run_args['mobility_encoding'],
                        run_args['compression'],
                        run_args['barebones_metadata'])

    # TSF ESI-MS Dataset
    elif schema == 'TSF' and 'MaldiApplicationType' not in data.analysis['GlobalMetadata'].keys():
        logging.info(get_iso8601_timestamp() + ':' + '.tsf file detected...')
        outfile = os.path.splitext(os.path.split(infile)[-1])[0] + '.mzML'
        logging.info(get_iso8601_timestamp() + ':' + 'Processing LC-MS data...')
        write_lcms_mzml(data,
                        infile,
                        run_args['outdir'],
                        outfile,
                        run_args['mode'],
                        run_args['ms2_only'],
                        run_args['exclude_mobility'],
                        run_args['profile_bins'],
                        run_args['mz_encoding'],
                        run_args['intensity_encoding'],
                        run_args['mobility_encoding'],
                        run_args['compression'],
                        run_args['barebones_metadata'])

    # TSF MALDI-qTOF Dried Droplet Dataset
    elif schema == 'TSF' \
            and 'MaldiApplicationType' in data.analysis['GlobalMetadata'].keys() \
            and data.analysis['GlobalMetadata']['MaldiApplicationType'] == 'SingleSpectra':
        logging.info(get_iso8601_timestamp() + ':' + '.tsf file detected...')
        outfile = os.path.splitext(os.path.split(infile)[-1])[0] + '.mzML'
        logging.info(get_iso8601_timestamp() + ':' + 'Processing MALDI dried droplet data...')
        write_maldi_dd_mzml(data,
                            infile,
                            run_args['outdir'],
                            outfile,
                            run_args['mode'],
                            run_args['ms2_only'],
                            run_args['exclude_mobility'],
                            run_args['profile_bins'],
                            run_args['mz_encoding'],
                            run_args['intensity_encoding'],
                            run_args['mobility_encoding'],
                            run_args['compression'],
                            run_args['maldi_output_mode'],
                            run_args['maldi_plate_map'],
                            run_args['barebones_metadata'])

    # TSF MALDI-qTOF MSI Dataset
    elif schema == 'TSF' \
            and 'MaldiApplicationType' in data.analysis['GlobalMetadata'].keys() \
            and data.analysis['GlobalMetadata']['MaldiApplicationType'] == 'Imaging':
        logging.info(get_iso8601_timestamp() + ':' + '.tsf file detected...')
        outfile = os.path.splitext(os.path.split(infile)[-1])[0] + '.imzML'
        logging.info(get_iso8601_timestamp() + ':' + 'Processing MALDI imaging mass spectrometry data...')
        write_maldi_ims_imzml(data,
                              run_args['outdir'],
                              outfile,
                              run_args['mode'],
                              run_args['exclude_mobility'],
                              run_args['profile_bins'],
                              run_args['imzml_mode'],
                              run_args['mz_encoding'],
                              run_args['intensity_encoding'],
                              run_args['mobility_encoding'],
                              run_args['compression'])

    # TDF ESI-TIMS-MS Dataset
    elif schema == 'TDF' \
            and 'MaldiApplicationType' not in data.analysis['GlobalMetadata'].keys():
        logging.info(get_iso8601_timestamp() + ':' + '.tdf file detected...')
        outfile = os.path.splitext(os.path.split(infile)[-1])[0] + '.mzML'
        logging.info(get_iso8601_timestamp() + ':' + 'Processing LC-TIMS-MS data...')
        write_lcms_mzml(data,
                        infile,
                        run_args['outdir'],
                        outfile,
                        run_args['mode'],
                        run_args['ms2_only'],
                        run_args['exclude_mobility'],
                        run_args['profile_bins'],
                        run_args['mz_encoding'],
                        run_args['intensity_encoding'],
                        run_args['mobility_encoding'],
                        run_args['compression'],
                        run_args['barebones_metadata'])

    # TDF MALDI-TIMS-qTOF Dried Droplet Dataset
    elif schema == 'TDF' \
            and 'MaldiApplicationType' in data.analysis['GlobalMetadata'].keys() \
            and data.analysis['GlobalMetadata']['MaldiApplicationType'] == 'SingleSpectra':
        logging.info(get_iso8601_timestamp() + ':' + '.tdf file detected...')
        outfile = os.path.splitext(os.path.split(infile)[-1])[0] + '.mzML'
        logging.info(get_iso8601_timestamp() + ':' + 'Processing MALDI-TIMS dried droplet data...')
        write_maldi_dd_mzml(data,
                            infile,
                            run_args['outdir'],
                            outfile,
                            run_args['mode'],
                            run_args['ms2_only'],
                            run_args['exclude_mobility'],
                            run_args['profile_bins'],
                            run_args['mz_encoding'],
                            run_args['intensity_encoding'],
                            run_args['mobility_encoding'],
                            run_args['compression'],
                            run_args['maldi_output_mode'],
                            run_args['maldi_plate_map'],
                            run_args['barebones_metadata'])

    # TDF MALDI-TIMS-qTOF MSI Dataset
    elif schema == 'TDF' \
            and 'MaldiApplicationType' in data.analysis['GlobalMetadata'].keys() \
            and data.analysis['GlobalMetadata']['MaldiApplicationType'] == 'Imaging':
        # Check if iprm-PASEF dataset.
        msms_mode_id = data.analysis['PropertyDefinitions'][data.analysis['PropertyDefinitions']['PermanentName'] ==
                                                            'Mode_ScanMode'].to_dict(orient='records')[0]['Id']
        properties_dicts = data.analysis['Properties'][data.analysis['Properties']['Property'] ==
                                                       msms_mode_id].to_dict(orient='records')
        # iprm-PASEF workflow writes different MGF, mzML, or imzML/ibd files for each precursor.
        if all(i['Value'] == 12 for i in properties_dicts):
            logging.info(get_iso8601_timestamp() + ':' + '.tdf file detected...')
            logging.info(
                get_iso8601_timestamp() + ':' + 'Processing MALDI-TIMS iprm-PASEF imaging mass spectrometry data...')
            if run_args['iprm_format'] == 'mgf':
                logging.info(get_iso8601_timestamp() + ':' + 'Writing data to MGF file(s)...')
                outfile = os.path.splitext(os.path.split(infile)[-1])[0] + '.mgf'
                write_maldi_ims_iprm_mgf(data,
                                         run_args['outdir'],
                                         outfile,
                                         run_args['mode'],
                                         run_args['exclude_mobility'],
                                         run_args['profile_bins'],
                                         run_args['mz_encoding'],
                                         run_args['intensity_encoding'],
                                         run_args['mobility_encoding'],
                                         run_args['iprm_output_mode'])
            elif run_args['iprm_format'] == 'mzml':
                logging.info(get_iso8601_timestamp() + ':' + 'Writing data to mzML file(s)...')
                outfile = os.path.splitext(os.path.split(infile)[-1])[0] + '.mzML'
                write_maldi_ims_iprm_mzml(data,
                                          infile,
                                          run_args['outdir'],
                                          outfile,
                                          run_args['mode'],
                                          run_args['ms2_only'],
                                          run_args['exclude_mobility'],
                                          run_args['profile_bins'],
                                          run_args['mz_encoding'],
                                          run_args['intensity_encoding'],
                                          run_args['mobility_encoding'],
                                          run_args['compression'],
                                          run_args['barebones_metadata'],
                                          run_args['iprm_output_mode'])
            elif run_args['iprm_format'] == 'imzml':
                logging.info(get_iso8601_timestamp() + ':' + 'Writing data to imzML file(s)...')
                outfile = os.path.splitext(os.path.split(infile)[-1])[0] + '.imzML'
                write_maldi_ims_iprm_imzml(data,
                                           run_args['outdir'],
                                           outfile,
                                           run_args['mode'],
                                           run_args['exclude_mobility'],
                                           run_args['profile_bins'],
                                           run_args['imzml_mode'],
                                           run_args['mz_encoding'],
                                           run_args['intensity_encoding'],
                                           run_args['mobility_encoding'],
                                           run_args['compression'])
        # Standard imzML export workflow.
        else:
            logging.info(get_iso8601_timestamp() + ':' + '.tdf file detected...')
            outfile = os.path.splitext(os.path.split(infile)[-1])[0] + '.imzML'
            logging.info(get_iso8601_timestamp() + ':' + 'Processing MALDI-TIMS imaging mass spectrometry data...')
            write_maldi_ims_imzml(data,
                                  run_args['outdir'],
                                  outfile,
                                  run_args['mode'],
                                  run_args['exclude_mobility'],
                                  run_args['profile_bins'],
                                  run_args['imzml_mode'],
                                  run_args['mz_encoding'],
                                  run_args['intensity_encoding'],
                                  run_args['mobility_encoding'],
                                  run_args['compression'])

    else:
        logging.warning(get_iso8601_timestamp() + ':' + 'Unable to determine acquisition mode using metadata for' +
                        infile + '...')
        logging.warning(get_iso8601_timestamp() + ':' + 'Exiting...')

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
