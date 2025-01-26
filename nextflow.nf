#!~/bin nextflow

nextflow.enable.dsl=2

// params

// required params
params.input = ''

// optional params
params.mode = 'centroid'  // mode can be 'centroid' or 'raw'
params.compression = 'zlib' // zlib or none

// timsconvert params
params.ms2_only = 'False'  // only convert ms2 spectra?
params.use_raw_calibration = 'False'  // whether or not to use raw or recalibrated data (if available)
params.pressure_compensation_strategy = 'global'  // choose whether global, per frame, or no TIMS pressure compensation is used
params.exclude_mobility = 'True'  // exclude mobility arrays from MS1 spectra?
params.mz_encoding = 64  // 64 or 32 bit encoding
params.intensity_encoding = 64  // 64 or 32 bit encoding
params.mobility_encoding = 64  // 64 or 32 bit encoding
params.barebones_metadata = 'False'  // only use barebones metadata if downstream tools are not compatible with timstof cv params
params.profile_bins = 0  // perform binning of profile mode data into N bins; set to 0 to disable
params.maldi_output_mode = 'combined' // choose whether MALDI spectra are output to individual files, a single combined file with multiple spectra, or grouped by sample via maldi_plate_map
params.maldi_plate_map = ''  // path to plate map containing spot metadata in CSV format
params.imzml_mode = 'processed'  // whether to export imzml data as processed or continuous m/z arrays
params.iprm_format = 'mzml'  // export iprm-PASEF ms2 data as mzml, mgf, or imzml files
params.iprm_output_mode = 'combined'  // export iprm-PASEF data to combined or individual files per precursor; imzml will always be individual

// timsconvert system params
params.verbose = 'True'

// Boiler Plate
TOOL_FOLDER = "$baseDir/bin"
params.publishdir = "nf_output"

// Process
process convert {
    publishDir "$params.publishdir", mode: 'copy'

    input:
    file input_file

    output:
    file "spectra/*"

    script:
    def ms2_flag = params.ms2_only == 'True' ? "--ms2_only" : ''
    def use_raw_calibration_flag = params.use_raw_calibration == 'True' ? "--use_raw_calibration" : ''
    def exclude_mobility_flag = params.exclude_mobility == 'True' ? "--exclude_mobility" : ''
    def barebones_metadata_flag = params.barebones_metadata == 'True' ? "--barebones_metadata" : ''
    def verbose_flag = params.verbose == 'True' ? "--verbose" : ''

    """
    mkdir spectra
    python3 $TOOL_FOLDER/cmd.py \
    --input $input_file \
    --outdir spectra \
    --mode ${params.mode} \
    --compression ${params.compression} \
    ${ms2_flag} \
    ${use_raw_calibration_flag} \
    --pressure_compensation_strategy ${params.pressure_compensation_strategy} \
    ${exclude_mobility_flag} \
    --mz_encoding ${params.mz_encoding} \
    --intensity_encoding ${params.intensity_encoding} \
    --mobility_encoding ${params.mobility_encoding} \
    ${barebones_metadata_flag} \
    --profile_bins ${params.profile_bins} \
    --maldi_output_mode ${params.maldi_output_mode} \
    --maldi_plate_map ${params.maldi_plate_map} \
    --imzml_mode ${params.imzml_mode} \
    --iprm_format ${params.iprm_format} \
    --iprm_output_mode ${params.iprm_output_mode} \
    ${verbose_flag}
    """
}

workflow {
    input_ch = Channel.fromPath(params.input, type:'dir', checkIfExists: true)
    convert(input_ch)
}