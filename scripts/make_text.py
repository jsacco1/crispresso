#!/bin/env python
# Title: make_text.py
# Author: James Sacco

import argparse
import logging
import os
import sys

import pandas as pd
from pathlib import Path

logging.basicConfig(filename='tmp/make_text_log.log', encoding='utf-8', level=logging.DEBUG, filemode='w', format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
logging.debug('DEBUG message.')
logging.info('INFO message.')
logging.warning('WARNING message.')
logging.error('ERROR message.')

# disable pandas warnings
pd.options.mode.chained_assignment = None 

def main():

    logging.info('make_text.py started in %s', os.getcwd())

    parser = argparse.ArgumentParser(description="Make input.txt file for CRISPResso2 run")  #formatter_class=argparse.RawTextHelpFormatter
    parser.add_argument('input_file', type = str,
                        help='''input_file: Path to Excel (.xlsx) file to convert to txt file.\nExample: ~/path/to/input_excel.xlsx\n\n''',
                        metavar='')
    parser.add_argument('--quantification_window_center', type = int, default = -3, dest='quantification_window_center',
                        help="Center of quantification window to use within respect to the 3' end of the provided sgRNA sequence.\nThe default is -3 and is suitable for the Cas9 system.",
                        metavar='')
    parser.add_argument('--quantification_window_size', type = int, default = 1, dest='quantification_window_size',
                        help='''Defines the size (in bp) of the quantification window extending from the position specified by the "--cleavage_offset" or "--quantification_window_center."\ndefault: 1''',
                        metavar='')
    parser.add_argument('--min_average_read_quality', type = int, default = 10, dest='min_average_read_quality',
                        help="Minimum average quality score (phred33).\ndefault: 10",
                        metavar='')
    parser.add_argument('fastq_dir', type = str,
                        help = '''fastq_dir: Path to fastq.gz files.\nExample: ~/path/to/Fastq/\n\n''',
                        metavar='')
    parser.add_argument('--min_frequency_alleles_around_cut_to_plot', type = float, default = 0.2, dest='min_frequency_alleles_around_cut_to_plot', 
                        help="Minimum percent reads required to report an allele in the alleles table plot.\ndefault: 0.2",
                        metavar='')
    parser.add_argument('--plot_window_size', type = int, default = 20, dest='plot_window_size',
                        help="Defines the size of the window extending from the quantification window center to plot. Nucleotides within plot_window_size of the quantification_window_center for each guide are plotted.\ndefault: 25",
                        metavar='')
    parser.add_argument('outfile', type = str, default = "input.txt",
                        help='''output_file: Path to output file.\nDefault: input.txt\nExample: ~/path/to/input.txt \n''', 
                        metavar='')
    args = parser.parse_args()
    print("CRISPResso2 parameters: \n", args)
    logging.info("params:\n %s", args)

    make_text_input(
                    args.input_file, 
                    args.quantification_window_center, 
                    args.quantification_window_size,
                    args.min_average_read_quality, 
                    args.fastq_dir, 
                    args.min_frequency_alleles_around_cut_to_plot, 
                    args.plot_window_size, 
                    args.outfile
                    )


def make_text_input(
        input_file, 
        quantification_window_center, 
        quantification_window_size, 
        min_average_read_quality,
        fastq_dir,
        min_frequency_alleles_around_cut_to_plot, 
        plot_window_size, 
        outfile
        ):

    # Create input.txt (tsv) for CRISPRessoBatch from input_excel.xlsx.
    logging.info('Making input.txt file')

    try:
        input_df = pd.read_excel(input_file)
        logging.info('Successfully opened %s', input_file)
    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception:
        logging.info.error('Failed to open file', exc_info=True)

    logging.info('Checking for null values in input_excel.xlsx...')
    assert(len(input_df[input_df.isnull().any(axis=1)]) == 0)

    out_df = input_df[["NGS_sample_name"]]

    out_df.columns = ['name']
  
    out_df['fastq_r1'] = out_df['name'] + '_R1_001.fastq.gz'
    out_df['fastq_r1'] = fastq_dir + out_df['fastq_r1'].astype(str)

    out_df['fastq_r2'] = out_df['name'] + '_R2_001.fastq.gz'
    out_df['fastq_r2'] = fastq_dir + out_df['fastq_r2'].astype(str)

    logging.info('Checking that all paths to fastq files exist...')

    missing_fastq_files = {x for x in out_df['fastq_r2'] if not Path(x).exists()}.union({x for x in out_df['fastq_r1'] if not Path(x).exists()})

    for file in missing_fastq_files:
        try:
            f = open(file, 'r')
        except IOError:
            logging.warning("File does not exist: {}".format(file))

    out_df.insert(1, 'quantification_window_center', quantification_window_center)

    out_df.insert(2, 'quantification_window_size', quantification_window_size)

    out_df.insert(3, 'min_average_read_quality', min_average_read_quality)

    fastq_dir = Path(fastq_dir)
    assert fastq_dir.exists()

    out_df['amplicon_seq'] = input_df[["Ref_sequence"]]
    out_df['guide_seq'] = input_df[["Guide_Sequence"]]
    out_df['min_frequency_alleles_around_cut_to_plot'] = min_frequency_alleles_around_cut_to_plot
    out_df['plot_window_size'] = plot_window_size

    assert(len(out_df[out_df.isnull().any(axis=1)]) == 0)

    try:
        out_df.to_csv(outfile, sep="\t", index = False)
    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception as exception:
        logging.error('Failed to write file', exc_info=True)
    
    logging.info("Successfully converted %s to %s for CRISPRessoBatch!", input_file, outfile)

    sys.exit(0)

if __name__ == '__main__':
    main()