from pathlib import Path
import csv
import logging
import re

# Controller function
def file_parser(master):
    try:
        if master.filename.suffix.lower().endswith('txt'):
            text_file_parser(master)
        elif master.filename.suffix.lower().endswith('arw'):
            arw_file_parser(master)
        elif master.filename.suffix.lower().endswith('csv'):
            csv_file_parser(master)
        else:
            logging.getLogger(__name__).warn('Unsopported filetype for '+str(master.filename))
    except Exception as e:
        logging.getLogger(__name__).error(e)

# Individual file format functions
def arw_file_parser(master):
    chrom_data = []
    with Path(master.filename).open() as fr:
        for line in fr:
            lines = line.split('\r')
            for line in lines:
                if line[0][0].isdigit() is False:
                    pass
                else:
                    line_chunks = line.rstrip()
                    line_chunks = line_chunks.split()
                    chrom_data.append((float(line_chunks[0]),
                                       float(line_chunks[1])))
    master.chrom_data = chrom_data

def csv_file_parser(master):
    chrom_data = []
    with Path(master.filename).open() as fr:
        csv_reader = csv.reader(fr)
        for row in csv_reader:
            chrom_data.append((float(row[0]), float(row[-1])))

    master.chrom_data = chrom_data

def text_file_parser(master):
    chrom_data = []
    with Path(master.filename).open() as fr:
        for line in fr:
            if line[0].isdigit() is True:
                line_chunks = line.strip().split()

                time_sep = re.sub(r'-?\d', '', line_chunks[0],
                                  flags=re.U)
                for sep in time_sep[:-1]:
                    line_chunks[0] = line_chunks[0].replace(sep, '')
                if time_sep:
                    line_chunks[0] = line_chunks[0].replace(
                        time_sep[-1], '.')

                int_sep = re.sub(r'-?\d', '', line_chunks[-1],
                                 flags=re.U)
                for sep in int_sep[:-1]:
                    line_chunks[-1] = line_chunks[-1].replace(
                        sep[-1], '')
                if int_sep:
                    line_chunks[-1] = line_chunks[-1].replace(
                        int_sep[-1], '.')

                chrom_data.append((float(line_chunks[0]),
                                   float(line_chunks[-1])))
    master.chrom_data = chrom_data
