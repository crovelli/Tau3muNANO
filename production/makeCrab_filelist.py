import os
import glob
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-d', '--input_dir', type = str, default = 'MCsamples_2022.yml', help = 'Directory storing crab input files')
parser.add_argument('-r', '--remove', default='/eos/cms', help = 'leading character to remove')
parser.add_argument('-o', '--output', default='filelist_crab', help = 'file name in output')
args = parser.parse_args()

file_list  = glob.glob(args.input_dir+"/*.root")
formatted_files = [f.replace(args.remove, '', 1) for f in file_list]
input_to_crab = args.output + '.txt'
with open(input_to_crab, 'w') as file:
    [file.write('%s\n'%f) for f in formatted_files]

#print(formatted_files)
#draft_file = 'ppW3MuNu_Run3Summer22EEDRPremix_cfg_draft.py' 
#final_file = 'ppW3MuNu_Run3Summer22EEDRPremix_cfg.py' 

#with open(input_to_crab, 'r') as file:
#    file_contents = file.read()
#    updated_contents = file_contents.replace(string_to_replace, formatted_files)
#
#    with open(final_file, 'w') as file:
#        file.write(updated_contents)
