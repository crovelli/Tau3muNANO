import glob
import subprocess
import re
import argparse
import pandas as pd
from color import color as c


def parse_CRABstatus(shell_output):
    # Define regex patterns to capture job statuses
    # typical output is 'finished     		 98.1% (1166/1188)'
    failed_pattern      = re.compile(r'failed\s+(\d+.\d+)%\s+\(\s*(\d+)/(\d+)\)')
    finished_pattern    = re.compile(r'finished\s+(\d+.\d+)%\s+\(\s*(\d+)/(\d+)\)')
    running_pattern     = re.compile(r'running\s+(\d+.\d+)%\s+\(\s*(\d+)/(\d+)\)')

    # Search for the patterns in the shell output
    failed_match = failed_pattern.search(shell_output)
    finished_match = finished_pattern.search(shell_output)
    running_match = running_pattern.search(shell_output)

    # Extract the values
    failed_jobs, completed_jobs, running_jobs, total_jobs = 0,0,0,0
    if failed_match:
        if (debug): print(f'[+] Failed match: {failed_match.groups()}')
        failed_jobs = int(failed_match.group(2))

    if finished_match:
        if (debug):print(f'[+] Finished match: {finished_match.groups()}')
        completed_jobs = int(finished_match.group(2))
        total_jobs     = int(finished_match.group(3))
    

    if running_match:
        if (debug):print(f'[+] Running match: {running_match.groups()}')
        running_jobs = int(running_match.group(2))

    # Total jobs should be consistent across all matches
    if finished_match and running_match:
        assert total_jobs == int(finished_match.group(3)) == int(running_match.group(3))

    return {
        'failed_jobs': failed_jobs,
        'completed_jobs': completed_jobs,
        'total_jobs': total_jobs
    }

def parse_CRABreport(shell_output):
    # Define regex patterns to capture job statuses
    # typical output is 'finished     		 98.1% (1166/1188)'
    file_processed    = re.compile(r'  Number of files processed: (\d+)')
    events_read       = re.compile(r'  Number of events read: (\d+)')
    events_saved      = re.compile(r'  Number of events written in EDM files: (\d+)')

    # Search for the patterns in the shell output
    file_processed_match = file_processed.search(shell_output)
    events_read_match = events_read.search(shell_output)
    events_saved_match = events_saved.search(shell_output)

    # Extract the values
    file_processed, events_read, events_saved = 0,0,0
    if file_processed_match:
        if (debug): print(f'[+] Files processed match: {file_processed_match.groups()}')
        file_processed = int(file_processed_match.group(1))
    
    if events_read_match:
        if (debug):print(f'[+] Events read match: {events_read_match.groups()}')
        events_read = int(events_read_match.group(1))
    
    if events_saved_match:
        if (debug):print(f'[+] Events saved match: {events_saved_match.groups()}')
        events_saved = int(events_saved_match.group(1))
    
    return {
        'file_processed': file_processed,
        'events_read': events_read,
        'events_saved': events_saved
    }

parser = argparse.ArgumentParser(description='Check the status of CRAB jobs')
parser.add_argument('--debug', action='store_true', help='Enable debugging mode')
parser.add_argument('--input_dir', type=str, help='Directory containing CRAB jobs')
parser.add_argument('--resub_threshold', type=float, default=0.05, help='Threshold for resubmitting jobs')
parser.add_argument('--produce_report', action='store_true', help='Produce a report of the job')
args = parser.parse_args()
debug = args.debug
print('\n')

# Set the output report if needed
if args.produce_report:
    columns = ['job_id', 'jobs_done', 'total_jobs', 'processed_files', 'events_read', 'events_saved']
    report_df = pd.DataFrame(columns=columns)


# List all the directories containing CRAB jobs
print(f'{c.BOLD}[+] cheking crab jobs in {args.input_dir}')
print(f'----------------------------------------------------------------------{c.END}\n')
dir_list = glob.glob(f'{args.input_dir}/*')
# Loop over the crab jobs irectories and check the status
for dir in dir_list:
#dir = dir_list[0]
    print(f" - {c.BOLD}{dir}{c.END}")
    is_repor_complete = False
    # status
    command = ['crab', 'status', '-d', dir]
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    shell_output = result.stdout
    if (debug) :print(shell_output)

    # Parse the shell output
    job_info = parse_CRABstatus(shell_output)

    # Print the retrieved information
    print(f"[i] Job status for {command[3]}")
    print(f"\tFailed jobs: {job_info['failed_jobs']}")
    print(f"\tCompleted jobs: {job_info['completed_jobs']}")
    print(f"\tTotal jobs: {job_info['total_jobs']}")
    if (job_info['total_jobs'] == 0):
        print(f"[status] No jobs found in the directory --> {c.RED}skipping{c.END}\n")
        continue
    # if above the threshold resubmit jobs
    failed_fraction = job_info['failed_jobs'] / job_info['total_jobs']
    if (failed_fraction > 0.05):
        print(f"[status] fraction of failed jobs is {failed_fraction} which is above the threshold of 0.05 --> r{c.RED}esubmitting jobs{c.END}")
        resubmit_command = ['crab', 'resubmit', '-d', dir]
        subprocess.run(resubmit_command)
    # else parse the job report
    else:
        print(f"[status] fraction of failed jobs is {failed_fraction} which is below the threshold of 0.05 --> {c.GREEN}not resubmitting jobs{c.END}")
        report_command = ['crab', 'report', '-d', dir]
        report = subprocess.run(report_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        report_output = report.stdout
        if (debug): print(report_output)
        
        job_report = parse_CRABreport(report_output)
        print(f"[i] Job report for {command[3]}")
        print(f"\tFiles processed: {job_report['file_processed']}")
        print(f"\tEvents read: {job_report['events_read']}")
        print(f"\tEvents saved: {job_report['events_saved']}")
        is_repor_complete = True    
    # save report only if the report is complete
    if (args.produce_report and is_repor_complete):
        tmp = pd.DataFrame({
            'job_id': dir,
            'jobs_done': job_info['completed_jobs'],
            'total_jobs': job_info['total_jobs'],
            'processed_files': job_report['file_processed'],
            'events_read': job_report['events_read'],
            'events_saved': job_report['events_saved']
        }, columns=columns, index=[0])
        report_df = pd.concat([report_df, tmp], ignore_index=True)

    print('\n')
    

if (args.produce_report):       
    print(f"[i] Saving the report to {args.input_dir}/report.csv")
    report_df.to_csv(f'{args.input_dir}/report.csv', index=False)

