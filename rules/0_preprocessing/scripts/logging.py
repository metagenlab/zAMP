### Functions to generate logs

from shutil import copyfile
import sys
import datetime
import os
import getpass


## Define log directory path
if "logging_folder" not in config.keys():
    logging_folder = "logs/"
else:
    logging_folder=config["logging_folder"]

if not logging_folder.endswith("/"):
    logging_folder = logging_folder + "/"

## Generate a day_hour folder to store execution specific parameters (command, config fie, sample sheet)
today = datetime.datetime.now()
date = today.strftime("%Y_%m_%d")
time = today.strftime('%H_%M_%S_%f')[:-4]
time_tag = str(date + "_" + time)
timed_logging_folder = logging_folder + "/" + time_tag + "/"

## Except for dryrun execution, create a log folder. This will be filled by execution parameters and be provided to all rules to record their standard ouput or specific log files. 
if "--dryrun" not in sys.argv and "-n" not in sys.argv:

    os.makedirs(logging_folder, exist_ok=True)

    ## Copy the command line, the config file and the sample sheet in a subdirectory of the log directory
    ### Log command line
    cmd_file = timed_logging_folder + "cmd.txt"
    with open(cmd_file, "w") as f:
        f.write(" ".join(sys.argv)+"\n")

    ### Log config file
    if workflow.overwrite_configfiles[0] is not None:
            copyfile(workflow.overwrite_configfiles[0], timed_logging_folder+"config.yaml")

    ### Log local or SRA sample sheet
    if "local_samples" in config.keys():
        copyfile(config["local_samples"], timed_logging_folder+"local_samples.tsv")
    if "sra_samples" in config.keys():
        copyfile(config["sra_samples"], timed_logging_folder+"sra_samples.tsv")

    ### Log git hash of the pipeline
    git_log_path = timed_logging_folder + "git.txt"
    with open(git_log_path, "w") as devnull:
        subprocess.run(args=["git -C " + workflow.basedir + " rev-parse HEAD"], shell=True, stdout=devnull)   

    ### Log user ID  
    user_path = timed_logging_folder + "user.txt"
    user_cred = getpass.getuser()
    with open(user_path, "w") as f:
        f.write(user_cred +"\n")

