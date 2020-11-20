### Functions to generate logs

from shutil import copyfile
import sys
import datetime
import os
import getpass

try:
    pipeline_path 
except NameError:
    pipeline_path = workflow.basedir + "/../"

if "logging_folder" not in config.keys():
    logging_folder = "logs/"
else:
    logging_folder=config["logging_folder"]

if "--dryrun" not in sys.argv and "-n" not in sys.argv:
    today = datetime.datetime.now()
    date = today.strftime("%Y/%m/%d")
    time = today.strftime('%H_%M_%S_%f')[:-4]
    logging_folder = str(logging_folder+"/"+date+"/"+time).replace("//", "/")
    if not logging_folder.endswith("/"):
        logging_folder = logging_folder + "/"
    os.makedirs(logging_folder, exist_ok=True)
    cmd_file = logging_folder + "/cmd.txt"
    with open(cmd_file, "w") as f:
        f.write(" ".join(sys.argv)+"\n")
    if workflow.overwrite_configfiles[0] is not None:
            copyfile(workflow.overwrite_configfiles[0], logging_folder+"/config.yaml")
    if "local_samples" in config.keys():
        copyfile(config["local_samples"], logging_folder+"/local_samples.tsv")
    if "sra_samples" in config.keys():
        copyfile(config["sra_samples"], logging_folder+"/sra_samples.tsv")
    git_log_path = logging_folder + "/git.txt"
    with open(git_log_path, "w") as devnull:
        subprocess.run(args=["git -C " + workflow.basedir + " rev-parse HEAD"], shell=True, stdout=devnull)     
    user_path = logging_folder + "/user.txt"
    user_cred = getpass.getuser()
    with open(user_path, "w") as f:
        f.write(user_cred +"\n")

