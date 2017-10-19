
import os
import re
import socket
import subprocess
from tqdm import tqdm


def copyToBaobabNodes(file_to_copy):
    """ 
    Copy a given file to each of the 6 baobab nodes
    """
    print "Copying %s to the 6 baobab nodes" % (file_to_copy)
    print "\tTip: If you are being asked for a password, setup ssh keys to allow passwordless ssh and scp"
    ip_template = "192.168.200.%s"
    # use the localdisk storage on the baobab nodes to speed up read/write times
    copy_to_dir = re.sub("^/data", "/localdisk", os.path.dirname(os.path.abspath(file_to_copy))) 
    # loop through the 6 nodes
    for i in tqdm(range(1,7)):
        command = "ssh %s 'mkdir -p %s'; scp %s %s:%s" % (ip_template%i, copy_to_dir, os.path.abspath(file_to_copy), ip_template%i, copy_to_dir)
        if 'baobab' in socket.gethostname():
            subprocess.check_call(command, shell=True)
        else:
            runCommandOnBaobab(command)


def runCommandOnBaobab(command):
    """ 
    Run a given command on baobab. Useful for copying a file from baobab to each of the nodes
    """
    command = "ssh -t baobab.cbb.lan \"%s\"" % (command)
    print "Running: %s" % (command)
    subprocess.check_call(command, shell=True)
