#!/usr/bin/python

import os
import re
import random
import shlex
import getpass

NSTARTJOB=0
NJOBS=100
FILENAME='FS_nx_0_ny_0_fx_0_fy_0_'
DIRNAME='root://eoscms//eos/cms/store/temp/user/cbaus/castorgeotest/minbias_new_overlapfix/'
CFG = 'FS_minbias'

for job in range(NSTARTJOB,NSTARTJOB+NJOBS):
  run = 'cmsRun ' + CFG + '.py'
  job_script_name = 'job-' + CFG + '-%d.sh' % job
  f = open(job_script_name, "w")
  f.write("#!/bin/bash\n")
  f.write("#\n#\n")
  f.write('echo --------------------------------------------------------------------------\n')
  f.write('echo    Begin ' + job_script_name + ' on `hostname` @ `date`\n')
  f.write('echo --------------------------------------------------------------------------\n')
  
  f.write("cd ~/Geometry/CMSSW_7_1_0_pre3/src/\n")
  f.write("eval `scramv1 runtime -sh`\n")
  f.write("cd ~/Geometry/CMSSW_7_1_0_pre3/test/\n")

#  f.write("mkdir -p " + os.path.join(PROGRAM_OUT_DIR) + "\n")  
#  f.write("cd " + os.path.join(PROGRAM_DIR)  + "\n")
#  outfile = os.path.join(PROGRAM_OUT_DIR,OUTPUT_ROOTFILE_NAME % job)

#  f.write('nsls ' + rootfile + '\n')
#  f.write('if [ $? -eq 0 ]; then exit 0; fi \n')
  tmpdir = os.path.join('/tmp/cbaus/', 'job_'+ CFG)
  f.write("mkdir -p " + tmpdir + "\n")  
  f.write(run + ' ' + tmpdir + '/' + FILENAME + '%d.root \n' % job)
  f.write('/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select cp ' + tmpdir + '/' + FILENAME + '%d.root' %job + ' ' + DIRNAME + '\n')
  
  f.write('echo -------------------------------------\n')
  f.write('echo            Job OK @ `date`\n')
  f.write('echo -------------------------------------\n')
  f.close()

  os.system("chmod +x " + job_script_name)
  os.system('bsub -q 8nh ' + job_script_name)
  #os.system("rm -f " + job_script_name) #doesn't work
