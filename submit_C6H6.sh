file="C6H6_WI.in"
Executable=switchblade.exe

#!/bin/bash	  
#$ -V	                   #Inherit the submission environment
#$ -cwd	            # Start job in submission directory
#$ -R y
#$ -N C6H6-WI	            # Job Name
#$ -j y	            # Combine stderr and stdout
#$ -o $JOB_NAME.o$JOB_ID  # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way 240	     # Requests 12 tasks/node, 24 cores total
#$ -q normal	            # Queue name normal
#$ -l h_rt=4:00:00       # Run time (hh:mm:ss) - 1.5 hours
	                  
set -x	                   # Echo commands, use set echo with csh

export MV2_USE_SHMEM_COLL=0
run="ibrun ./$Executable $file"	          
echo $run
$run
