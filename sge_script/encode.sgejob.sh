# Written by Duy Le - ledduy@ieee.org
# Last update Jun 26, 2012
#!/bin/sh
# Force to use shell sh. Note that #$ is SGE command
#$ -S /bin/sh
# Force to limit hosts running jobs
#$ -q all.q@@bc4hosts,all.q@@bc5hosts
# Log starting time
date 
# include your library here
export LD_LIBRARY_PATH=/net/per900a/raid0/plsang/usr.local/lib:/usr/local/lib:$LD_LIBRARY_PATH
# display your command here
echo [$HOSTNAME] [$JOB_ID] [matlab -nodisplay -r "encode_block( '$1', '$2', $3, $4 )"]
# change to your code directory here
cd /net/per900a/raid0/plsang/tools/kaori-secode-med
# Log info of current dir
pwd
# run your command with parameters ($1, $2,...) here, string variable is put in ' '
matlab -nodisplay -r "encode_block( '$1', '$2', $3, $4 )"
# Log ending time
date
