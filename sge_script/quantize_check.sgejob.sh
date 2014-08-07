# Written by Duy Le - ledduy@ieee.org
# Last update Jun 26, 2012
#!/bin/sh
# Force to use shell sh. Note that #$ is SGE command
#$ -S /bin/sh
# Force to limit hosts running jobs
#$ -q all.q@@bc4hosts
# Log starting time
date 
# change to your code directory here
cd /net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/code
# Log info of current dir
pwd
# run your command with parameters ($1, $2,...) here, string variable is put in ' '
matlab -nodisplay -r "quantize_check( $1, $2 )"
# Log ending time
date
