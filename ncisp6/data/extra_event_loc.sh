# This bash script is to extra events location from SAC files
# Mao 20230301
#!/bin/bash

export SAC_DISPLAY_COPYRIGHT=0

homedir=$(pwd)
echo $homedir

for file in `ls $homedir/raw/`
do
    echo $file
    cd $homedir/raw/$file

    sac << EOF
    r NE030.01.1*sac
    lh gcarc 
    q
EOF
done

