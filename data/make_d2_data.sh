#!/bin/sh

dir_d2="./deshima2.0"
# D2_supplements=' ./D2_supplements'
data_d2=`ls -1 $dir_d2`

for d in $data_d2; do
    cd $dir_d2/$d
    obsid=`echo $d | awk '{sub(/cosmos_/, "")}1' `
    ant_time_data=`less ${obsid}.ant | awk '!/^#/ && !/^2090/ {print $1}'`
    # make .misty
    t_max=`echo "$ant_time_data" | sort -nr | head -n 1`
    t_min=`echo "$ant_time_data" | sort -n | head -n 1`
    NR=`echo "$ant_time_data" | wc -l `
    echo "$obsid $t_max $t_min $NR"
    echo "# MiSTY test data" > $obsid.misty
    echo '#X Y timestamp framecnt Z flag' > ${obsid}.skychop
    echo '# /data/spacekids/local/bin/des_chopper_logger.py /data/spacekids/data/log_chile/des_chopper/chopper --ip 192.168.2.181 --port 4001 --interval 1 --n_heaters 0 --n_motors 1 --rotate_daily# /data/spacekids/local/bin/des_chopper_logger.py /data/spacekids/data/log_chile/des_chopper/chopper --ip 192.168.2.181 --port 4001 --interval 1 --n_heaters 0 --n_motors 1 --rotate_daily\n# [ 1] C Time\n# [ 2] motor0: 126 127    ready\n# [ 3] motor0: 126 127    move\n# [ 4] motor0: 126 127    inposition\n# [ 5] motor0: 128 129    current alarm code\n# [ 6] motor0: 194 195    current selected  data No.\n# [ 7] motor0: 196 197    current operation data No.\n# [ 8] motor0: 198 199    command position\n# [ 9] motor0: 200 201    command speed [r/min]\n# [10] motor0: 202 203    command speed [Hz]\n# [11] motor0: 204 205    feedback position\n# [12] motor0: 206 207    feedback speed [r/min]\n# [13] motor0: 208 209    feedback speed [Hz]\n# [14] motor0: 210 211    remaining dwell time [msec]\n# [15] motor0: 212 213    direct I/O\n# [16] motor0: 214 215    torque monitor\n# [17] motor0: 218 219    cumulative load monitor\n# [18] motor0: 248 249    driver temperature [degC]\n# [19] motor0: 250 251    motor temperature  [degC]\n# [20] motor0: 252 253    odometer  [kRev]\n# [21] motor0: 254 255    tripmeter [kRev]\n' > ${obsid}.bbchop
    echo '# cabin temperature' > ${obsid}.cabin
 
    for t in ${ant_time_data}; do
        ut=`date -d "$(echo $t | sed 's/\(....\)\(..\)\(..\)\(..\)\(..\)\(..\)\(.*\)/\1-\2-\3 \4:\5:\6\7/')" +%s.%2N`
	
	echo $t | awk '{printf("%s %s %s %s %s %05.2f  73.886400 25.740000 -24.491380  87.137089   7.452227   6.156458 1\n", substr($0,1,4), substr($0,5,2), substr($0,7,2), substr($0,9,2), substr($0,11,2), substr($0,13,5), $0)}' >> $obsid.misty
	
	cabin_data=`echo $t | awk '{printf("%s/%s/%s %s:%s  13.2   16.6 9999.0 9999.0 9999.0 9999.0 9999.0 9999.0 9999.0 9999.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0\n", substr($0,1,4), substr($0,5,2), substr($0,7,2), substr($0,9,2), substr($0,11,2), $0)}'`
        if [ "$cabin_data" != "$cabin_data2" ]; then
            echo "$cabin_data" >> $obsid.cabin
	fi
        cabin_data2="$cabin_data"
	
	echo "103.0 82.6 $ut 29929 0.0 0" >> $obsid.skychop
	
	echo "$ut 0 1 0 00 210 210 294.20 3000 10.000 289.60 3000 10.000 0 80340000 17.5   36 26.80 37.30 223543 223543" >> $obsid.bbchop
    done 
    cd ../../
done


