#!/bin/sh

dir_d2="./deshima2.0"
# D2_supplements=' ./D2_supplements'
data_d2=`ls -1 $dir_d2`

for d in $data_d2; do
    cd $dir_d2/$d
    obsid=`echo -e $d | awk '{sub(/cosmos_/, "")}1' `
    ant_time_data=`less ${obsid}.ant | awk '!/^#/ && !/^2090/ {print $1}'`
    # make .tsky
    t_max=`echo -e "$ant_time_data" | sort -nr | head -n 1`
    t_min=`echo -e "$ant_time_data" | sort -n | head -n 1`
    NR=`echo -e "$ant_time_data" | wc -l `
    echo -e "$obsid $t_max $t_min $NR"
    echo -e "# MiSTi tsky test data" > $obsid.tsky
    echo -e "# [0] YYYY \n# [1] MM\n# [2] DD" >> $obsid.tsky
    echo -e "# [3] hh\n# [4] mm\n# [5] ss.ss" >> $obsid.tsky
    echo -e "# [6] az (deg)\n# [7] el (deg)" >> $obsid.tsky
    echo -e "# [8] ??? power (dBm)" >> $obsid.tsky
    echo -e "# [9]  ??? Hot load tempearture (deg): heated TK-RAM" >> $obsid.tsky
    echo -e "# [10] ??? Receiver room temperature (deg): Room absorber" >> $obsid.tsky
    echo -e "# [11] ??? Primary mirror room temperature (deg): Inside GORE-TEX membrane housing" >> $obsid.tsky
    echo -e "# [12] ??? Chopper mirror status, 0: Sky, 1: Room, 2: Hot" >> $obsid.tsky
    echo -e "# [13] ???" >> $obsid.tsky
    echo -e "# [14] ???" >> $obsid.tsky
    # make .pwv
    echo -e "# MiSTi pwv test data" > $obsid.pwv
    echo -e "# [0] (UTC) YYYY/MM/DD \n# [1] (UTC) hh:mm:ss.ss \n# [2] unixtime\n# [3] Az (deg)\n# [4] El (deg)\n# [5] PWV(um)\n# [6] Tground(K)" >> $obsid.pwv

    # make Skychop
    echo -e '#set velocity (rpm): 87.000000\n#set smapling rate (Hz): 1000.000000\n#ts state' > ${obsid}.skychop
    echo -e '# /data/spacekids/local/bin/des_chopper_logger.py /data/spacekids/data/log_chile/des_chopper/chopper --ip 192.168.2.181 --port 4001 --interval 1 --n_heaters 0 --n_motors 1 --rotate_daily# /data/spacekids/local/bin/des_chopper_logger.py /data/spacekids/data/log_chile/des_chopper/chopper --ip 192.168.2.181 --port 4001 --interval 1 --n_heaters 0 --n_motors 1 --rotate_daily\n# [ 1] C Time\n# [ 2] motor0: 126 127    ready\n# [ 3] motor0: 126 127    move\n# [ 4] motor0: 126 127    inposition\n# [ 5] motor0: 128 129    current alarm code\n# [ 6] motor0: 194 195    current selected  data No.\n# [ 7] motor0: 196 197    current operation data No.\n# [ 8] motor0: 198 199    command position\n# [ 9] motor0: 200 201    command speed [r/min]\n# [10] motor0: 202 203    command speed [Hz]\n# [11] motor0: 204 205    feedback position\n# [12] motor0: 206 207    feedback speed [r/min]\n# [13] motor0: 208 209    feedback speed [Hz]\n# [14] motor0: 210 211    remaining dwell time [msec]\n# [15] motor0: 212 213    direct I/O\n# [16] motor0: 214 215    torque monitor\n# [17] motor0: 218 219    cumulative load monitor\n# [18] motor0: 248 249    driver temperature [degC]\n# [19] motor0: 250 251    motor temperature  [degC]\n# [20] motor0: 252 253    odometer  [kRev]\n# [21] motor0: 254 255    tripmeter [kRev]' > ${obsid}.roomchop
    echo -e '# cabin temperature' > ${obsid}.cabin
    echo -e '# YYYY/mm/dd HH:MM UpperCabinTemperature MainCabinTemperature' > ${obsid}.cabin

    ut_start=`date -u -d "$(echo -e $t_min | sed 's/\(....\)\(..\)\(..\)\(..\)\(..\)\(..\)\(.*\)/\1-\2-\3 \4:\5:\6\7/')" +%s.%6N`
    ut_end=`date -u -d "$(echo -e $t_max | sed 's/\(....\)\(..\)\(..\)\(..\)\(..\)\(..\)\(.*\)/\1-\2-\3 \4:\5:\6\7/')" +%s.%6N`
	
    for t in ${ant_time_data}; do
        ut=`date -u -d "$(echo -e $t | sed 's/\(....\)\(..\)\(..\)\(..\)\(..\)\(..\)\(.*\)/\1-\2-\3 \4:\5:\6\7/')" +%s.%2N`
	
	echo -e $t | awk '{printf("%s %s %s %s %s %05.2f 180.000 90.000  69.360 -11.952 -12.047  4.3210e+03 0.0000e+00  4.3679e+01 0.0000e+00\n", substr($0,1,4), substr($0,5,2), substr($0,7,2), substr($0,9,2), substr($0,11,2), substr($0,13,5), $0)}' >> $obsid.tsky
	echo -e $t | awk '{printf("%s/%s/%s %s:%s:%05.2f  180.000 90.000   610.0 258.0\n", substr($0,1,4), substr($0,5,2), substr($0,7,2), substr($0,9,2), substr($0,11,2), substr($0,13,5), $0)}' >> $obsid.pwv
	
	cabin_data=`echo -e $t | awk '{printf("%s/%s/%s %s:%s  13.2   16.6 9999.0 9999.0 9999.0 9999.0 9999.0 9999.0 9999.0 9999.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0\n", substr($0,1,4), substr($0,5,2), substr($0,7,2), substr($0,9,2), substr($0,11,2), $0)}'`
        if [ "$cabin_data" != "$cabin_data2" ]; then
            echo -e "$cabin_data" >> $obsid.cabin
	fi
        cabin_data2="$cabin_data"
	
	
	echo -e "$ut 0 1 0 00 210 210 294.20 3000 10.000 289.60 3000 10.000 0 80340000 17.5   36 26.80 37.30 223543 223543" >> $obsid.roomchop
    done 
   
    for ut in `seq $ut_start 0.0010 $ut_end`; do
	echo -e "$ut 1" >> $obsid.skychop
    done
    #gzip $obsid.skychop    
    
    cd ../../
done


