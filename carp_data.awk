BEGIN {FS=","
       REDUCE=0
       STEP=10.0
       TSYS=120
       TMIN=10
       }
/./ {lmst_hour=$4
     lmst_min=$5
     lmst_sec=$6
     
     lmst = lmst_hour + (lmst_min/60) + (lmst_sec/3600.0)
     ndx = (lmst*3600.0)/STEP
     ndx = int(ndx)
     
     value = 0.0
     for (i = 9+REDUCE; i <= (NF-REDUCE); i++)
     {
		value += $i
	 }
	 values[ndx] += value
	 valuecnt[ndx] += 1
	 }
END {
        av = -1
        a = 0.3
        b = 1.0 - a
        minval = 99999.0
        for (i = 0; i < int(86400/STEP); i++)
        {
            if (valuecnt[i] >= 1)
            {
                 if (av == -1)
                 {
					av = values[i]/valuecnt[i]
				 }
                 av = a*(values[i]/valuecnt[i])+(b*av)
                 newvalues[i] = av
                 if (av < minval)
                 {
					minval = av
				 }
            }
        }
        for (i = 0; i < int(86400/STEP); i++)
        {
			if (valuecnt[i] >= 1)
			{
			    tv = (newvalues[i]/minval)
			    tv *= (TSYS+TMIN)
			    tv -= TSYS
				printf ("%f %.5e\n", (i*STEP)/3600.0, tv)
			}
	    }
     }
