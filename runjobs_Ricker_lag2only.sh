# Ricker2 lag2 only : y-length 1,000, X is 1,000 by 10

nice julia fitGPMTD.jl 427191 100 5 Ricker_lag2only &> logs/r1.txt &
nice julia fitGPMTD.jl 427192 100 5 Ricker_lag2only &> logs/r2.txt &
nice julia fitGPMTD.jl 427193 100 5 Ricker_lag2only &> logs/r3.txt &

#wait

nice julia fitGPMTD.jl 427191 500 2 Ricker_lag2only &> logs/r4.txt &
nice julia fitGPMTD.jl 427192 500 2 Ricker_lag2only &> logs/r5.txt &
nice julia fitGPMTD.jl 427193 500 2 Ricker_lag2only &> logs/r6.txt &

wait
