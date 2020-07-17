# Ricker2d : y-length 1,000, X is 1,000 by 10

nice julia fitGPMTD.jl 506191 100 5 Ricker2d_log &> logs/r2d1.txt &
nice julia fitGPMTD.jl 506192 100 5 Ricker2d_log &> logs/r2d2.txt &
nice julia fitGPMTD.jl 506193 100 5 Ricker2d_log &> logs/r2d3.txt &

nice julia fitGPMTD.jl 506191 500 5 Ricker2d_log &> logs/r2d4.txt &
nice julia fitGPMTD.jl 506192 500 5 Ricker2d_log &> logs/r2d5.txt &
nice julia fitGPMTD.jl 506193 500 5 Ricker2d_log &> logs/r2d6.txt &

wait
