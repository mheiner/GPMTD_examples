## source this script with a trailing integer argument to identify model runs for post-processing

julia postsim_to_R.jl $1 &
wait
Rscript --vanilla postProcess_loop.R $1 &
wait
