# Scripts instructions

1. Choose an example
2. Run script 01 to generate benchmark data (requires large amount of disk space)
3. Run script 02 to discretize data 
4. Run script 03 to generate MSMs for MSM/RD (require pyemma and it can be memory intensive)
5. Run the rest of the scripts (CPU and memory intensive, might take a long time)
6. Once all of the scripts are run, you can go to the notebooks on  `examples/benchmark_plots/` and produce the final plots.
7. Note some plots might require you to do small modifictions of scripts and run them. For instance, the scripts with `bound2bound` suffix need that you specify initial and final bound states.
