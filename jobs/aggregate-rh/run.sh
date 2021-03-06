#! /bin/bash
src/springs \
   --journal=results/aggregate-rh \
   --schira-param-file=jobs/aggregate-rh/initial-schira-params.txt \
   --schira-minimize-stepsize=0 \
   --data-minimize-stepsize=0.005 \
   --schira-cutoff=0.25 \
   --spring-cutoff=0.015 \
   --schira-strength=5.0 \
   --spring-strength=1.0 \
   --initial-energy=5.0 \
   --dampen=0.999 \
   --annealing-iterations=4 \
   --steps=5000 \
   "$@" \
   jobs/aggregate-rh/start.bin \
   results/aggregate-rh/result.bin
