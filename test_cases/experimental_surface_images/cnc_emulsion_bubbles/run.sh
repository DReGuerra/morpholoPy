#!/bin/bash

python src/characteristic_length.py \
--file 30_-S-In-H-Batch1-Day7-1.jpg \
--dconv True \
--pop_num 2 \
--bar_len 50 \
--dof_lo_sigma 3.0 \
--dof_hi_sigma 12 \
--canny_sigma 1.4 \
--lo_len_lim_pop1 0.8 \
--hi_len_lim_pop1 1.5 \
--lo_len_lim_pop2 0.16 \
--hi_len_lim_pop2 0.4
