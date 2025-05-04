#!/bin/bash

python src/characteristic_length.py \
--file a1XG3p_1000x_003.tif \
--dconv True \
--pop_num 1 \
--bar_len 50 \
--bar_pxl 346 \
--dof_lo_sigma 1.0 \
--dof_hi_sigma 10 \
--canny_sigma 2.0 \
--lo_len_lim_pop1 0.4 \
--hi_len_lim_pop1 1.3 \
--lo_len_lim_pop2 3.4 \
--hi_len_lim_pop2 4.5
