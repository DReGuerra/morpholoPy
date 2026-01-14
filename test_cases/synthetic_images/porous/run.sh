#!/bin/bash

# python3 src/porous_media.py

python3 src/characteristic_length.py \
--file bubbles_pores_1000x1000.png \
--dconv True \
--pop_num 1 \
--bar_len 1 \
--bar_pxl 100 \
--dof_lo_sigma 3.0 \
--dof_hi_sigma 6 \
--canny_sigma 1.4 \
--lo_len_lim_pop1 0 \
--hi_len_lim_pop1 1.5 \
--lo_len_lim_pop2 0 \
--hi_len_lim_pop2 1.5
