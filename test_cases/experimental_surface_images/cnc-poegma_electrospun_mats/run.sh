#!/bin/bash

python src/characteristic_length.py \
--file poegma-cnc_1.65wtpct_elecspun_aligned_200x.png \
--dconv True \
--pop_num 2 \
--bar_len 20 \
--bar_pxl 82.625 \
--dof_lo_sigma 1.0 \
--dof_hi_sigma None \
--canny_sigma 0.5 \
--lo_len_lim_pop1 5 \
--hi_len_lim_pop1 10 \
--lo_len_lim_pop2 0 \
--hi_len_lim_pop2 2.5 \
--theta_lims 90 180
