#!/bin/bash

python src/characteristic_length.py \
--file kdf_biaxial_20um.tif \
--dconv True \
--pop_num 1 \
--bar_len 20 \
--bar_pxl 170 \
--dof_lo_sigma 1.1 \
--dof_hi_sigma None \
--canny_sigma 1.4 \
--lo_len_lim_pop1 1.5 \
--hi_len_lim_pop1 3 \
--lo_len_lim_pop2 3.4 \
--hi_len_lim_pop2 4.5
