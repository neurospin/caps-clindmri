#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
FreeSurfer Color Look Up Tables for the regions of interest, extracted from
FreeSurferColorLUT.txt
"""

# -----------------------------------------------------------------------------
# Subcortical regions of interest, to read from "aseg.mgz"

ASEG_RAW_LUT = \
"""
10  Left-Thalamus-Proper                    0   118 14  0
11  Left-Caudate                            122 186 220 0
12  Left-Putamen                            236 13  176 0
13  Left-Pallidum                           12  48  255 0
17  Left-Hippocampus                        220 216 20  0
18  Left-Amygdala                           103 255 255 0
26  Left-Accumbens-area                     255 165 0   0
28  Left-VentralDC                          165 42  42  0
49  Right-Thalamus-Proper                   0   118 14  0
50  Right-Caudate                           122 186 220 0
51  Right-Putamen                           236 13  176 0
52  Right-Pallidum                          13  48  255 0
53  Right-Hippocampus                       220 216 20  0
54  Right-Amygdala                          103 255 255 0
58  Right-Accumbens-area                    255 165 0   0
60  Right-VentralDC                         165 42  42  0
""".strip().split("\n")

# -----------------------------------------------------------------------------
# Cortical regions of interest, for 2 atlases:
# - Desikan (aparc+aseg.mgz)
# - Destrieux (aparc.a2009s+aseg.mgz)

# Only the labels for the left hemisphere cortex.
# The other labels are inferred according to the rules:
# <right cortex label> = <left cortex label> + 1000
# <white matter label> = <cortex label> + 2000

APARC_ASEG_RAW_LUT = \
"""
1001    ctx-lh-bankssts                     25  100 40  0
1002    ctx-lh-caudalanteriorcingulate      125 100 160 0
1003    ctx-lh-caudalmiddlefrontal          100 25  0   0
1004    ctx-lh-corpuscallosum               120 70  50  0
1005    ctx-lh-cuneus                       220 20  100 0
1006    ctx-lh-entorhinal                   220 20  10  0
1007    ctx-lh-fusiform                     180 220 140 0
1008    ctx-lh-inferiorparietal             220 60  220 0
1009    ctx-lh-inferiortemporal             180 40  120 0
1010    ctx-lh-isthmuscingulate             140 20  140 0
1011    ctx-lh-lateraloccipital             20  30  140 0
1012    ctx-lh-lateralorbitofrontal         35  75  50  0
1013    ctx-lh-lingual                      225 140 140 0
1014    ctx-lh-medialorbitofrontal          200 35  75  0
1015    ctx-lh-middletemporal               160 100 50  0
1016    ctx-lh-parahippocampal              20  220 60  0
1017    ctx-lh-paracentral                  60  220 60  0
1018    ctx-lh-parsopercularis              220 180 140 0
1019    ctx-lh-parsorbitalis                20  100 50  0
1020    ctx-lh-parstriangularis             220 60  20  0
1021    ctx-lh-pericalcarine                120 100 60  0
1022    ctx-lh-postcentral                  220 20  20  0
1023    ctx-lh-posteriorcingulate           220 180 220 0
1024    ctx-lh-precentral                   60  20  220 0
1025    ctx-lh-precuneus                    160 140 180 0
1026    ctx-lh-rostralanteriorcingulate     80  20  140 0
1027    ctx-lh-rostralmiddlefrontal         75  50  125 0
1028    ctx-lh-superiorfrontal              20  220 160 0
1029    ctx-lh-superiorparietal             20  180 140 0
1030    ctx-lh-superiortemporal             140 220 220 0
1031    ctx-lh-supramarginal                80  160 20  0
1032    ctx-lh-frontalpole                  100 0   100 0
1033    ctx-lh-temporalpole                 70  70  70  0
1034    ctx-lh-transversetemporal           150 150 200 0
1035    ctx-lh-insula                       255 192 32  0
""".strip().split("\n")


APARC_A2009S_ASEG_RAW_LUT = \
"""
11101  ctx_lh_G_and_S_frontomargin             23 220  60   0
11102  ctx_lh_G_and_S_occipital_inf            23  60 180   0
11103  ctx_lh_G_and_S_paracentral              63 100  60   0
11104  ctx_lh_G_and_S_subcentral               63  20 220   0
11105  ctx_lh_G_and_S_transv_frontopol         13   0 250   0
11106  ctx_lh_G_and_S_cingul-Ant               26  60   0   0
11107  ctx_lh_G_and_S_cingul-Mid-Ant           26  60  75   0
11108  ctx_lh_G_and_S_cingul-Mid-Post          26  60 150   0
11109  ctx_lh_G_cingul-Post-dorsal             25  60 250   0
11110  ctx_lh_G_cingul-Post-ventral            60  25  25   0
11111  ctx_lh_G_cuneus                        180  20  20   0
11112  ctx_lh_G_front_inf-Opercular           220  20 100   0
11113  ctx_lh_G_front_inf-Orbital             140  60  60   0
11114  ctx_lh_G_front_inf-Triangul            180 220 140   0
11115  ctx_lh_G_front_middle                  140 100 180   0
11116  ctx_lh_G_front_sup                     180  20 140   0
11117  ctx_lh_G_Ins_lg_and_S_cent_ins          23  10  10   0
11118  ctx_lh_G_insular_short                 225 140 140   0
11119  ctx_lh_G_occipital_middle              180  60 180   0
11120  ctx_lh_G_occipital_sup                  20 220  60   0
11121  ctx_lh_G_oc-temp_lat-fusifor            60  20 140   0
11122  ctx_lh_G_oc-temp_med-Lingual           220 180 140   0
11123  ctx_lh_G_oc-temp_med-Parahip            65 100  20   0
11124  ctx_lh_G_orbital                       220  60  20   0
11125  ctx_lh_G_pariet_inf-Angular             20  60 220   0
11126  ctx_lh_G_pariet_inf-Supramar           100 100  60   0
11127  ctx_lh_G_parietal_sup                  220 180 220   0
11128  ctx_lh_G_postcentral                    20 180 140   0
11129  ctx_lh_G_precentral                     60 140 180   0
11130  ctx_lh_G_precuneus                      25  20 140   0
11131  ctx_lh_G_rectus                         20  60 100   0
11132  ctx_lh_G_subcallosal                    60 220  20   0
11133  ctx_lh_G_temp_sup-G_T_transv            60  60 220   0
11134  ctx_lh_G_temp_sup-Lateral              220  60 220   0
11135  ctx_lh_G_temp_sup-Plan_polar            65 220  60   0
11136  ctx_lh_G_temp_sup-Plan_tempo            25 140  20   0
11137  ctx_lh_G_temporal_inf                  220 220 100   0
11138  ctx_lh_G_temporal_middle               180  60  60   0
11139  ctx_lh_Lat_Fis-ant-Horizont             61  20 220   0
11140  ctx_lh_Lat_Fis-ant-Vertical             61  20  60   0
11141  ctx_lh_Lat_Fis-post                     61  60 100   0
11142  ctx_lh_Medial_wall                      25  25  25   0
11143  ctx_lh_Pole_occipital                  140  20  60   0
11144  ctx_lh_Pole_temporal                   220 180  20   0
11145  ctx_lh_S_calcarine                      63 180 180   0
11146  ctx_lh_S_central                       221  20  10   0
11147  ctx_lh_S_cingul-Marginalis             221  20 100   0
11148  ctx_lh_S_circular_insula_ant           221  60 140   0
11149  ctx_lh_S_circular_insula_inf           221  20 220   0
11150  ctx_lh_S_circular_insula_sup            61 220 220   0
11151  ctx_lh_S_collat_transv_ant             100 200 200   0
11152  ctx_lh_S_collat_transv_post             10 200 200   0
11153  ctx_lh_S_front_inf                     221 220  20   0
11154  ctx_lh_S_front_middle                  141  20 100   0
11155  ctx_lh_S_front_sup                      61 220 100   0
11156  ctx_lh_S_interm_prim-Jensen            141  60  20   0
11157  ctx_lh_S_intrapariet_and_P_trans       143  20 220   0
11158  ctx_lh_S_oc_middle_and_Lunatus         101  60 220   0
11159  ctx_lh_S_oc_sup_and_transversal         21  20 140   0
11160  ctx_lh_S_occipital_ant                  61  20 180   0
11161  ctx_lh_S_oc-temp_lat                   221 140  20   0
11162  ctx_lh_S_oc-temp_med_and_Lingual       141 100 220   0
11163  ctx_lh_S_orbital_lateral               221 100  20   0
11164  ctx_lh_S_orbital_med-olfact            181 200  20   0
11165  ctx_lh_S_orbital-H_Shaped              101  20  20   0
11166  ctx_lh_S_parieto_occipital             101 100 180   0
11167  ctx_lh_S_pericallosal                  181 220  20   0
11168  ctx_lh_S_postcentral                    21 140 200   0
11169  ctx_lh_S_precentral-inf-part            21  20 240   0
11170  ctx_lh_S_precentral-sup-part            21  20 200   0
11171  ctx_lh_S_suborbital                     21  20  60   0
11172  ctx_lh_S_subparietal                   101  60  60   0
11173  ctx_lh_S_temporal_inf                   21 180 180   0
11174  ctx_lh_S_temporal_sup                  223 220  60   0
11175  ctx_lh_S_temporal_transverse           221  60  60   0
""".strip().split("\n")
