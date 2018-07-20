#!/bin/bash

python opal.py -g cami_ii_mg/gs.profile \
cami_ii_mg/commonkmers.profile \
cami_ii_mg/focus.profile \
cami_ii_mg/metaphlan2.profile \
cami_ii_mg/metaphyler.profile \
cami_ii_mg/quikr.profile \
cami_ii_mg/tipp.profile \
cami_ii_mg/motu.profile \
../phat_unconstr.profile \
../phat_constr.profile \
-l "CommonKmers, FOCUS, Metaphlan, MetaPhyler, Quikr, TIPP, mOTU, PhAT (U), PhAT (C)" \
-o opal_output \
-d "2nd CAMI Challenge Mouse Gut Toy Dataset"

#python opal.py -g cami_ii_mg/gs.profile \
#cami_ii_mg/commonkmers.profile \
#cami_ii_mg/focus.profile \
#cami_ii_mg/metaphlan2.profile \
#cami_ii_mg/metaphyler.profile \
#cami_ii_mg/quikr.profile \
#cami_ii_mg/tipp.profile \
#cami_ii_mg/motu.profile \
#../phat_unconstr.profile \
#-l "CommonKmers, FOCUS, Metaphlan, MetaPhyler, Quikr, TIPP, mOTU, PhAT (U)" \
#-o opal_output_U_new \
#-d "2nd CAMI Challenge Mouse Gut Toy Dataset"
