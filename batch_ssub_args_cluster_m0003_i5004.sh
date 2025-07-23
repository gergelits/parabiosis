#!/bin/sh

# set "TEST=1" to run a test
TEST=0

# load modules if you are on cluster
if [ "$(whoami)" == "vaclavgergelits" ]
then
    echo "Job run locally on Desktop."
else
    module load R/4.4.0
    module load ssub
fi

R_SCRIPT_NAME="batch_args_cluster.r"


CELLTYPE_ALL_TMP=( "Treg" "Tconv" "CD8" "NK" "B" )

TISSUE_ALL_TMP=( "Adrenals" "BoneMarrow" "Brain" "FRT" "IEL" "Kidney" "Liver" "LPL" \
                 "LN" "Lung" "MLN" "Muscle" "Pancreas" "PP" "Skin" "Spleen" "WAT" "Thymus" )

ARG_MCMC="i5004"

MODEL_ID="m0003"
CORES=$((${ARG_MCMC:1:10} % 100))
echo "$CORES"

# TEST DATA:
CELLTYPE_TEST=( "Tconv" )
TISSUE_TEST=( "Brain" )


if [ $TEST -eq 1 ]
then
    echo "TEST EXAMPLE"
    CELLTYPE_ALL=$CELLTYPE_TEST
    TISSUE_ALL=$TISSUE_TEST
else
    CELLTYPE_ALL=${CELLTYPE_ALL_TMP[*]}
    TISSUE_ALL=${TISSUE_ALL_TMP[*]}
fi

echo "${CELLTYPE_ALL[*]}"
echo "${TISSUE_ALL[*]}"

# calculate parabiosis:
for TIS_0 in "NULL"
do
    for ARG_CELLTYPE in ${CELLTYPE_ALL[@]}
    do
        for ARG_TISSUE in ${TISSUE_ALL[@]}
        do
            if [ "$(whoami)" == "vaclavgergelits" ]
            then
                Rscript --quiet --vanilla $R_SCRIPT_NAME \
                $ARG_CELLTYPE $ARG_TISSUE $TIS_0 $ARG_MCMC $MODEL_ID &
            else
                ssub -o "${ARG_CELLTYPE}_${ARG_TISSUE}_${TIS_0}_${ARG_MCMC}_${MODEL_ID}.log" --mem=20G -c $CORES --email \
                Rscript --quiet --vanilla $R_SCRIPT_NAME \
                $ARG_CELLTYPE $ARG_TISSUE $TIS_0 $ARG_MCMC $MODEL_ID &
            fi
        done
    done
done
# 
# 
# # OR SIMPLY:
# 
# # Rscript --quiet --vanilla "batch_args_cluster.r" "Treg" "Pancreas" "NULL" &

# # in PowerShell:
# $TISSUES = @( 'Adrenals', 'BoneMarrow', 'Brain', 'IEL', 'Kidney', 'Liver',
#               'LPL', 'LN', 'Lung', 'Spleen',
#               'MLN', 'Muscle', 'Pancreas', 'PP', 'Skin', 'WAT' )
# foreach ( $tis_i in $TISSUES )
# {
#     Rscript --quiet --vanilla "batch_args_cluster.r" "NK" $tis_i "Spleen" &
# }
