# Source the setup.sh script before running this.

# Run 4-way comparisons for all interactions. The
# trajectories for each interaction are in the
# work/<celltype> directories.

TDIR=${PROJ_DIR}/work

GENE_INTER_FILE=data/interactions-present-by-gene.txt

CELL_TYPES="EC EE EB ISC"

grepElement() {
    local e
    for e in "${@:2}"; do 
        echo $e | grep -q $1 && return 0 
    done
    return 1
}


nn=0
declare -a I1GENES
while read line
do
    I1GENES[$nn]=$line
    (( ++nn ))
done < <(cut -d, -f1 < ${GENE_INTER_FILE})

nn=0
declare -a I2GENES
while read line
do
    I2GENES[$nn]=$line
    (( ++nn ))
done < <(cut -d, -f2 < ${GENE_INTER_FILE})

nI1GENES=${#I1GENES[*]}
nI2GENES=${#I2GENES[*]}

if [ ! $nI1GENES -eq $nI2GENES ] ; then
    echo ERROR: unexpected number of genes in all interactions: $nI1GENES vs $nI2GENES
    exit 1
fi

CPX2_NOTHING=10
OUT_FILE=${PROJ_DIR}/out/4way-by-cell-type.out
if [ -e ${OUT_FILE} ] ; then rm -f ${OUT_FILE} ; fi

echo Relation4,FBidA,FBidB
for ((gg=0; $gg<${nI1GENES}; ++gg)) ; do
    G1=${I1GENES[$gg]}
    G2=${I2GENES[$gg]}

    # We will used the trajectory-list file option of CPX2.
    echo -n "" > tflist.txt

    for CT in $CELL_TYPES ; do
        echo work/${CT}/${G1}-${G2}-${CT}.trj >> tflist.txt
    done

    echo "PERFORMING ANALYSIS ON TRAJECTORIES:"
    cat tflist.txt

    declare -a CPX2_OUT
    CPX2_OUT=()
    nn=0
    while read line ; do
        CPX2_OUT[$nn]="$line"
        (( ++nn ))
    #done < <(echo Nothing to see here.)
    done < <(${PROJ_DIR}/bin/glnsp -M comparison  -p 1 -g 1 -K 0 -J 0 -T tflist.txt)
    if [ ${#CPX2_OUT[*]} -gt ${CPX2_NOTHING} ] ; then
        FBIDS="${G1},${G2}"
        if grepElement CONSERVED "$(echo ${CPX2_OUT[@]})" ; then
            echo "############################################" >> $OUT_FILE
            echo CONSERVED,$FBIDS
            echo  @@@@ $(basename ${T1FILES[$ff]})   VS  $(basename ${T2FILES[$ff]}) >> $OUT_FILE
            echo "!!!CONSERVED" >> $OUT_FILE
            echo "############################################" >> $OUT_FILE
            echo "" >> $OUT_FILE
            echo "" >> $OUT_FILE
        else
            REL_TYPE=RELATIVE
            if grepElement ABSOLUTE "$(echo ${CPX2_OUT[@]})" ; then
                REL_TYPE=ABSOLUTE
            fi
            echo ${REL_TYPE}_DIFFERENTIAL,$FBIDS
            echo "############################################" >> $OUT_FILE
            echo  @@@@ $(basename ${T1FILES[$ff]})   VS  $(basename ${T2FILES[$ff]}) >> $OUT_FILE
            for (( ii=0; $ii<${#CPX2_OUT[*]}; ++ii )) ; do
                echo ${CPX2_OUT[$ii]} >> $OUT_FILE
            done
            echo "############################################" >> $OUT_FILE
            echo "" >> $OUT_FILE
            echo "" >> $OUT_FILE
        fi
    fi
done

