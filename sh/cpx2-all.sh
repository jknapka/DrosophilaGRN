# Source the setup.sh script before running this.

# Run all the comparisons in a given directory.
# arg[0] is the name of the directory under work/
# that contains the # trajectory files, such as "R1"
# if the directory is work/R1. Each trajectory file whose
# name fits the pattern FBgbX-FBgnY-arg[0].trj will
# be compared with the dual trajectory file named
# FBgnX-FBgnY-__X__arg[0].trj. When CPX2 outputs
# "interesting" (non-conserved) results, those
# results are written to out/arg[0].out

TDIR_LEAF=$1
TDIR=${PROJ_DIR}/work/${TDIR_LEAF}

PMODE=3
if [ ! "" == "$2" ] ; then PMODE=$2 ; fi

grepElement() {
    local e
    for e in "${@:2}"; do 
        echo $e | grep -q $1 && return 0 
    done
    return 1
}


nn=0
declare -a T1FILES
export T1FILES
while read line
do
    T1FILES[$nn]=${TDIR}/"$line"
    (( ++nn ))
done < <(ls -1 ${TDIR} | grep -v __X__)

nn=0
declare -a T2FILES
export T2FILES
while read line
do
    T2FILES[$nn]="$line"
    (( ++nn ))
done < <(ls -1 ${TDIR}/*__X__*)


nT1FILES=${#T1FILES[*]}
nT2FILES=${#T2FILES[@]}

if [ ! $nT1FILES -eq $nT2FILES ] ; then
    echo ERROR: unexpected number of input files in each condition: $nT1FILES vs $nT2FILES
    exit 1
fi

CPX2_NOTHING=10
OUT_FILE=${PROJ_DIR}/out/${TDIR_LEAF}.out
if [ -e ${OUT_FILE} ] ; then rm -f ${OUT_FILE} ; fi

echo Relation,FBidA,FBidB,Celltype,pd,pc,pt,pz\(1\),pz\(2\)
for ((ff=0; $ff<${nT1FILES}; ++ff)) ; do
    T2FILES[$ff]=${TDIR}/$(basename ${T1FILES[$ff]} "${TDIR_LEAF}.trj")__X__${TDIR_LEAF}.trj
    #echo $ff : Comparing ${T1FILES[$ff]} vs ${T2FILES[$ff]}
    declare -a CPX2_OUT
    CPX2_OUT=()
    nn=0

    # Collect glnsp output into a temp file.
    TEMP_FILE=$(tempfile -d .)
    ${PROJ_DIR}/bin/glnsp -M comparison -P $PMODE -p 1 -g 1 -K 0 -J 0   -1 ${T1FILES[$ff]} -2 ${T2FILES[$ff]} > ${TEMP_FILE} 2>&1

    # If there is more output than the "no significant results" boilerplate,
    # grep the array for the relationships of interest.
    if [  $(wc -l ${TEMP_FILE} | awk '{print $1}') -gt ${CPX2_NOTHING} ] ; then
        FBIDS="$(basename ${T1FILES[$ff]} .trj | tr - , )"
        PVALS=($(egrep -o "p[a-z]+=[0-9.e-]+" ${TEMP_FILE} | head -5 | egrep -o "[-.0-9e]+" | tr \  ,))
        PVALS=$(echo ${PVALS[@]} | tr \  ,)
        if grep -q CONSERVED ${TEMP_FILE} ; then
            echo CONSERVED,$FBIDS,$PVALS
            echo "############################################" >> $OUT_FILE
            echo  @@@@ $(basename ${T1FILES[$ff]})   VS  $(basename ${T2FILES[$ff]}) >> $OUT_FILE
            echo "!!!CONSERVED" >> $OUT_FILE
            echo "############################################" >> $OUT_FILE
            echo "" >> $OUT_FILE
            echo "" >> $OUT_FILE
        else
            REL_TYPE=RELATIVE
            if grep -q ABSOLUTE ${TEMP_FILE} ; then
                REL_TYPE=ABSOLUTE
            fi
            echo ${REL_TYPE}_DIFFERENTIAL,$FBIDS,$PVALS
            echo "############################################" >> $OUT_FILE
            echo  @@@@ $(basename ${T1FILES[$ff]})   VS  $(basename ${T2FILES[$ff]}) >> $OUT_FILE
            echo cat ${TEMP_FILE} >> $OUT_FILE
            echo "############################################" >> $OUT_FILE
            echo "" >> $OUT_FILE
            echo "" >> $OUT_FILE
        fi
    fi
    rm $TEMP_FILE
done

