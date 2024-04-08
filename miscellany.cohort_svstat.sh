#!/bin/bash


# Check if the input file argument is provided
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <inputdir> <outputfile>"
    exit 1
fi

workdir="$1"
#workdir="/home/hyunwoo/data_set/project/svProject/wgs/ceph/sv_call/"
filtersrc="/home/hyunwoo/src/python/svProject/dnsv_workflow/step1-1.filter_contig_50bp_lcr_intraBND.py"
lcrbed="/home/hyunwoo/data_set/project/svProject/genomic_context/low_complexity/LCR-hs38.bed"

cd $workdir

# Check if the input directory exists
if [ ! -d "$workdir" ]; then
    echo "Error: Input directory not found."
    exit 1
fi

# Iterate over subdirectories under the input directory
statfile="$2"

if [ ! -f "$statfile" ]; then
    touch $statfile
else
    rm $statfile
    touch $statfile
fi


printf "sample\tfamily\tDEL\tDUP\tINV\tBND\n" >> $statfile

for family in "$workdir"/*/; do
    # Check if the path is a directory
    if [ -d "$family" ]; then
        echo "Processing subdirectory: $family"
        targetdir="${family}/etching_call_v2/"
        for subject in "$targetdir"/*/; do
            #echo "$family : $subject"
            if [ -d "$subject" ]; then
                familynum=$(basename ${family})
                subjectnum=$(basename ${subject})
                #targetfile="${subject}/${familynum}_${subjectnum}.scored.filtered.typed.vcf"
                targetfile="${subject}/${familynum}_${subjectnum}.scored.filtered.typed.contigFilt_50bp_PASS_nonLCR_nonIntraBND.addVAF.dnsv.vcf"
                #outputfile="${subject}/${familynum}_${subjectnum}.scored.filtered.typed.contigFilt_50bp_PASS_nonLCR_nonIntraBND.vcf"
                filtercmd="python ${filtersrc} -v ${targetfile} -b ${lcrbed} -o ${outputfile}"

                echo "$targetfile"
                #echo $filtercmd
                #eval $filtercmd
                if [ -f "$targetfile" ]; then
                    del=$(grep -v '^#' ${targetfile} | grep 'SVTYPE=DEL' -c)
                    dup=$(grep -v '^#' ${targetfile} | grep 'SVTYPE=DUP' -c)
                    inv=$(grep -v '^#' ${targetfile} | grep 'SVTYPE=INV' -c)
                    bnd=$(grep -v '^#' ${targetfile} | grep 'SVTYPE=BND' -c)
                    printf "${subjectnum}\t${familynum}\t${del}\t${dup}\t${inv}\t${bnd}\n" >> $statfile
                fi
            fi
        done
        # Add your commands here to process each subdirectory
    fi
done
