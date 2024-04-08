#!/bin/bash


# Check if the input file argument is provided
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <inputdir> <outputfile>"
    exit 1
fi

workdir="$1"
runmode="$2"
#workdir="/home/hyunwoo/data_set/project/svProject/wgs/ceph/sv_call/"

step1="/home/hyunwoo/src/python/svProject/dnsv_workflow/step1-1.filter_contig_50bp_lcr_intraBND.py"
step2="/home/hyunwoo/src/python/svProject/dnsv_workflow/step2-1.get_roi_ref_kmer.py"
step3="/home/hyunwoo/src/python/svProject/dnsv_workflow/step2-2.get_roi_repeatkmer.py"
step4="/home/hyunwoo//src/python/svProject/dnsv_workflow/step3.estimate_vaf.py"

genome="/home/hyunwoo/data_set/genome/hg38/gatk_bundle/Homo_sapiens_assembly38.fasta"
lcrbed="/home/hyunwoo/data_set/project/svProject/genomic_context/low_complexity/LCR-hs38.bed"
repeatkmerdb="/home/hyunwoo/data_set/project/svProject/genomic_context/repeatMasker/Homo_sapiens_assembly38_repeatmaskerkmers"
regressor="/home/hyunwoo/data_set/project/svProject/wgs/simulation/rf_train_seed240312/rf_model_gridCV/randomForestRegressor_bestParameter.joblib"

cd $workdir

# Check if the input directory exists
if [ ! -d "$workdir" ]; then
    echo "Error: Input directory not found."
    exit 1
fi

# Iterate over subdirectories under the input directory

for family in "$workdir"/*/; do
    # Check if the path is a directory
    if [ -d "$family" ]; then
        echo "Processing subdirectory: $family"
        targetdir="${family}/etching_call_v2/"
        for subject in "$targetdir"/*/; do
            echo "$family : $subject"
            if [ -d "$subject" ]; then
                cd $subject
                if [ ! -d "vaf_filter" ]; then
                    mkdir "vaf_filter"
                fi
                cd "vaf_filter"
                outputdir="${subject}/vaf_filter/"

                familynum=$(basename ${family})
                subjectnum=$(basename ${subject})
                targetfile="${subject}/${familynum}_${subjectnum}.scored.filtered.typed.vcf"
                bam="${subject}/${familynum}_${subjectnum}.sort.bam"
                kmer="${subject}/${familynum}_${subjectnum}.sample"
                outputfile="${subject}/${familynum}_${subjectnum}.scored.filtered.typed.contigFilt_50bp_PASS_nonLCR_nonIntraBND.vcf"
                step1cmd="python ${step1} -v ${targetfile} -b ${lcrbed} -o ${outputfile}"

                echo "VCF FILTERING | INPUT FILE : " $targetfile
                echo $step1cmd
                eval $step1cmd
                echo $(basename $targetfile) STEP1 DONE

                inputfile=$outputfile
                output1="${familynum}_${subjectnum}_ROI_refkmer.fasta"
                output2="${familynum}_${subjectnum}_ROI_refkmer.pickle"
                step2cmd="python ${step2} -v ${inputfile} -g ${genome} -f ${output1} -k ${kmer} -o ${outputdir} -p ${output2}"
                echo "ROI REFERENCE KMER EXTRACTION | INPUT FILE : " $inputfile
                echo $step2cmd
                eval $step2cmd #extract ROI kmer
                echo $(basename $inputfile) STEP2-ROI REF KMER EXTRACTION DONE

                step3cmd="python ${step3} -i ${outputdir}/roi_refkmer_intersect -r ${repeatkmerdb} -o ${outputdir}"
                echo "ROI REPEAT KMER EXTRACTION | INPUT FILE : " $inputfile
                echo $step3cmd
                eval $step3cmd #extract ROI repeat kmer
                echo $(basename $inputfile) STEP3-ROI REPEAT KMER EXTRACTION DONE

                refpickle=$output2
                repeatpickle="${outputdir}/roi_repeatkmer_intersect.pickle"
                outputvaf="${outputdir}/${familynum}_${subjectnum}_eVAF_table.txt"
                step4cmd="python ${step4} -i ${bam} -v ${inputfile} -g ${genome} -m ${runmode} -k ${refpickle} -R ${repeatpickle} -r ${regressor} -o ${outputvaf}"
                echo "eVAF CALCULATION | INPUT FILE : " $inputfile
                echo $step4cmd
                eval $step4cmd #calculate eVAF
                echo $(basename $inputfile) STEP4 DONE
            fi
        done
    fi
done
