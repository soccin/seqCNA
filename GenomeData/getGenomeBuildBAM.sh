#!/bin/bash

if [ "$#" != "1" ]; then
    echo usage getGenomeBuild.sh BAM
    exit
fi

SAMTOOLS=/mnt/ceph/users/info/software/samtools-1.10/samtools

GENOME_MD5=$($SAMTOOLS view -H $1 | egrep "^@SQ" | cut -f-3 | sort  | md5sum - | awk '{print $1}')

case $GENOME_MD5 in
    b879c678e7fd80718cf20d10c6b846e4)
    # b37 gatk /ifs/depot/assemblies/H.sapiens/b37/b37.dict
    echo "b37"
    ;;

    117fce86b797081e0af6d69cbd94dcde)
    # b37 version used by DMP pipeline
    echo "b37_dmp"
    ;;

    5b4e380a6b4fc3494cfc66c917d41b37)
    # UCSC hg19 /ifs/depot/assemblies/H.sapiens/hg19/hg19.dict
    echo "hg19"
    ;;

    5322924312cb71bfc76e9018d9313676)
    # Xenograph genome (b37+mm10)
    echo "hg19+mm10"
    ;;

    3d72c6961689390556ed2d5a33e66e17)
    # Main chromosomes only (used by cfDNA collaboration)
    echo "hg19-mainOnly"
    ;;

    933b376d936c265fc6b44c8bd19fc66d)
    # TCGA BAMs UR:ftp://ftp.ncbi.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/special_requests/GRCh37-lite.fa.gz
    # AS:GRCh37-lite (b37-ish)
    echo "GRCh37-lite"
    ;;

    7f8c5323ff7e0ff6b5d40efe53eaf050)
    # BIC Xeno-graph genome
    echo "b37+mm10"
    ;;

    213e3e9448956f0de53bf8e9b3d6b7c1)
    # Human b38 (hg38)
    echo "b38"
    ;;

    d660fd17a979374182d3ba8b6d76cac0)
    # UCSC mm10 /ifs/depot/assemblies/M.musculus/mm10/mm10.dict
    echo "mm10"
    ;;

    34839afd79d8b772037ed7a0e0a4f9c3)
    # UCSC mm10
    echo "mm10_hBRAF_V600E"
    ;;

    f9cd233a3d5c9540eece434c65f84f1c)
    # mm9 Full
    echo "mm9Full"
    ;;

    0835b244e0adb20253e5fa1c4ee58ec4)
    # mouse_GRCm38
    echo "GRC_m38"
    ;;

    # Transgenic model_MellingI_1911b
    6e3786e28e9bb34a0c760142e32f835f)
    echo "mm10-MellingI_1911b"
    ;;

    # Relabling of mm10 chr1=>1

    # WT MM10
    5328e9d5744c2ad31c3c442028520510)
    echo "mm10-relabel"
    ;;

    # Transgenic cmo_hBRAF_Ex15-18
    6a3985296926fb48fd65ba5feae4e640)
    echo "mm10-relabel"
    ;;

    # Transgenic model_MellingI_1911b
    ec15a0439c7d7d47e079fbb4c4158267)
    echo "mm10-relabel"
    ;;


    *)
    echo "unknown" $GENOME_MD5
    ;;
esac

