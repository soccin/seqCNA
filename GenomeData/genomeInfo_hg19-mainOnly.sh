# Need to use the version which has lower case ==> upper case
# to fix issue with VEP

echo "***** WARNING *****"
echo "* Should only be used for FACETS or other methods"
echo "* that do not care about exact genome match."
echo "* hg19-mainOnly only has main chromosomes but we"
echo "* are returning links for full build here."
echo "*****"

GENOME=/ifs/work/socci/Depot/Genomes/H.sapiens/hg19_Fixed/hg19_Fixed.fasta
GENOMEFAI=/ifs/work/socci/Depot/Genomes/H.sapiens/hg19_Fixed/hg19_Fixed.fasta.fai
EXACDB=/ifs/work/socci/Depot/Pipelines/Variant/PostProcess/db/ExAC.r0.3.sites.pass.minus_somatic.vcf.gz
FACETSNPS=/home/socci/Code/Pipelines/FACETS/FACETS.app/lib/dbsnp_137.hg19__RmDupsClean__plusPseudo50__DROP_SORT.vcf.gz