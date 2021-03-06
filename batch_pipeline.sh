fastqDir=/data/fastq/tilapia/Huang/sex_difference/trimming
genomeDir=/data/tilapia_source_file
genomeFile=tilapia.fa
genomeName=${genomeFile%.*}
gffFile=tilapia_20171005_with_yp_genes.gff3
gffName=${gffFile%.*}
project_species=tilapia
R1_Append=_R1.gz
R2_Append=_R2.gz

cpu_threads=$(nproc)
gatk_path=$(echo $(which gatk)/$(readlink $(which gatk)) | sed 's/bin\/gatk\/\.\.\///g' | sed 's/gatk$//g')$(ls $(echo $(which gatk)/$(readlink $(which gatk)) | sed 's/bin\/gatk\/\.\.\///g' | sed 's/gatk$//g') | grep -P -o '.*.jar$')
picard_path=$(echo $(which picard)/$(readlink $(which picard)) | sed 's/bin\/picard\/\.\.\///g').jar
java_path=$(which java)
trimmomatic_path=$(echo $(which trimmomatic)/$(readlink $(which trimmomatic)) | sed 's/bin\/trimmomatic\/\.\.\///g' | sed 's/\/trimmomatic$//g')

cd ${fastqDir}
for ID in $(ls | grep -P -o '.*(?=(_R1.gz))')
do
	if [ ! -f ${ID}_R1_paired.fastq ] && [ ! -f ${ID}_sambamba.bam ];then
		echo process trimming ${ID}
		nohup trimmomatic PE -threads ${cpu_threads} ${ID}${R1_Append} ${ID}${R2_Append} ${ID}_R1_paired.fastq.gz ${ID}_R1_unpaired.fastq.gz ${ID}_R2_paired.fastq.gz ${ID}_R2_unpaired.fastq.gz ILLUMINACLIP:${trimmomatic_path}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 &
		pids="$pids $!"
		echo
		echo start waitting for finish ${ID} trimming
		wait $pids
	fi
done

echo ====== Finish Trimming ======
echo
echo ====== Start Mapping ======

for ID in $(ls | grep -P -o '.*(?=(_R1_paired.fastq.gz))')
do
	if [ ! -f ${ID}.sam ] && [ ! -f ${ID}_rg_added_sorted.bam ];then

		### Check if splicesites file exist or not ###
		if [ ! -f ${genomeDir}/${gffName}_splicesites.txt ]; then
			cd ${genomeDir}
			if [ ! ${gffName}.gtf ]; then
				gffread ${gffFile} -T -o ${gffName}.gtf
				pids="$pids $!"
				echo
				echo start waitting Transform GFF to GTF
				wait $pids
			fi

			hisat2_extract_splice_sites.py ${gffName}.gtf > ${gffName}_splicesites.txt
			pids="$pids $!"
			echo
			echo start waitting for creating Splicesites
			wait $pids

			cd ${fastqDir}
		fi

		### Create HISAT2 index ###
		if [ ! -f ${genomeDir}/${genomeFile}.1.ht2 ]; then
			cd ${genomeDir}
			hisat2-build -p ${cpu_threads} ${genomeFile} ${genomeFile}
			pids="$pids $!"
			echo
			echo start waitting for Building HISAT2 index
			wait $pids
			cd ${fastqDir}
		fi

		echo process mapping ${ID}
		nohup hisat2 -p ${cpu_threads} -x ${genomeDir}/${genomeFile} -1 ${ID}_R1_paired.fastq.gz -2 ${ID}_R2_paired.fastq.gz --known-splicesite-infile ${genomeDir}/${gffName}_splicesites.txt -S ${ID}.sam &
		pids="$pids $!"
		echo
		echo start waitting for finish ${ID} mapping
		wait $pids
		if [ -f ${ID}_R1.fastq ];then
			echo ====== Deleting Ori FASTQ File ======
			rm ${ID}${R1_Append} 
			rm ${ID}${R2_Append} 
			rm ${ID}_trim.log
		fi
	fi
	if [ ! -f ${ID}_sambamba.bam ];then
		echo process Sambamba ${ID}
		nohup sambamba view -S ${ID}.sam -t ${cpu_threads} -f bam -o ${ID}_sambamba.bam &
		pids="$pids $!"
		echo
		echo start waitting for finish ${ID} Sambamba Convert
		wait $pids

		echo process Sambamba Sort ${ID}
		nohup sambamba sort ${ID}_sambamba.bam -t ${cpu_threads} --tmpdir ./tmp &
		pids="$pids $!"
		echo
		echo start waitting for finish ${ID} Sambamba Sort
		wait $pids
		if [ -f ${ID}_R1_paired.fastq.gz ];then
			echo ====== Deleting FASTQ File ======
		#	rm ${ID}_R1_paired.fastq
		#	rm ${ID}_R2_paired.fastq
			rm ${ID}_R1_unpaired.fastq.gz
			rm ${ID}_R2_unpaired.fastq.gz
		fi
	fi
done
for ID in $(ls | grep -P -o '.*(?=(_sambamba.sorted.bam))')
do
	if [ ! -d ${ID}_Ballgown_Table ];then
		echo process StringTie ${ID}
		nohup stringtie ${ID}_sambamba.sorted.bam -eG ${genomeDir}/${gffFile} -p ${cpu_threads} -b ${ID}_Ballgown_Table > ${ID}_stringtie.log &
		pids="$pids $!"
	fi
done
	echo
	echo start waitting for finish ${ID} StringTie
	wait $pids
	echo finish ${ID} StringTie

echo ====== Finish Sambamba-StringTie Pipeline ======
echo
echo ====== Start Variant Calling =======

for ID in $(ls | grep -P -o '.*(?=(.sam))')
do
	### Check fai index exist ###
	if [ ! -f ${genomeDir}/${genomeFile}.fai ];then
		cd ${genomeDir}/
		samtools faidx ${genomeFile}
		cd ${fastqDir}/
	fi
	### Check dict exist ###
	if [ ! -f ${genomeDir}/${genomeName}.dict ];then
		cd ${genomeDir}
		picard CreateSequenceDictionary R= ${genomeFile} O= ${genomeName}.dict
		cd ${fastqDir}
	fi
	if [ ! -f ${ID}_rg_added_sorted.bam ];then
		echo process AddOrReplaceReadGroups ${ID}
		RGID=$(grep -v @ ${ID}.sam | head -1 | cut -f 1 | cut -d ':' -f 3 | cut -c -5)
		RGPU=$(grep -v @ ${ID}.sam | head -1 | cut -f 1 | cut -d ':' -f 3)
		RGLB=$(echo ${ID} | cut -d "_" -f 1)
		nohup ${java_path} -Djava.io.tmpdir=/data/tmp -Xms1g -Xmx4g -jar ${picard_path} AddOrReplaceReadGroups I=${ID}.sam O=${ID}_rg_added_sorted.bam SO=coordinate RGID=${RGID} RGLB=${RGLB} RGPL=illumina RGPU=${RGPU} RGSM=${ID} &
		pids="$pids $!"
		echo
		echo start waitting for finish ${ID} AddOrReplaceReadGroups
		wait $pids
	fi
done

for ID in $(ls | grep -P -o '.*(?=(_rg_added_sorted.bam))')
do
	if [ ! -f ${ID}_dedupped.bam ];then
		echo process MarkDuplicates ${ID}
		nohup ${java_path} -Djava.io.tmpdir=/data/tmp -Xms1g -Xmx4g -jar ${picard_path} MarkDuplicates I=${ID}_rg_added_sorted.bam O=${ID}_dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics &
		pids="$pids $!"
		echo
		echo start waitting for finish ${ID} MarkDuplicates
		wait $pids
		if [ -f ${ID}.sam ];then
			echo ====== Deleting SAM File ======
			rm ${ID}.sam
		fi
	fi
done

for ID in $(ls | grep -P -o '.*(?=(_dedupped.bam))')
do
	if [ ! -f ${ID}_split.bam ];then
		echo process SplitNCigarReads ${ID}
		nohup gatk --java-options "-Djava.io.tmpdir=/data/tmp -Xmx20g" SplitNCigarReads -R ${genomeDir}/${genomeFile} -I ${ID}_dedupped.bam -O ${ID}_split.bam &
		pids="$pids $!"
		echo
		echo start waitting for finish ${ID} SplitNCigarReads
		wait $pids
	fi
done

for ID in $(ls | grep -P -o '.*(?=(_split.bam))')
do
	if [ ! -f ${ID}.vcf ];then
		echo process Variant Calling ${ID}
		nohup gatk --java-options "-Djava.io.tmpdir=/data/tmp -Xms1g -Xmx20g" HaplotypeCaller -R ${genomeDir}/${genomeFile} -I ${ID}_split.bam -stand-call-conf 20.0 -O ${ID}.vcf &
		pids="$pids $!"
		echo
		echo start waitting for finish ${ID} Variant Calling
		wait $pids
	fi
done
echo ====== Finish Variant Calling ======
echo
echo ====== Start VEP Discovering ======

for ID in $(ls | grep -P -o '.*(?=(.vcf))')
do
	if [ ! -f ${ID}_variant_effect.txt ];then
		### Check GFF index exist ###
		if [ ! -f ${genomeDir}/${gffName}.gff.gz.tbi ];then
			cd ${genomeDir}
			### Check GFF gz file exist ###
			if [ ! -f ${genomeDir}/${gffName}.gff.gz ];then
				grep -v "#" ${genomeDir}/${gffFile} | sort -k1,1 -k4,4n -k5,5n | bgzip -c >${genomeDir}/${gffName}.gff.gz
				pids="$pids $!"
				wait $pids
			fi
			tabix -p gff ${genomeDir}/${gffName}.gff.gz
			pids="$pids $!"
			wait $pids
			cd ${fastqDir}
		fi
		### Check fasta index exist ###
		if [ ! -f ${genomeDir}/${genomeName}.gz ];then
			cd ${genomeDir}
			bgzip ${genomeDir}/${genomeFile}
			pids="$pids $!"
			wait $pids
			cd ${fastqDir}
		fi
		echo process Variant Discovering ${ID}
		vep -i ${ID}.vcf -gff ${genomeDir}/${gffName}.gff.gz -fasta ${genomeDir}/${genomeFile}.gz -o ${ID}_variant_effect.txt
		echo
		echo start waitting for finish Variant ${ID} Discovering
		pids="$pids $!"
		wait $pids
		echo
	fi
done

echo ====== Finish VEP Discovering ======
