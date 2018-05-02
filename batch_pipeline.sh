fastqDir=/data/fastq/tilapia/Huang/sex_difference
genomeDir=/data/tilapia_source_file
genomeFile=tilapia.fa
genomeName=${genomeFile%.*}
gffFile=tilapia_20171005_with_yp_genes.gff3
gffName=${gffFile%.*}
project_species=tilapia
R1_Append=_R1.fastq
R2_Append=_R2.fastq
cpu_threads=grep -c ^processor /proc/cpuinfo
cd ${fastqDir}
for ID in $(ls | grep -P -o '.*(?=(_R1.fastq))')
do
	if [ ! -f ${ID}_R1_paired.fastq ] && [ ! -f ${ID}_sambamba.bam ];then
		echo process trimming ${ID}
		nohup trimmomatic PE -threads ${cpu_threads} -trimlog ${ID}_trim.log ${ID}${R1_Append} ${ID}${R2_Append} ${ID}_R1_paired.fastq ${ID}_R1_unpaired.fastq ${ID}_R2_paired.fastq ${ID}_R2_unpaired.fastq ILLUMINACLIP:/opt/conda/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &
		pids="$pids $!"
		echo
		echo start waitting for finish ${ID} trimming
		wait $pids
	fi
done

echo ====== Finish Trimming ======
echo
echo ====== Start Mapping ======

for ID in $(ls | grep -P -o '.*(?=(_R1_paired.fastq))')
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
		nohup hisat2 -p ${cpu_threads} -x ${genomeDir}/${genomeFile} -1 ${ID}_R1_paired.fastq -2 ${ID}_R2_paired.fastq --known-splicesite-infile ${genomeDir}/${gffName}_splicesites.txt -S ${ID}.sam &
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
		if [ -f ${ID}_R1_paired.fastq ];then
			echo ====== Deleting FASTQ File ======
		#	rm ${ID}_R1_paired.fastq
		#	rm ${ID}_R2_paired.fastq
			rm ${ID}_R1_unpaired.fastq
			rm ${ID}_R2_unpaired.fastq
		fi
	fi
one
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
		nohup /opt/conda/bin/java -Djava.io.tmpdir=/data/tmp -Xms1g -Xmx4g -jar /opt/conda/share/gatk4-4.0.4.0-0/gatk-package-4.0.4.0-local.jar AddOrReplaceReadGroups I=${ID}.sam O=${ID}_rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample >${ID}_addorresortgroup.log &
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
		nohup /opt/conda/bin/java -Djava.io.tmpdir=/data/tmp -Xms1g -Xmx4g -jar /opt/conda/share/gatk4-4.0.4.0-0/gatk-package-4.0.4.0-local.jar MarkDuplicates I=${ID}_rg_added_sorted.bam O=${ID}_dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics > ${ID}_markduplicates.log &
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
		nohup java -Djava.io.tmpdir=/data/tmp -Xmx20g -jar /opt/conda/share/gatk4-4.0.4.0-0/gatk-package-4.0.4.0-local.jar -T SplitNCigarReads -R ${genomeDir}/${genomeFile} -I ${ID}_dedupped.bam -o ${ID}_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS &
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
		nohup java -Djava.io.tmpdir=/data/tmp -Xmx60g -jar /opt/conda/share/gatk4-4.0.4.0-0/gatk-package-4.0.4.0-local.jar HaplotypeCaller -R ${genomeDir}/${genomeFile} -I ${ID}_split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -O ${ID}.vcf &
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
