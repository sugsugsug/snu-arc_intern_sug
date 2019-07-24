## Copyright Broad Institute, 2019
## 
## This WDL pipeline implements data pre-processing according to the GATK Best Practices.  
##
## Requirements/expectations :
## - Pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile 
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
##
## Output :
## - A clean BAM file and its index, suitable for variant discovery analyses.
##
## Software version requirements 
## - GATK 4 or later
## - BWA 0.7.15-r1140
## - Picard 2.16.0-SNAPSHOT
## - Samtools 1.3.1 (using htslib 1.3.1)
## - Python 2.7
##
## Cromwell version support 
## - Successfully tested on v37
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.

# WORKFLOW DEFINITION 
workflow PreProcessingForVariantDiscovery_GATK4 {

  String sample_name
  String ref_name

  File flowcell_unmapped_bams_list
  String unmapped_bam_suffix
  
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices

  String? bwa_commandline_override
  String bwa_commandline = select_first([bwa_commandline_override, "bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"]) 
  Int compression_level
  
  String? gatk_docker_override
  String gatk_docker = select_first([gatk_docker_override, "broadinstitute/gatk:4.1.0.0"])
  String? gatk_path_override
  String gatk_path = select_first([gatk_path_override, "/arc-local/gatk/Profiling__/gatk/build/bundle-files-collected/gatk"])

  String? gotc_docker_override
  String gotc_docker = select_first([gotc_docker_override, "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"])
  String? gotc_path_override
  String gotc_path = select_first([gotc_path_override, "/usr/gitc/"])

  String bwa_path = "/arc-local/gatk/Profiling__/bwa/"

  String picard_path = "/arc-local/gatk/Profiling__/picard/build/libs/"

  String? python_docker_override
  String python_docker = select_first([python_docker_override, "python:2.7"])  

  Int flowcell_small_disk
  Int flowcell_medium_disk
  Int agg_small_disk
  Int agg_medium_disk
  Int agg_large_disk

  String? preemptible_tries_override
  Int preemptible_tries = select_first([preemptible_tries_override, "3"])

  String base_file_name = sample_name + "." + ref_name

  Array[File] flowcell_unmapped_bams = read_lines(flowcell_unmapped_bams_list)

  # Get the version of BWA to include in the PG record in the header of the BAM produced 
  # by MergeBamAlignment. 
  Array[File] MergeBamAlignment_output_bam = read_lines("/arc-local/gatk/inputs/ForChip/markduplicates_spark_sort_fix.txt")

# Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  call MarkDuplicates {
    input:
      TI = 0,
      input_bams = MergeBamAlignment_output_bam,
      output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = base_file_name + ".duplicate_metrics",
      docker_image = gatk_docker,
      gatk_path = gatk_path,
        gotc_path = picard_path,
      disk_size = agg_large_disk,
      compression_level = compression_level,
      preemptible_tries = preemptible_tries
  }

  call TS as S2{
	input:
		St = "After MarkDuplicates:\n",
		I = MarkDuplicates.T
  }
  # Sort aggregated+deduped BAM file and fix tags
  call SortAndFixTags {
    input:
      TI = S2.TSR,
      input_bam = MarkDuplicates.output_bam,
      output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      docker_image = gatk_docker,
      gatk_path = gatk_path,
        gotc_path = picard_path,
      disk_size = agg_large_disk,
      preemptible_tries = 0,
      compression_level = compression_level
  }
  call TS as S3{
	input:
		St = "After SortAndFixTags:\n",
		I = SortAndFixTags.T
  }

  # Sort aggregated+deduped BAM file and fix tags
  output {
    File duplication_metrics = MarkDuplicates.duplicate_metrics
  } 
}

# TASK DEFINITIONS

# Get version of BWA
task GetBwaVersion {
  
  Int TI  

  Int preemptible_tries
  String mem_size

  String docker_image
  String bwa_path

  command {
    # Not setting "set -o pipefail" here because /bwa has a rc=1 and we don't want to allow rc=1 to succeed 
    # because the sed may also fail with that error and that is something we actually want to fail on.
    ${bwa_path}bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
  runtime {
    preemptible: preemptible_tries
    memory: mem_size
  }
  output {
    String version = read_string(stdout())
    Int T = 0
  }
}

# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment
task SamToFastqAndBwaMem {
  Int TI

  File input_bam
  String bwa_commandline
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit), 
  # listing the reference contigs that are "alternative". Leave blank in JSON for legacy 
  # references such as b37 and hg19.
  File? ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa

  Int compression_level
  Int preemptible_tries
  Int disk_size
  String mem_size
  String num_cpu

  String docker_image
  String bwa_path
  String gotc_path
  String java_opt

  command <<<
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}

		java -Dsamjdk.compression_level=${compression_level} ${java_opt} -jar ${gotc_path}picard.jar \
      SamToFastq \
			INPUT=${input_bam} \
			FASTQ=/dev/stdout \
			INTERLEAVE=true \
			NON_PF=true \
    | \
		${bwa_path}${bwa_commandline} /dev/stdin -  2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) \
    > tmp.txt
		(time (cat tmp.txt | (samtools view -1 - > ${output_bam_basename}.bam))) 2> "/arc-local/gatk/Profiling__/SamtoolsTime.txt"

  >>>
  runtime {
    preemptible: preemptible_tries
    memory: mem_size
    cpu: num_cpu
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    Int T = 0
    File output_bam = "${output_bam_basename}.bam"
    File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
  }
}

# Merge original input uBAM file with BWA-aligned BAM file
task MergeBamAlignment {
  Int TI

  File unmapped_bam
  String bwa_commandline
  String bwa_version
  File aligned_bam
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  Int compression_level
  Int preemptible_tries
  Int disk_size
  String mem_size

  String docker_image
  String gatk_path
  String gotc_path
  String java_opt

  command {
    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    java -Dsamjdk.compression_level=${compression_level} ${java_opt} -jar ${gotc_path}picard.jar\
      MergeBamAlignment \
      VALIDATION_STRINGENCY=SILENT \
      EXPECTED_ORIENTATIONS=FR \
      ATTRIBUTES_TO_RETAIN=X0 \
      ALIGNED_BAM=${aligned_bam} \
      UNMAPPED_BAM=${unmapped_bam} \
      OUTPUT=${output_bam_basename}.bam \
      REFERENCE_SEQUENCE=${ref_fasta} \
      PAIRED_RUN=true \
      SORT_ORDER="unsorted" \
      IS_BISULFITE_SEQUENCE=false \
      ALIGNED_READS_ONLY=false \
      CLIP_ADAPTERS=false \
      MAX_RECORDS_IN_RAM=2000000 \
      ADD_MATE_CIGAR=true \
      MAX_INSERTIONS_OR_DELETIONS=-1 \
      PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
      PROGRAM_RECORD_ID="bwamem" \
      PROGRAM_GROUP_VERSION="${bwa_version}" \
      PROGRAM_GROUP_COMMAND_LINE="${bwa_commandline}" \
      PROGRAM_GROUP_NAME="bwamem" \
      UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
      ALIGNER_PROPER_PAIR_FLAGS=true \
      UNMAP_CONTAMINANT_READS=true
  }
  runtime {
    preemptible: preemptible_tries
    memory: mem_size
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    Int T = 0
    File output_bam = "${output_bam_basename}.bam"
  }
}

task Sort {
  Int TI

  File input_bam
  String output_bam_basename
  
  Int compression_level
  Int preemptible_tries
  Int disk_size
  String mem_size

  String gatk_path
  String gotc_path
  String java_opt_sort

  command {
    set -o pipefail

    java -Dsamjdk.compression_level=${compression_level} ${java_opt_sort} -jar ${gotc_path}picard.jar\
      SortSam \
      --INPUT ${input_bam} \
      --OUTPUT ${output_bam_basename}.bam \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false \
  }
  runtime {
    preemptible: preemptible_tries
    memory: mem_size
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    Int T = 0
    File output_bam = "${output_bam_basename}.bam"
  }
}
# Sort BAM file by coordinate order and fix tag values for NM and UQ
task SortAndFixTags {
  Int TI

  File input_bam
  String output_bam_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  
  Int compression_level
  Int preemptible_tries
  Int disk_size
  String mem_size

  String docker_image
  String gatk_path
  String gotc_path
  String java_opt_sort
  String java_opt_fix

  command {
    set -o pipefail

	sudo time -ap -o "/arc-local/gatk/LocalTool/surfaceTime/sortspark.txt" \
    ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt_sort}" \
      SortSamSpark \
      -I ${input_bam} \
      -O tmp.txt \
      --conf 'spark.executor.cores=8' \
      --sort-order "coordinate" \
      --create-output-bam-index false 
	sudo time -ap -o "/arc-local/gatk/LocalTool/surfaceTime/fixtags.txt" \
    java -Dsamjdk.compression_level=${compression_level} ${java_opt_fix} -jar ${gotc_path}picard.jar \
      SetNmMdAndUqTags \
      I=tmp.txt \
      O=${output_bam_basename}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      R=${ref_fasta}
  }
  runtime {
    preemptible: preemptible_tries
    memory: mem_size
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    Int T = 0
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}

task MarkDuplicates {
  Int TI

  Array[File] input_bams
  String output_bam_basename
  String metrics_filename
  
  Int compression_level
  Int preemptible_tries
  Int disk_size
  String mem_size

  String docker_image
  String gatk_path
  String gotc_path
  String java_opt

 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
 # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
	sudo time -ap -o "/arc-local/gatk/LocalTool/surfaceTime/markduplicates.txt" \
    java -Dsamjdk.compression_level=${compression_level} ${java_opt} -jar ${gotc_path}picard.jar \
      MarkDuplicates \
      INPUT=${sep=' INPUT=' input_bams} \
      OUTPUT=${output_bam_basename}.bam \
      METRICS_FILE=${metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname" \
      CREATE_MD5_FILE=true
  }
  runtime {
    preemptible: preemptible_tries
    memory: mem_size
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    Int T = 0
    File output_bam = "${output_bam_basename}.bam"
    File duplicate_metrics = "${metrics_filename}"
  }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicatesSparkandFix {
  Int TI

  Array[File] input_bams
  File ref_fasta
  String output_bam_basename
  String metrics_filename
  
  Int compression_level
  Int preemptible_tries
  Int disk_size
  String mem_size

  String docker_image
  String gatk_path
  String gotc_path
  String java_opt
  String java_opt_fix
  String java_opt_sort

 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
 # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    set -o pipefail

    ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" \
      MarkDuplicatesSpark \
      -I ${sep=' -I ' input_bams} \
      -O /dev/stdout \
      --conf 'spark.executor.cores=8' \
      -M ${metrics_filename} \
      -VS SILENT \
      --optical-duplicate-pixel-distance 2500 \
      --treat-unsorted-as-querygroup-ordered true \
	| \
    java -Dsamjdk.compression_level=${compression_level} ${java_opt_fix} -jar ${gotc_path}picard.jar\
        SetNmMdAndUqTags \
	I=/dev/stdin \
	O=${output_bam_basename}.bam \
      R=${ref_fasta} \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true 
  }
  runtime {
    preemptible: preemptible_tries
    memory: mem_size
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    Int T = 0
    File output_bam = "${output_bam_basename}.bam"
    File duplicate_metrics = "${metrics_filename}"
  }
}

# Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {
  Int TI
  File ref_dict  
  
  Int preemptible_tries
  String mem_size

  String docker_image

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter. 
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("${ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in 
    # some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: mem_size
  }
  output {
    Int T = 0
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  Int TI

  File input_bam
  File input_bam_index
  String recalibration_report_filename
  Array[String] sequence_group_interval
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  
  Int preemptible_tries
  Int disk_size
  String mem_size

  String docker_image
  String gatk_path
  String java_opt

  command { 
    ${gatk_path} --java-options "${java_opt}" \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --use-original-qualities \
      -O ${recalibration_report_filename} \
      --known-sites ${dbSNP_vcf} \
      --known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
      -L ${sep=" -L " sequence_group_interval}
  }
  runtime {
    preemptible: preemptible_tries
    memory: mem_size
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    Int T = 0
    File recalibration_report = "${recalibration_report_filename}"
  }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
# Note that when run from GATK 3.x the tool is not a walker and is invoked differently.
task GatherBqsrReports {
  Int TI

  Array[File] input_bqsr_reports
  String output_report_filename

  Int preemptible_tries
  Int disk_size
  String mem_size

  String docker_image
  String gatk_path
  String java_opt

  command {
    ${gatk_path} --java-options "${java_opt}" \
      GatherBQSRReports \
      -I ${sep=' -I ' input_bqsr_reports} \
      -O ${output_report_filename}
  }
  runtime {
    preemptible: preemptible_tries
    memory: mem_size
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    Int T = 0
    File output_bqsr_report = "${output_report_filename}"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  Int TI

  File input_bam
  File input_bam_index
  String output_bam_basename
  File recalibration_report
  Array[String] sequence_group_interval
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  Int preemptible_tries
  Int disk_size 
  String mem_size

  String docker_image
  String gatk_path
  String java_opt

  command {  
    ${gatk_path} --java-options "${java_opt}" \
      ApplyBQSR \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${output_bam_basename}.bam \
      -L ${sep=" -L " sequence_group_interval} \
      -bqsr ${recalibration_report} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities
  }
  runtime {
    preemptible: preemptible_tries
    memory: mem_size
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    Int T = 0
    File recalibrated_bam = "${output_bam_basename}.bam"
  }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  Int TI

  Array[File] input_bams
  String output_bam_basename

  Int compression_level
  Int preemptible_tries
  Int disk_size
  String mem_size

  String docker_image
  String gatk_path
  String java_opt

  command {
    ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" \
      GatherBamFiles \
      --INPUT ${sep=' --INPUT ' input_bams} \
      --OUTPUT ${output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true
  }
  runtime {
    preemptible: preemptible_tries
    memory: mem_size
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    Int T = 0
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}

task TS {
	Int? I
	Array[Int]? AI
	String St
	command { echo "${St}" >> "/arc-local/gatk/TSresult.txt"
		date >> "/arc-local/gatk/TSresult.txt" }
	output {
		Int TSR = 0
	}
}


