version 1.0

task bwa_mem {
  input {
    String sample_name
    File fastq1_gz
    File fastq2_gz
    File genome_index_tgz
    Int min_seed_length = 19
    String? read_group
    String docker_image = "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    Int cpu = 4
    Int memory_gb = 8
    Int? disk_gb
  }

  String genome_index_basename = basename(genome_index_tgz, ".tar.gz")
  String actual_read_group = select_first([
    read_group,
    "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tLB:${sample_name}\\tPL:ILLUMINA"
 ])
  Int actual_disk_gb = select_first([
    disk_gb,
    ceil(2 * (size(genome_index_tgz, "G") + size(fastq1_gz, "G") + size(fastq2_gz, "G")))
  ])

  command <<<
    set -eux
    tar xzvf ~{genome_index_tgz}
    bwa mem \
    -M \
    -t ~{cpu} \
    -R "~{actual_read_group}" \
    -k ~{min_seed_length} \
    ~{genome_index_basename}.fa \
    ~{fastq1_gz} ~{fastq2_gz} | \
    samtools view -Sb > ~{sample_name}.bam
  >>>

  output {
    File bam = "${sample_name}.bam"
  }

  runtime {
    docker: docker_image
    cpu: "${cpu}"
    memory: "${memory_gb} GB"
    disks: "local-disk ${actual_disk_gb} SSD"
    dx_timeout: "1D"
    dx_restart: object { max: 3 }
  }

  meta {
    title: "BWA-MEM"
    description: "Align paired-end reads using BWA MEM"
    details: { upstreamLicenses: "GPLv3" }
  }

  parameter_meta {
    sample_name: {
      label: "Sample Name",
      help: "Name of the sample; used to prefix output files"
    }
    fastq1_gz: {
      label: "FASTQ 1 (gzipped)",
      description: "Gzipped fastq file of first paired-end reads",
      stream: true
    }
    fastq2_gz: {
      label: "FASTQ 2 (gzipped)",
      description: "Gzipped fastq file of second paired-end reads",
      stream: true
    }
    genome_index_tgz: {
      label: "Genome Index (.tgz)",
      description: "Tarball of the reference genome and BWA index",
      stream: true
    }
    min_seed_length: {
      label: "Minimum Seed Length",
      help: "Matches shorter than INT will be missed.",
      group: "Advanced",
      default: 19
    }
    read_group: {
      label: "Read Group",
      help: "(Optional) the read group to add to aligned reads",
      group: "Advanced"
    }
    docker_image: {
      label: "Docker Image",
      help: "Name of the docker image to use",
      group: "Resources",
      default: "broadinstitute/baseimg"
    }
    cpu: {
      label: "CPUs",
      help: "Minimum number of CPUs to use",
      group: "Resources",
      default: 4
    }
    memory_gb: {
      label: "Memory (GB)",
      help: "Minimum amount of memory required",
      group: "Resources",
      default: 8
    }
    disk_gb: {
      label: "Disk Space (GB)",
      help: "Minimum amount of disk space required (in GB); by default this is calculated from the inputs",
      group: "Resources"
    }
  }
}