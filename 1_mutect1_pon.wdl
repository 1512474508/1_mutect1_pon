workflow M1_PON {
	# inputs	
	Array[File] normal_bams
	Array[File] normal_bais
	File ref_fasta
	File ref_fai
	File ref_dict
	String mutect1_docker   
	String gatk_docker  
	
	Array[Pair[File,File]] normal_bam_pairs = zip(normal_bams, normal_bais)

    scatter (normal_bam_pair in normal_bam_pairs) {
        File normal_bam = normal_bam_pair.left
        File normal_bai = normal_bam_pair.right
        
		call M1 {
			input:
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                mutect1_docker = mutect1_docker
		}		
		
		call SelectVariants {
		input:
            input_vcf = M1.output_pon_vcf,
            input_vcf_idx = M1.output_pon_vcf_index,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            gatk_docker = gatk_docker,
            VCF_pre = M1.fname
	    }
	}
    
    call CreateFoFN {
        input:
            array_of_files = SelectVariants.filtered_vcf,
    }
    
	call CombineVariants {
		input:
            filtered_vcfs_list = CreateFoFN.fofn_list,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            gatk_docker = gatk_docker
	}

    output {
        Array[File] output_pon_vcf = M1.output_pon_vcf
        Array[File] output_pon_vcf_index = M1.output_pon_vcf_index        
        Array[File] output_pon_stats_txt = M1.output_pon_stats_txt
        Array[File] filtered_vcf = SelectVariants.filtered_vcf
        File filtered_vcfs_list = CreateFoFN.fofn_list
        File pon = CombineVariants.pon
    }

	meta {
		author: "Sehyun Oh"
        email: "shbrief@gmail.com"
        description: "Build Pool of Normal (PoN) using MuTect v.1.1.7 for PureCN"
    }
}

task M1 {
	# input
	File normal_bam
	File normal_bai
	File dbsnp_vcf
	File dbsnp_vcf_idx
	File cosmicVCF
	File cosmicVCF_idx
	File ref_fasta
	File ref_fai
	File ref_dict

	String BAM_pre = basename(normal_bam, ".bam")

	# runtime
	String mutect1_docker
    Int disk_size = ceil(size(normal_bam, "GB")) + 30

	command <<<
		java -jar -Xmx4g /home/mutect-1.1.7.jar \
		--analysis_type MuTect \
		-R ${ref_fasta} \
		--artifact_detection_mode \
		--dbsnp ${dbsnp_vcf} \
		--cosmic ${cosmicVCF} \
		-dt None \
		-I:tumor ${normal_bam} \
		-o ${BAM_pre}_pon_stats.txt \
		-vcf ${BAM_pre}_pon.vcf
	>>>

	runtime {
		docker: mutect1_docker   # jvivian/mutect
		memory: "32 GB"
		disks: "local-disk " + disk_size + " HDD"
	}

    output {
        String fname = "${BAM_pre}"
        File output_pon_vcf = "${BAM_pre}_pon.vcf"
        File output_pon_vcf_index = "${BAM_pre}_pon.vcf.idx"        
        File output_pon_stats_txt = "${BAM_pre}_pon_stats.txt"
    }
}

task SelectVariants {
    # input
    File input_vcf
    File input_vcf_idx
    File ref_fasta
    File ref_fai
    File ref_dict
    String VCF_pre

    # runtime
    String gatk_docker
    
    command <<<
        java -jar -Xmx4g /usr/GenomeAnalysisTK.jar \
        --analysis_type SelectVariants \
        -R ${ref_fasta} \
        -V ${input_vcf} \
        -o ${VCF_pre}_pon.vcf
    >>>
        
        runtime {
            docker: gatk_docker
            memory: "8 GB"
        }
    
    output {
        File filtered_vcf = "${VCF_pre}_pon.vcf"
    }
}

task CreateFoFN {
    # Command parameters
    Array[String] array_of_files
    String fofn_name
    
    command <<<
        mv ${write_lines(array_of_files)}  ${fofn_name}.list
    >>>
    
    output {
        File fofn_list = "${fofn_name}.list"
    }
    
    runtime {
        docker: "ubuntu:latest"
    }
}

task CombineVariants {
    # input
    File filtered_vcfs_list
    Array[File] input_vcfs = read_lines(filtered_vcfs_list)
    File ref_fasta
    File ref_fai
    File ref_dict
    
    # runtime
    String gatk_docker
    
    command <<<
        java -jar -Xmx24g /usr/GenomeAnalysisTK.jar \
        -T CombineVariants \
        -nt 4 --minimumN 5 --genotypemergeoption UNSORTED \
        -R ${ref_fasta} \
        -V ${sep=' -V ' input_vcfs} \
        -o "normals.merged.min5.vcf"
    >>>
        
        runtime {
            docker: gatk_docker
            memory: "32 GB"
        }
    
    output {
        File pon = "normals.merged.min5.vcf"
    }
}