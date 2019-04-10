workflow R2D2_Allelic_Imbalance_Workflow {
	File dnaGermline
	File dnaTumor
	File rnaGermline
	File rnaTumor
    Int diskSpace
    String pairId
    Int numThreads=8

    call AllelicImbalanceTask {
    	input:
           	dnaGermline=dnaGermline,
           	dnaTumor=dnaTumor,
           	rnaGermline=rnaGermline,
           	rnaTumor=rnaTumor,
           	pairId=pairId,
           	numThreads=numThreads
    }

    output {
    	AllelicImbalanceTask.categorizedVariants
    }
}

task AllelicImbalanceTask {
	File dnaGermline
	File dnaTumor
	File rnaGermline
	File rnaTumor
    String pairId
    Int numThreads

    command <<<
    	# log resource usage for debugging purposes
       	function runtimeInfo() {
        	echo [$(date)]
        	echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
        	echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
        	echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
        }
        while true;
        	do runtimeInfo;
           	sleep 15;
       	done &

        echo "pwd"
        pwd

        echo "ls -lh"
        ls -lh

        echo "python /r2d2v2.py ${dnaGermline} ${dnaTumor} ${rnaGermline} ${rnaTumor} ${pairId} --num_threads ${numThreads}"
        python /r2d2v2.py ${dnaGermline} ${dnaTumor} ${rnaGermline} ${rnaTumor} ${pairId} --num_threads ${numThreads}
    >>>

    output {
    	File categorizedVariants="${pairId}.categorized.maf"
    }

    runtime {
    	docker: "vanallenlab/r2d2_allelic_imbalance:1.0"
        memory: "3 GB"
        disks: "local-disk 2 HDD"
        bootDiskSizeGb: 20
        preemptible: 5
    }
}