awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' coronary_consensus_peaks.bed > coronary_consensus_peaks.saf
