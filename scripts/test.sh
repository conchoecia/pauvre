set -ev

git clone https://github.com/wdecoster/nanotest.git

pauvre -h

pauvre marginplot -h
pauvre marginplot -f nanotest/reads.fastq.gz
pauvre marginplot -f nanotest/reads.fastq.gz -t "my title" --filt_maxlen 30000
pauvre marginplot -f nanotest/reads.fastq.gz -y --filt_minqual 7 --fileform svg

pauvre stats -h
pauvre stats -f nanotest/reads.fastq.gz
pauvre stats -f nanotest/reads.fastq.gz -H --filt_maxqual 12
