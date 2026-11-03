cat ../tss_enrichment/tss_time.txt | grep "wall" | awk -F ' ' '{print "tss\ttotal\t"$NF}' > tss_cov_time.txt
cat ../tss_enrichment/tss_io_time.txt | grep "wall" | awk -F ' ' '{print "tss\tio\t"$NF}' >> tss_cov_time.txt
cat ../pseudodobulk_fragment_coverage/coverage_time.txt | grep "wall" | awk -F ' ' '{print "cov\ttotal\t"$NF}' >> tss_cov_time.txt
cat ../pseudodobulk_fragment_coverage/coverage_io_time.txt | grep "wall" | awk -F ' ' '{print "cov\tio\t"$NF}' >> tss_cov_time.txt

cat ../tss_enrichment/tss_time.txt | grep "Maximum resident set size" | awk -F ' ' '{print "tss\ttotal\t"$NF}' > tss_cov_mem.txt
cat ../tss_enrichment/tss_io_time.txt | grep "Maximum resident set size" | awk -F ' ' '{print "tss\tio\t"$NF}' >> tss_cov_mem.txt
cat ../pseudodobulk_fragment_coverage/coverage_time.txt | grep "Maximum resident set size" | awk -F ' ' '{print "cov\ttotal\t"$NF}' >> tss_cov_mem.txt
cat  ../pseudodobulk_fragment_coverage/coverage_io_time.txt | grep "Maximum resident set size" | awk -F ' ' '{print "cov\tio\t"$NF}' >> tss_cov_mem.txt
