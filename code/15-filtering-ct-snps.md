Filtering CT SNPs
================
Steven Roberts
07 August, 2023

- <a href="#1-locating-5x" id="toc-1-locating-5x">1 Locating 5x</a>

Need to take all 5x bismark tab outputs and filter to remove putative CT
SNPs..

# 1 Locating 5x

Grabbing from
<https://gannet.fish.washington.edu/seashell/bu-mox/scrubbed/032120-fds/>

``` bash
wget -r \
--no-directories --no-parent \
-P ../data \
-A "*_5x.tab" https://gannet.fish.washington.edu/seashell/bu-mox/scrubbed/032120-fds/
```

``` bash
head -20 ../data/*104*
```

``` bash
grep $'C\tT' ../output/SNP.vcf > ../output/CT-SNP.vcf
wc -l ../output/CT-SNP.vcf
```

``` bash
head ../output/CT-SNP.vcf
```

Need to be aware of potential difference in 0 base and 1 base. Will read
both in– in this case will just do a proof of concept on 104

``` r
test <- read.csv("../data/EPI-104_S28_L005_5x.tab", header = FALSE, sep = "\t")
```

``` r
test2 <- test %>%
  mutate(loci = paste0(V1, "_", V2))
```

``` r
ct <- read.csv("../output/CT-SNP.vcf", header = FALSE, sep = "\t")
```

``` r
ct2 <- ct %>%
  mutate(loci = paste0(V1, "_", V2))
```

``` r
inner_join(test2, ct2, by = "loci")
```

``` r
test2 %>%
  anti_join(ct2, by = "loci")
```

``` r
# 1. List all files with _5x.tab suffix
files <- list.files(path = "../data/", pattern = "_5x.tab$", full.names = TRUE)

# 2. Iterate over each file
for(file in files) {
  
  # Extract base filename without the directory for naming purposes
  base_name <- basename(file)
  
  # Read the file
  data <- read.csv(file, header = FALSE, sep = "\t")
  
  # Modify the data
  modified_data <- data %>%
    mutate(loci = paste0(V1, "_", V2)) %>%
    anti_join(ct2, by = "loci") %>%
    select(-loci)
    
  
  # Write the modified data to an output file
  output_file <- paste0("../output/f", base_name)
  write.table(modified_data, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
}
```
