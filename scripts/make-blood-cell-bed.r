# Read blood panel data
library(data.table) 

blood_panel <- data.frame(fread("https://github.com/MRCIEU/dnam-lung-cancer-pipeline/blob/main/data/blood_cell_types_extended.bed")) 
colnames(blood_panel) <- c("chr", "start", "end")

# Iterate through each row of the panel to extract regions
for (i in 1:nrow(blood_panel )) {
  # Define region
  region <- paste0(blood_panel$chr[i], ":", blood_panel$start[i], "-", blood_panel$end[i])
  # Execute system command to extract data and append to file
  system(paste("./wgbstools convert -r", region, "| head -n1 >> panel-BS-sites.txt"))
}

# Extract relevant columns from the file and create a BED file
awk 'BEGIN {FS="[:-]"; OFS="\t"} {start=$1+$2; end=$3+$4; print $1, start, end, $5, $6}' panel-BS-sites.txt | tail -n +1 > blood-cell-sites.bed 

# Convert beta values to table format and save to CSV file
./wgbstools beta_to_table  blood-cell-sites.bed --betas *.beta | column -t > blood_cell_types.csv
