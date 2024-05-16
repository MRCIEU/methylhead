### Preparation of Beta Cell Types for Estimating cell counts

*The script reads blood panel data from https://github.com/MRCIEU/dnam-lung-cancer-screening-panel/blob/master/panel-reduced.csv
*It writes the unique blood panel data to a CSV file
*The script iterates through each row of the blood panel to extract genomic regions.
*It uses the 'wgbstools' command-line tool to extract data for each region and appends it to a file.
*After that, it processes the extracted data to create a BED file containing relevant columns.
*Finally, it converts beta values to table format using 'wgbstools' and saves the result to a CSV file.
'wgbstools' tool was used for data extraction and conversion. https://github.com/nloyfer/wgbs_tools

# This script prepares beta cell types files for analysis.

# Read blood panel data
blood_panel <- fread("panel-reduced_BS.csv")

# Remove rows where start equals end
blood_panel <- subset(blood_panel, blood_panel$start != blood_panel$end)

# Write unique blood panel data to a CSV file
write.csv(blood_panel, "blood_panel_unique.csv", row.names = FALSE)

# Assign blood panel to another variable
panel <- blood_panel

# Iterate through each row of the panel to extract regions
for (i in 1:nrow(panel)) {
  # Define region
  region <- paste0("chr", panel$chr[i], ":", panel$start[i], "-", panel$end[i])
  
  # Execute system command to extract data and append to file
  system(paste("./wgbstools convert -r", region, "| head -n1 >> panel-reduced-BS-sites.txt"))
}

# Extract relevant columns from the file and create a BED file
awk 'BEGIN {FS="[:, -]"; OFS="\t"} {start=$1+$2; end=$3+$4; print $1, start,end , $11,$12}' panel-reduced-BS-sites.txt | tail -n +1 > blood_cell.bed

# Convert beta values to table format and save to CSV file
./wgbstools beta_to_table blood_cell.bed --betas *.beta | column -t > blood_cell_types.csv
