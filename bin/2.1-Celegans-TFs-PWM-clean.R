pwm.lines <- read_lines(file = "experiments/2018-09-02-TFs-motifs/raw/CisBP_2018_09_02_12_51_pm/PWM.txt") # Read each line of the motifs file

# Extract TF names-------------------------------------------------------------
# Subset element that match pattern (grep), and remove the pattern from each row (gsub)
TF_Names <- gsub("TF Name\t","",pwm.lines[grep("TF Name\t",pwm.lines)])

# Remove non-PWM related elements----------------------------------------------

# Remove empty character elements
pwm.lines <- pwm.lines[!(pwm.lines == "")]

# Set row patterns to remove
toMatch <- c("TF\t","TF Name\t","Gene\t","Motif\t","Family\t","Species\t")
# Convert to grep expression format
toMatch <- paste(toMatch,collapse="|")
# Remove all pattern matches
pwm.lines <- pwm.lines[grep(toMatch,pwm.lines,invert = TRUE)]

# Separate PWM tables----------------------------------------------------------

# Obtain position of each PWM table's columns.
# Index different chunks of rows corresponding to the same table
group <- cumsum(grepl("Pos\tA\tC\tG\tT", pwm.lines))

# Split each PWM tables according to the above indexing (split) -- list object
# Read of the list items as a tab-separated table (read.table + lapply)
PWM <- lapply(split(pwm.lines, group), function(x) read.table(text=x, sep="\t",header = TRUE))

# Data cleaning----------------------------------------------------------------

# Set TF names to each PWM
names(PWM) <- TF_Names

# Set Pos column as the row names of each PWM, and remove its values.
PWM <- lapply(PWM, function(x) { row.names(x) <- as.character(x[,'Pos']); x[,2:5]} )

# Remove *tra-1* as it is an empty PWM
PWM <- PWM[TF_Names != "tra-1"]

# Transpose the matrices, as that is the acceptable input format 
# (rows = bases, columns = positions).
PWM <- lapply(PWM,t)

# Save PWM as .rdata-----------------------------------------------------------
dir.create("experiments/2018-09-02-TFs-motifs/data")
save(PWM,file = "experiments/2018-09-02-TFs-motifs/data/TFs-PWM.rdata")