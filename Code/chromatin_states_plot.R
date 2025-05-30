# chromatin_states_plot.R
# Visualize chromatin state distribution with pie and bar plots

# Install and load required package
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
library(ggplot2)

# Define input file
input_file <- "chromatin_states_count"

# Check file exists
if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

# Load data
data <- read.table(input_file, header = FALSE, col.names = c("Counts", "Chromatin_States"))

# Define chromatin state categories
chromatin_states <- list(
  "TSS" = c("Tss"),
  "Flanking TSS" = c("TssFlnk", "TssFlnkD", "TssFlnkU"),
  "Transcription" = c("Tx", "TxWk"),
  "Enhancers" = c("Enh1", "Enh2", "EnhG1", "EnhG2"),
  "ZNF genes & repeats" = c("ZNF/Rpts"),
  "Heterochromatin" = c("Het"),
  "Bivalent/Poised TSS" = c("Biv"),
  "Repressed PolyComb" = c("ReprPC"),
  "Quiescent/Low" = c("Quies")
)

# Map each chromatin state to its functional category
data$Category <- sapply(data$Chromatin_States, function(state) {
  match <- names(chromatin_states)[sapply(chromatin_states, function(x) state %in% x)]
  if (length(match) > 0) return(match) else return("Unknown")
})

# Aggregate counts by category
category_data <- aggregate(Counts ~ Category, data = data, sum)

# Calculate percentages and labels
category_data$Percentages <- round((category_data$Counts / sum(category_data$Counts)) * 100, 1)
category_data$Labels <- paste(category_data$Category, category_data$Percentages, "%")

# Generate a pastel color palette
pastel_colors <- colorRampPalette(c("#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF", "#D7BAFF"))(nrow(category_data))

# --- PIE CHART ---
pie(
  category_data$Counts,
  labels = category_data$Labels,
  col = pastel_colors,
  main = "Distribution of Chromatin State Categories",
  cex.main = 1.5
)

# --- BARPLOT ---
par(mar = c(8, 5, 5, 2))  # Margins: bottom, left, top, right

bar_heights <- barplot(
  category_data$Percentages,
  names.arg = category_data$Category,
  col = pastel_colors,
  las = 2,  # Rotate x-axis labels
  main = "Distribution of Chromatin State Categories",
  ylab = "Percentage (%)",
  ylim = c(0, max(category_data$Percentages) + 10),
  cex.names = 0.8
)

# Add percentage labels above bars
text(
  x = bar_heights,
  y = category_data$Percentages + 1,
  labels = paste0(category_data$Percentages, "%"),
  pos = 3,
  cex = 0.9,
  col = "black"
)

# --- LEGEND ---
par(mar = c(0, 0, 0, 0))  # No margins
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")  # Blank canvas

# Prepare legend text and colors
legend_text <- sapply(sort(names(chromatin_states)), function(cat) {
  paste(cat, ":", paste(chromatin_states[[cat]], collapse = ", "))
})
legend_colors <- pastel_colors[1:length(legend_text)]

legend(
  "center",
  legend = legend_text,
  fill = legend_colors,
  bty = "n",
  title = "Chromatin States",
  cex = 0.9
)

