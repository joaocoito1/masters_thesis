
# Generates an annotated phylogenetic tree plot (colored by group, with 
# bootstrap support) for each IQ-TREE output. 
# Outputs one PNG + PDF per tree, and a list of groups without a tree.




##############
### Inputs ###
##############

inputs <- commandArgs(trailingOnly = TRUE)


# Path to folder with the original FASTA groups
groups_path <- inputs[1]

# Path to folder with tree files (*.treefile)
iqtree_path <- inputs[2]

# Path to folder where the output trees will be saved
trees_path <- inputs[3]

# Name of the file to be created with the groups without trees
groups_wo_trees <- inputs[4]

# Vector of group (clades) names
group_names <- trimws(strsplit(inputs[5], ",")[[1]])

# Vector of group (clades) colors
group_colors <- trimws(strsplit(inputs[6], ",")[[1]])

# Color of Recipient lineage
recipient_color <- inputs[7]

# Vector of HGT candidates
all_candidate_hgt_file <- inputs[8]
text <- paste(readLines(all_candidate_hgt_file), collapse = " ") # Remove all group IDs
all_candidate_hgt <- strsplit(text, "\\s+")[[1]] # Split on whitespace 
all_candidate_hgt <- all_candidate_hgt[nzchar(all_candidate_hgt)] # drop empty strings

# Path to log file
log_file <- inputs[9]


##############
###  Main  ###
##############

library(ggplot2)
library(ggtree)
library(ape)
library(phangorn)


##### Function to plot the trees #####

plot_tree_from_file<-function(tree_file, output_folder, group_names,
                              group_colors, recipient_label, recipient_color,
                              all_candidate_hgt, label_col_width = 0.01, 
                              base_size = 3, bootstrap_cutoff = 95){
  
  # Name colors with groups
  group_colors <- setNames(group_colors, group_names)
  recipient_color <- setNames(c(recipient_color), c(recipient_label))
  group_colors <- c(group_colors, recipient_color)
  
  # Read tree
  tree <- read.tree(tree_file)
  
  # Build group df
  tree_group <- data.frame(label=tree[["tip.label"]], group=tree[["tip.label"]])
  
  tree_group$group <- ifelse(
    tree_group$label %in% all_candidate_hgt,
    recipient_label,
    sapply(strsplit(tree_group$group, "\\|"), `[`, 2))
  
  
  # Root tree
  tree <- midpoint(tree)
  
  # Build ggtree object
  p <- ggtree(tree) %<+% tree_group
  
  # Change size based on label length and tip count
  max_label_len <- max(nchar(tree$tip.label))
  tree_depth <- max(p$data$x)
  x_max <- tree_depth * (1 + label_col_width * max_label_len)
  n_tips <- length(tree$tip.label)
  
  # Build the plot
  final_plot <- p +
    geom_tiplab(aes(fill = group),
                geom = "label",
                label.size = 0,
                color = "black",
                size = base_size) +
    geom_point2(aes(subset = as.numeric(label) > bootstrap_cutoff),
                size = 2, shape = 16, color = "black") +
    geom_treescale(x = 0, y = n_tips, fontsize = 3, linesize = 0.5) +
    scale_fill_manual(values = group_colors) +
    labs(fill = NULL) +
    guides(fill = guide_legend(override.aes = list(label = ""))) +
    theme(legend.position = "top",
          legend.text = element_text(margin = margin(l = 1, r = 10))) +
    coord_cartesian(clip = "off") +
    xlim(0, x_max)
  
  # Dynamic output dimensions
  out_width  <- max_label_len * 0.25
  out_height <- n_tips * 0.2

  dpi <- min(300, 30000 / max(out_width, out_height))
  dpi <- max(dpi, 30) # ensure that it does not go over the limit

  out_base <- paste0(sub("\\.treefile$", "", tree_file), "_tree")
  
  ggsave(paste0(out_base,".png"), path = output_folder, plot = final_plot,
         width = out_width, height = out_height, dpi = dpi, limitsize = FALSE)
  
  ggsave(paste0(out_base,".pdf"), path = output_folder, plot = final_plot,
         width = out_width, height = out_height, limitsize = FALSE)
}

##### Generate tree plots #####

setwd(iqtree_path)

tree_files <- grep("\\.treefile$", list.files(iqtree_path), value = TRUE)

n_trees <- length(tree_files) # total number of trees
n <- 0 # counter

# Run for every tree
for(file in tree_files){
  group <- sub("\\..*", "", file)
  plot_tree_from_file(file, trees_path, group_names, group_colors,
                      recipient_label="recipient",recipient_color,
                      all_candidate_hgt)

  n <- n + 1
  if (n %% 25 == 0) {
    message <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " : ", n, "/", n_trees, " trees created.")
    cat(message, "\n")
    cat(message, "\n", file = log_file, append = TRUE)
  }  

}

message <- paste0("\n",format(Sys.time(), "%Y-%m-%d %H:%M:%S")," : ","Tree plotting finished.")
cat(message, "\n")
cat(message, "\n", file = log_file, append = TRUE)


##### Check which groups do not have a tree #####

all_groups <- sub("\\.fasta$", "", list.files(groups_path, pattern = "\\.fasta$"))

groups_with_tree <- sub("\\..*", "", tree_files)

missing_groups <- setdiff(all_groups, groups_with_tree)

# Write one per line
writeLines(missing_groups, groups_wo_trees)

message <- paste0("\n",format(Sys.time(), "%Y-%m-%d %H:%M:%S")," : ",
                  length(missing_groups), " groups without a tree written to ", groups_wo_trees)
cat(message, "\n")
cat(message, "\n", file = log_file, append = TRUE)
