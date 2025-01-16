library(tidyverse)
library(ggtree)
library(treeio)
library(phangorn)
library(castor)
library(ggnewscale)

setwd("C:\\Users\\cassp\\Box Sync\\Feaga Lab\\Cassidy Prince\\snor")
file = "snor_HMMER_5e2.csv"

# Read in HMMER output as lines so the descriptions don't get cut off. 
df = data.frame(readLines(file)) %>%
  slice(4:n()) %>%
  separate(readLines.file., into = c(as.character(1:34)), sep = " +") %>%
  select(-2,-3,-4) %>%
  replace(is.na(.), " ") %>% #Replace NAs with spaces.
  unite("description", 16:31, sep = " ") %>%
  mutate_if(is.character, str_trim) %>% #Trim spaces off description.
  filter(`1` != "#") %>%
  mutate(species = str_extract(description, "(?<=\\[).*(?=\\])")) %>% #Extract species name from between [] in description. 
  `colnames<-`(c("accession", "full_eval", "full_score", "full_bias", "dom_eval", "dom_score", "dom_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description", "species")) %>% #Rename columns.
  mutate_at(vars(full_eval, full_score, full_bias, dom_eval, dom_score, dom_bias), as.numeric)

df_len = read_tsv("lengths.txt") %>%
  rename(accession = `#name`) %>%
  mutate(accession = str_replace(accession, "(?s) .*", ""))%>%
  left_join(df)

df_filt = df_len %>%
  filter(length < 200) %>%
  filter(full_score >= mean(full_score) - sd(full_score))

# Read in taxonomy and other metadata from NCBI Datasets.
df_md = left_join(read_tsv("fungi_tax_metadata.tsv"), read_tsv("fungi_metadata.tsv"), by = join_by(`Taxid` == `Organism Taxonomic ID`, `Tax name` == `Organism Name`)) %>%
  select(-Authority, -Rank, -`Basionym authority`, -`Curator common name`, -`Has type material`, -`Superkingdom name`, -`Superkingdom taxid`, -`Kingdom taxid`, -`Phylum taxid`, -`Class taxid`, -`Order taxid`, -`Family taxid`, -`Genus taxid`, -`Species taxid`, -`CheckM completeness`)

colnames(df_md)=str_replace_all(colnames(df_md), c(" " = "_" , "," = "" )) #Replace spaces in column names with underscores.

# Join HMMER data with metadata.
df = data.frame(cbind(species = df_md$Tax_name, snor = (df_md$Tax_name %in% df_filt$species))) %>%
  left_join(df_md, by = join_by(species == Tax_name))


table(df$Phylum_name, df$snor)
table(df$Group_name, df$snor)

ggplot(df_filt, aes(x = full_score)) +
  geom_histogram()


df_GCF_nc = read.csv("GCF_accessions_fungi_clean.csv") %>%
  select(-X)

df_blast = read_tsv("28S_blast_50.csv", skip = 5, col_names = c("query_acc", "subj_acc", "perc_id", "ali_len", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "eval", "bit_score")) %>%
  distinct(subj_acc, .keep_all = TRUE) %>%
  drop_na()
  
df_blast_GCF = left_join(df_blast, df_GCF_nc, by = join_by(`subj_acc` == `nuccore`)) %>%
  distinct(assembly, .keep_all = TRUE) 

write_tsv(df_blast_GCF, file = "fungi_blast_nd.tsv")




df_blast_md = left_join(df_blast_GCF, df, by = join_by(`assembly` == `Assembly Accession`))


table(df_blast_md$`Phylum name`, df_blast_md$snor)
table(df$`Phylum name`, df$snor)


dups = df_filt %>% 
  group_by(species) %>% 
  filter(n()>1)


#write.table(df_filt2$accession, file = "filt_accs.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)



###############treee

tree = read.newick("28S_fasttree.tre")



tree_tip_GCF = left_join(data.frame(tree$tip.label), df_GCF_nc, by = join_by("tree.tip.label" == "nuccore"), multiple = "any")

tree$tip.label = tree_tip_GCF$assembly

tree = get_subtree_with_tips(tree, only_tips = df$`Assembly Accession`)$subtree

tree_mid = midpoint(tree)

df1 = data.frame(df$snor)
rownames(df1) = df$`Assembly Accession`
df2 = data.frame(df$`Phylum name`)
rownames(df2) = df$`Assembly Accession`

p = ggtree(tree_mid, layout='circular', size=0.2) 

p1 = gheatmap(p, df1, offset=-0.05, width=0.1, font.size=1, colnames = FALSE, color=NA) + scale_fill_manual(values=c("TRUE" = "#961415", "FALSE" = "white", na.value = "white"), name = "SNOR presence") + 
  theme(text=element_text(size=18)) +
  new_scale_fill()

p2 = gheatmap(p1, df2, offset=0.2, width=0.1, font.size=1, colnames = FALSE, color=NA) +
  theme(text=element_text(size=15)) +
  scale_fill_discrete(name = "Phylum", na.value = "white")


p2



############### collapsed tree

df_small = left_join(tree_tip_GCF, df, by = join_by(assembly == `Assembly Accession`))

random = df_small %>%
  group_by(`Phylum name`) %>%
  sample_n(1) %>%
  rename(phylum = `Phylum name`) %>%
  select(-tree.tip.label)%>%
  as.data.frame()


write.csv(random, "random_genomes_for_tree.csv")
random = read.csv("random_genomes_for_tree.csv")

rownames(random) = random$assembly

col_tree = get_subtree_with_tips(tree, only_tips = random$assembly)$subtree

col_tree_mid = midpoint(col_tree)


table = as.data.frame(prop.table(table(df$`Phylum name`, df$snor), margin = 1)*100)

table_new = data.frame(cbind(snor = table$Freq[table$Var2 == "TRUE"]))
rownames(table_new) = random$assembly


table_for_tree(df$`Phylum name`, df$snor, random$assembly)





p = ggtree(col_tree_mid, size = 0.8) %<+% random + 
  xlim(0, 3) + 
  geom_tiplab(aes(label=phylum), align = TRUE, size = 4) + 
  geom_nodepoint(aes(fill = as.numeric(label)*100), size = 2, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage") + 
  new_scale_fill()

p


###Current problem is that I don't have the right number of phyla in my big tree (missing microsporidia and olpidiomycota [only one genome so this is fine] completely). maybe do a hand drawn cladogram? instead of a tree? idk
