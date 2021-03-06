# Analysis scripts for Baez-Ortega et al., 2019
# Step 5. Definition of phylogenetic tumour groups and group-specific variants

# Adrian Baez-Ortega, 2018


# TO RUN THIS SCRIPT IN THE TERMINAL
# ----------------------------------
# Run the commands below in the terminal, replacing '/path/to/TCG2019' with 
# the path to the TCG2019 directory.
#
#    cd /path/to/TCG2019
#    Rscript scripts/5_PhyloGroups.R


# TO RUN THIS SCRIPT IN RSTUDIO
# -----------------------------
# Before starting, run the line below in RStudio, replacing '/path/to/TCG2019' with 
# the path to the TCG2019 directory.
#
#    setwd("path/to/TCG2019")


# If the paths to the input or output files differ from the ones 
# used in this script, they may be updated by modifying the lines below.

INPUT = list(
    
    # Path to input file containing the ML phylogenetic tree inferred by RAxML
    # (see Supplementary Methods in Baez-Ortega et al., 2018)
    TREE = file.path("data", "original", "RAxML_CTVT_MLTree.newick.tree"),
    
    # Path to input file containing country and town information per sample
    COUNTRIES = file.path("data", "original", "Sample_Country_Town.tsv"),
    
    # Path to input file generated by script 2_ImportVariants.R
    VAR.TABLES = file.path("data", "processed", "Variant_Tables.RData")
    
)

OUTPUT = list(
    
    # Path to output tree plot file
    TREE.PDF = file.path("output", "CTVT_Tree_RAxML.pdf"),
    
    # Path to output mutational spectra file
    SPECTRA.PDF = file.path("output", "Phylo_Groups_Spectra.pdf"),
    
    # Path to output RData file
    RDATA = file.path("data", "processed", "Phylo_Tree_Groups.RData")
    
)


# Print input and output file paths
cat("\nInput files:")
for (name in INPUT) {
    cat("\n  ", name)
}
cat("\n\nOutput files:")
for (name in OUTPUT) {
    cat("\n  ", name)
}
cat("\n\n")


# Load packages
PACKAGES = c("ape", "sigfit", "stringr")
cat("Loading packages:", paste(PACKAGES, collapse=", "), "\n")
for (package in PACKAGES) {
    suppressWarnings(library(package, character.only=TRUE, quietly=TRUE))
}


# Constant: number of 'normal' phylogenetic groups
NG = 58


# Load input data
cat("Loading data...\n")
load(INPUT$VAR.TABLES)
tree = read.tree(INPUT$TREE)
country.town = read.table(INPUT$COUNTRIES, sep="\t", header=T, stringsAsFactors=F)
stopifnot(all.equal(country.town$Sample, samples[tumours]))


# Make country labels for every tumour (including towns only for India)
countries = ifelse(country.town$Country == "India",
                   paste0(country.town$Country, " (", country.town$Town, ")"),
                   country.town$Country)


# Replace tree node labels
node.idx = match(tree$tip.label, samples[tumours])
tree$tip.label = paste(samples[tumours], countries, sep=" - ")[node.idx]
tree$tip.label[length(tree$tip.label)] = "REF (CanFam3.1)"


# Plot tree with country labels
cat("\nPlotting phylogenetic tree to output directory...\n")
pdf(OUTPUT$TREE.PDF, 20, 100)
plot(tree, font=1, label.offset=0.0000005)
title(paste0("CTVT phylogram (RAxML)\n", length(tree$tip.label)-1, " tumours"), 
      cex.main=2.5, line=-8)
invisible(dev.off())


# Define phylogenetic tumour groups
# Groups were manually defined based in the tumour tree, based on the distribution of
# somatic variants; every group should have at least ~700 unique SNVs.
# (See Supplementary Methods in Baez-Ortega et al., 2018.)
# We define two types of phylogenetic groups:
#  - Normal groups consist of tree branches and their corresponding tips (tumours),
#    which should have at least ~700 unique variants and, if possible, represent a 
#    single geographical location. Identified by a number (e.g. [15]).
#  - Ancestral groups are composed of ancestral variant sets, i.e. the roots of
#    certain branches in the tree. Identified by 'A' and a number (e.g. [A2]).
cat("\nDefining phylogenetic tumour groups:\n")

# Group identifiers (normal and ancestral)
group.ids = c(1:NG, paste0("A", 1:3))

# Group names (normal, extra, and excluded)
group.names = paste0(c("India (Kolkata)", "Nicaragua", "Mexico", "Belize-Mexico", "Mexico-USA (Tucson)", "Nicaragua", "Nicaragua", "Nicaragua", "Belize", 
                       "Nicaragua-Guatemala-Honduras", "Mexico-USA (Tucson)", "Colombia-Panama", "Ecuador", "Nicaragua-Colombia (San Andres)", "Ecuador", 
                       "Chile", "Paraguay-Uruguay", "The Gambia", "Suriname", "Grenada", "Grenada", "India (Jaipur)", "Pakistan", "India (Leh)", 
                       "India (Jaipur)", "Thailand", "India (Tamil Nadu)", "India (Bodhgaya-Bylakuppe-Tamil Nadu)", "Brazil", "Brazil-Uruguay", "Brazil", 
                       "Brazil", "South Africa", "South Africa", "South Africa", "Lesotho-South Africa", "Italy", "Venezuela", "Venezuela", "Grenada", 
                       "Reunion-Senegal", "Paraguay", "Greece-Turkey-Tanzania-Kenya-Uganda", "Mauritius", "Suriname-Cape Verde-Portugal", 
                       "Tanzania-Turkey-Romania-Russia", "USA (St Louis)", "Australia", "Australia", "Colombia", "Malawi", "Ukraine-Romania", 
                       "Russia (Yoshkar-Ola)", "Russia-Ukraine", "Ukraine", "Armenia-Turkey", "India (Jaipur)", "India (Jaipur)",
                       "Basal trunk (ancestral)", "Non-India (ancestral)", "India (ancestral)"),
                     " [", 
                     group.ids, 
                     "]")

# Select tumours present in the tree
tumours.tree = samples %in% str_split_fixed(tree$tip.label, " - ", 2)[,1]

# Initialise group index matrix
phylo.groups.idx = matrix(FALSE, nrow=sum(tumours.tree), ncol=length(group.names),
                          dimnames=list(samples[tumours.tree], group.names))

# Assign tumours to groups through group indices (normal groups)
phylo.groups.idx[c("1322T", "13T", "12T"), 
                 1] = TRUE
phylo.groups.idx[c("1274T", "1277T"), 
                 2] = TRUE
phylo.groups.idx[c("49T", "53T", "56T", "1485T", "1263T"), 
                 3] = TRUE
phylo.groups.idx[c("1340T", "740T", "839T1", "464T", "782T", "739T", "738T", "59T", "170T"), 
                 4] = TRUE
phylo.groups.idx[c("176T", "175T", "1264T", "1262T", "57T", "54T", "55T", "60T", "171T", "180T", "1387T1"), 
                 5] = TRUE
phylo.groups.idx[c("550T", "572T", "1298T", "1297T", "564T", "555T", "1304T", "551T", "540T", "1283T", 
                   "1284T", "1290T", "1311T", "1288T", "1280T", "571T", "1300T", "1272T", "1312T", 
                   "1273T", "93T", "1281T", "534T", "92T", "350T1", "1216T", "558T"), 
                 6] = TRUE
phylo.groups.idx[c("559T", "1315T", "1275T", "575T", "1314T", "556T"), 
                 7] = TRUE
phylo.groups.idx[c("100T", "98T", "99T", "96T", "95T"), 
                 8] = TRUE
phylo.groups.idx[c("815T", "820T", "823T", "452T", "817T", "511T", "453T", "459T", "824T", "819T", "821T", "816T", "814T"), 
                 9] = TRUE
phylo.groups.idx[c("525T", "527T", "523T", "644T", "650T", "649T", "645T", "647T", "432T", "535T", "530T1", "94T", "349T2"), 
                 10] = TRUE
phylo.groups.idx[c("1389T", "1391T", "928T1"), 
                 11] = TRUE
phylo.groups.idx[c("1353T", "1530T", "1539T", "1535T", "1538T", "1250T", "1534T", "1533T", "1251T", "1536T", "463T", "414T"), 
                 12] = TRUE
phylo.groups.idx[c("425T", "421T", "326T", "1104T", "377T", "379T", "621T", "422T"), 
                 13] = TRUE
phylo.groups.idx[c("1217T", "1213T", "1210T", "1218T", "1215T", "319T", "532T", "531T"), 
                 14] = TRUE
phylo.groups.idx[c("625T", "624T", "623T", "622T", "423T", "1106T"), 
                 15] = TRUE
phylo.groups.idx[c("1160T", "582T", "961T", "962T", "372T", "420T", "1159T1", "959T", "1161T"), 
                 16] = TRUE
phylo.groups.idx[c("360T", "364T", "361T", "1245T", "1247T"), 
                 17] = TRUE
phylo.groups.idx[c("635T", "1173T", "666T", "813T", "469T", "1172T", "812T", "1024T1", "609T", "608T", "1023T", "667T", "468T", "639T1"), 
                 18] = TRUE
phylo.groups.idx[c("388T", "501T", "383T", "499T", "393T", "507T", "495T", "504T", "391T", "387T", "509T", "390T", "502T1"), 
                 19] = TRUE
phylo.groups.idx[c("124T", "126T"), 
                 20] = TRUE
phylo.groups.idx[c("862T", "859T", "861T", "342T", "196T", "131T"), 
                 21] = TRUE
phylo.groups.idx[c("980T", "411T", "313T", "485T", "986T", "400T", "410T", "1397T", "396T", "1398T", "317T", 
                   "981T", "484T", "1393T", "1399T", "310T", "409T", "308T", "974T", "315T", "1395T1", 
                   "475T", "1396T", "479T", "309T", "481T", "404T", "482T", "1318T", "407T", "985T"), 
                 22] = TRUE
phylo.groups.idx[c("714T", "711T1", "704T1", "710T", "706T1"), 
                 23] = TRUE
phylo.groups.idx[c("154T", "146T", "103T", "106T", "279T", "101T", "331T", "234T", "149T", "102T", "108T", 
                   "109T", "104T", "275T", "152T", "153T", "226T", "151T"), 
                 24] = TRUE
phylo.groups.idx[c("402T", "975T", "406T"), 
                 25] = TRUE
phylo.groups.idx[c("1259T", "1261T", "1224T1", "1224T2", "1260T", "1225T"), 
                 26] = TRUE
phylo.groups.idx[c("831T", "832T", "441T", "438T", "827T", "443T", "829T", "439T", "828T", "435T"), 
                 27] = TRUE
phylo.groups.idx[c("830T", "1178T", "1522T", "1523T", "281T"), 
                 28] = TRUE
phylo.groups.idx[c("757T", "756T", "769T1", "86T", "772T1", "758T", "82T", "84T", "80T", "775T2", "771T1"), 
                 29] = TRUE
phylo.groups.idx[c("773T1", "1246T"), 
                 30] = TRUE
phylo.groups.idx[c("1051T", "1049T", "1050T", "764T", "1077T", "763T", "83T", "774T1"), 
                 31] = TRUE
phylo.groups.idx[c("770T1", "766T1"), 
                 32] = TRUE
phylo.groups.idx[c("1080T", "791T", "691T", "1240T", "680T", "1517T", "1443T", "684T", "1206T", "867T", "1207T", "1424T", 
                   "1235T", "873T", "1083T", "1088T", "1082T", "1512T", "681T", "1448T", "934T", "1430T1", "1354T", "1444T"), 
                 33] = TRUE
phylo.groups.idx[c("1515T", "939T", "1358T", "1230T", "1504T", "1362T", "1422T", "1086T", "1421T", "1500T", "1508T"), 
                 34] = TRUE
phylo.groups.idx[c("1359T", "1089T", "693T", "937T", "1439T1", "1090T1", "1190T", "793T", "792T", "874T", "687T", "866T", "694T", 
                   "1084T", "1509T", "1447T", "1365T", "1510T", "1496T", "1513T", "875T", "795T", "1445T", "790T", "1438T1", 
                   "1195T", "933T", "1433T1", "1433T3", "798T", "1503T", "1099T", "1094T"), 
                 35] = TRUE
phylo.groups.idx[c("642T1", "1423T"), 
                 36] = TRUE
phylo.groups.idx[c("18T1", "17T1", "15T1"), 
                 37] = TRUE
phylo.groups.idx[c("595T1", "584T", "598T", "596T", "590T", "591T", "594T", "597T"), 
                 38] = TRUE
phylo.groups.idx[c("586T", "585T"), 
                 39] = TRUE
phylo.groups.idx[c("850T", "339T", "853T", "130T", "134T", "847T", "855T", "128T", "344T", "335T", 
                   "341T", "858T", "856T", "343T", "851T", "338T", "345T"), 
                 40] = TRUE
phylo.groups.idx[c("462T", "996T2"), 
                 41] = TRUE
phylo.groups.idx[c("366T", "363T", "367T", "365T"), 
                 42] = TRUE
phylo.groups.idx[c("415T", "604T", "605T", "1326T1", "416T", "925T", "488T2", "631T", "431T1", "122T", "209T", "120T1", "634T"), 
                 43] = TRUE
phylo.groups.idx[c("354T", "357T", "356T1", "355T", "351T", "352T", "353T"), 
                 44] = TRUE
phylo.groups.idx[c("33T", "35Ta", "386T", "385T", "506T1", "500T", "382T", "505T", "32T"), 
                 45] = TRUE
phylo.groups.idx[c("324T2", "321T2", "1333T1", "801T", "1335T", "1167T1"), 
                 46] = TRUE
phylo.groups.idx[c("1208T"), 
                 47] = TRUE
phylo.groups.idx[c("512T", "514T", "89T", "44T", "23T", "43T", "1041T", "1042T", "1040T", "136T", "137T1", "140T", 
                   "143T", "1002T", "265T", "286T1", "515T", "1026T", "513T", "1000T", "135T", "39T", "25T", "40T", 
                   "22T", "37T", "1001T", "90T", "45T1", "267T", "139T", "142T", "48T1", "47T", "41T", "29T1"), 
                 48] = TRUE
phylo.groups.idx[c("1468T", "1471T", "1469T", "28T", "1470T", "516T", "285T", "26T"), 
                 49] = TRUE
phylo.groups.idx[c("1351T", "1529T", "1537T", "1531T", "1532T"), 
                 50] = TRUE
phylo.groups.idx[c("614T", "613T", "656T", "615T", "611T"), 
                 51] = TRUE
phylo.groups.idx[c("833T1", "800T", "1166T", "799T", "460T", "1334T", "493T1", "492T1", "616T", "629T"), 
                 52] = TRUE
phylo.groups.idx[c("834T1"), 
                 53] = TRUE
phylo.groups.idx[c("1337T", "1336T", "835T", "602T", "1518T", "660T", "809T", "941T"), 
                 54] = TRUE
phylo.groups.idx[c("491T1", "446T", "450T", "1482T", "490T"), 
                 55] = TRUE
phylo.groups.idx[c("1521T", "1525T", "1332T"), 
                 56] = TRUE
phylo.groups.idx[c("399T", "982T"), 
                 57] = TRUE
phylo.groups.idx[c("401T", "976T", "977T"), 
                 58] = TRUE

# For the ancestral groups, select all samples descending from each, 
# including those that are not included in any of the normals group
A3.tumours = rowSums(phylo.groups.idx[, 57:58]) | rownames(phylo.groups.idx) == "270T"
phylo.groups.idx[, match("A3", group.ids)] = A3.tumours
phylo.groups.idx[, match("A2", group.ids)] = !A3.tumours
phylo.groups.idx[, match("A1", group.ids)] = rep(TRUE, sum(tumours.tree))

print(structure(cbind(group.ids, group.names, colSums(phylo.groups.idx)),
                dimnames=list(NULL, c("ID", "Name", "Samples"))), quote=FALSE)


# Identify group-unique variants
# We define group-unique variants of a group as those tumour-only variants that
# show support in ≥1 sample within the group, and no support outside the group
cat("\nIdentifying group-unique variants...\n")

group.unique.idx = matrix(FALSE, nrow=sum(tumour.only.idx), ncol=length(group.names),
                          dimnames=list(NULL, colnames(phylo.groups.idx)))
# Normal groups
for (j in 1:NG) {
    grp.idx = phylo.groups.idx[, j]
    if (sum(grp.idx) > 1) {
        group.unique.idx[, j] = 
            rowSums(snvs.nv[tumour.only.idx, tumours.tree][, grp.idx] >= MIN.READS) > 0 &
            rowSums(snvs.nv[tumour.only.idx, tumours.tree][, !grp.idx] >= MIN.READS) == 0
    }
    else if (sum(grp.idx) == 1) {
        group.unique.idx[,j] = 
            snvs.nv[tumour.only.idx, tumours.tree][, grp.idx] >= MIN.READS &
            rowSums(snvs.nv[tumour.only.idx, tumours.tree][, !grp.idx] >= MIN.READS) == 0
    }
    cat("  ", group.names[j], ": ", sum(group.unique.idx[,j]), " unique variants\n", sep="")
}

# Ancestral groups: we define a min. number of samples required for ancestral variants
# Basal trunk (ancestral) [A1]
# Min. samples for ancestral vars.: 534/539
MIN.SAMPLES = 534
j = match("A1", group.ids)
grp.idx = phylo.groups.idx[, j]
group.unique.idx[, j] = 
    rowSums(snvs.nv[tumour.only.idx, tumours.tree][, grp.idx] >= MIN.READS) >= MIN.SAMPLES
cat("  ", group.names[j], ": ", sum(group.unique.idx[, j]), " unique variants\n", sep="")

# India (ancestral) [A3]
# Min. samples for ancestral vars.: 4/6
MIN.SAMPLES = 4
j = match("A3", group.ids)
grp.idx = phylo.groups.idx[, j]
group.unique.idx[, j] = 
    rowSums(snvs.nv[tumour.only.idx, tumours.tree][, grp.idx] >= MIN.READS) >= MIN.SAMPLES &
    rowSums(snvs.nv[tumour.only.idx, tumours.tree][, !grp.idx] >= MIN.READS) == 0
cat("  ", group.names[j], ": ", sum(group.unique.idx[, j]), " unique variants\n", sep="")

# Non-India (ancestral) [A2] (special case)
# Min. samples for ancestral vars.: ≥1 sample in [1] + ≥1 sample in [54-56] + 0 samples in [A3]
j1 = 1
j2 = 54:56
j3 = match("A3", group.ids)
k = match("A2", group.ids)
grp.idx1 = phylo.groups.idx[, j1]
grp.idx2 = as.logical(rowSums(phylo.groups.idx[, j2]))
grp.idx3 = phylo.groups.idx[, j3]
ancs.idx = 
    rowSums(snvs.nv[tumour.only.idx, tumours.tree][, grp.idx1] >= MIN.READS) > 0 &
    rowSums(snvs.nv[tumour.only.idx, tumours.tree][, grp.idx2] >= MIN.READS) > 0 &
    rowSums(snvs.nv[tumour.only.idx, tumours.tree][, grp.idx3] >= MIN.READS) == 0
group.unique.idx[, k] = ancs.idx
cat("  ", group.names[k], ": ", sum(group.unique.idx[, k]), " unique variants\n", sep="")


# Identify and discard variants that are ancestral to groups with long trunks
# (this is to ensure that variants occurred after the group's arrival to its
# country; see Supplementary Methods)
cat("\nDiscarding ancestral variants from some groups...\n")

# List of groups to correct and min. #samples required for ancestral variants:
# The Gambia [18] - Min. #samples for ancestral vars.: 11/14
# Suriname [19]   - Min. #samples for ancestral vars.: 11/13
# Italy [37]      - Min. #samples for ancestral vars.: 3/3
# Grenada [40]    - Min. #samples for ancestral vars.: 10/17
# Mauritius [44]  - Min. #samples for ancestral vars.: 7/7
# Malawi [51]     - Min. #samples for ancestral vars.: 5/5
# Ukraine-Romania [52] - Min. #samples for ancestral vars.: 7/10
group.list = list(c(18, 11), c(19, 11), c(37, 3), c(40, 10),
                  c(44, 7), c(51, 5), c(52, 7))

for (grp in group.list) {
    j = grp[1]
    min.samples = grp[2]
    grp.idx = phylo.groups.idx[, j]
    ancs.idx = 
        rowSums(snvs.nv[tumour.only.idx, tumours.tree][, grp.idx] >= MIN.READS) >= min.samples &
        rowSums(snvs.nv[tumour.only.idx, tumours.tree][, !grp.idx] >= MIN.READS) == 0
    group.unique.idx[, j] = group.unique.idx[, j] & !ancs.idx
    cat("  ", group.names[j], ": ", sum(ancs.idx), " ancestral variants discarded\n", sep="")
}


# Use sigfit to plot spectra of group-unique variant sets
cat("\nPlotting mutational spectra of group-unique variants to output directory...\n")

unique.vars.idx = unlist(apply(group.unique.idx, 2, which))
trinucs = substr(str_split_fixed(snvs.metadata$INFO, ";", 20)[, 12], 13, 15)
unique.vars.table = data.frame(group = rep(group.names, colSums(group.unique.idx)),
                               ref = snvs.metadata$REF[tumour.only.idx][unique.vars.idx],
                               alt = snvs.metadata$ALT[tumour.only.idx][unique.vars.idx],
                               trinuc = trinucs[tumour.only.idx][unique.vars.idx],
                               stringsAsFactors=F)
unique.vars.cat = build_catalogues(unique.vars.table)
invisible(plot_spectrum(unique.vars.cat, pdf_path=OUTPUT$SPECTRA.PDF))


cat("\nSaving generated objects to file ", OUTPUT$RDATA, "...\n", sep="")
save(countries, tree, tumours.tree, phylo.groups.idx, group.unique.idx, group.names, group.ids,
     file=OUTPUT$RDATA)

cat("\nDone\n\n")

