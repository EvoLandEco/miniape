dummy = read.tree("./test/dummy.tre")

plot(dummy)

nodelabels()
tiplabels()
edgelabels()

edgelabels(text = round(dummy$edge.length, 2), cex = 0.8, bg = "lightyellow", adj = c(0.5, -3), frame = "rect", col = "blue")

# reorder the tree
dummy = reorder(dummy, "cladewise")
View(dummy)
plot(dummy)
nodelabels()
tiplabels()
edgelabels()
edgelabels(text = round(dummy$edge.length, 2), cex = 0.8, bg = "lightyellow", adj = c(0.5, -3), frame = "rect", col = "blue")


# read a complex tree

complex_tree = read.tree("./complex_tree.tre")
plot(complex_tree)

# Drop tips "B", "D", and "Q" from the tree.
pruned_complex = drop.tip(complex_tree, c("B", "D", "Q"))
plot(pruned_complex)

js_pruned_complex = read.tree("./pruned_complex_tree.tre")
plot(js_pruned_complex)

# Plot reference and js pruned trees side by side
par(mfrow = c(1, 2))

js_pruned = read.tree("./test/out/pruned_complex_tree.tre")

plot(js_pruned, main = "js pruned")
# Add labels to the nodes
nodelabels()
tiplabels()
edgelabels()
edgelabels(text = round(pruned_complex$edge.length, 2), cex = 0.8, bg = "lightyellow", adj = c(0.5, -3), frame = "rect", col = "blue")

plot(pruned_complex, main = "R pruned")
# Add labels to the nodes
nodelabels()
tiplabels()
edgelabels()
edgelabels(text = round(pruned_complex$edge.length, 2), cex = 0.8, bg = "lightyellow", adj = c(0.5, -3), frame = "rect", col = "blue")

# Use ape to generate a birth-death tree, ultrametric.
big_tree = rlineage(
  Tmax = 10,
  birth = 0.5,
  death = 0.01
)

big_tree = drop.fossil(big_tree)

# write the tree to newick
write.tree(big_tree, file = "./test/big_tree.tre")

# read the tree pruned by js
js_pruned_big_tree = read.tree("./test/out/pruned_big_tree.tre")

# read the txt list of dropped tips
dropped_tips = readLines("./test/out/tips_to_drop.txt")

# Drop these tips from the tree
pruned_big_tree = drop.tip(big_tree, dropped_tips)

# Plot reference and js pruned trees side by side
par(mfrow = c(1, 2))
plot(js_pruned_big_tree, main = "js pruned")
plot(pruned_big_tree, main = "R pruned")

# Check if tip labels are the same
identical(js_pruned_big_tree$tip.label, pruned_big_tree$tip.label)

# Check if edge lengths are the same
all.equal(js_pruned_big_tree$edge.length, pruned_big_tree$edge.length)

# Check if the tree topology is the same
identical(js_pruned_big_tree$edge, pruned_big_tree$edge)