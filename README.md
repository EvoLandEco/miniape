# miniape

**miniape** is a lightweight JavaScript library for phylogenetic tree manipulation, **ported from** the [APE](https://cran.r-project.org/package=ape) package in R.

It supports a simplified `phylo` object format and provides essential utilities for reading, transforming, and exporting phylogenetic trees in Newick format.

---

## Features

- 📖 Parse Newick strings into JavaScript `phylo` objects
- ✂️ Drop one or more tips from a tree (`dropTip`)
- 🌱 Root and unroot trees (`rootPhylo`, `unrootPhylo`)
- 🔁 Reorder trees (`reorder`)
- 🔽 Collapse internal single-child nodes (`collapseSingles`)
- 🌐 Extract topology partitions (`propPart`)
- 📝 Write `phylo` objects back to Newick format (`writeNewick`)

---

## Installation

Clone this repository:

```bash
git clone https://github.com/EvoLandEco/miniape.git
