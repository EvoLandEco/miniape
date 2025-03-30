/**
 * @file dist-topo.js
 * @module distTopo
 *
 * This module implements functions for calculating the 
 * “proportional partitions” of a set of phylogenetic trees,
 * porting the R code in prop.part and its utilities.
 *
 * The main user-facing function is `propPart`, which mimics R’s
 * prop.part. It takes one or more phylogenetic tree objects (each having
 * an edge matrix and tip labels), reorders them in postorder, computes 
 * their bipartitions, and returns an object containing:
 *   - partitions: an array of integer arrays (each representing a clade)
 *   - number: an array of counts for each partition
 *   - labels: the tip labels of the tree
 *   - class: the string "prop.part"
 *
 * Note: This module uses the utility function `reorder` from reorder.js.
 */

import { reorder } from './reorder.js';

/**
 * Compare two arrays for equality (element-wise).
 *
 * @param {Array<number>} a - The first array.
 * @param {Array<number>} b - The second array.
 * @returns {boolean} True if arrays are of the same length and have equal elements.
 */
function arrayEquals(a, b) {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) return false;
  }
  return true;
}

/**
 * Computes the bipartitions (clades) for a tree based on its edge matrix.
 *
 * @param {Array<Array<number>>} orig - The edge matrix (each row is [parent, child])
 *                                      with nodes numbered in a 1-indexed fashion.
 * @param {number} nTips - The number of tip nodes.
 * @returns {Array<Array<number>>} An array of integer arrays, each representing a clade.
 */
function bipartition2(orig, nTips) {
  // Extract parent and child vectors.
  const parent = orig.map(row => row[0]);
  const children = orig.map(row => row[1]);
  const m = Math.max(...parent);
  const nnode = m - nTips;
  // Create an array of empty arrays for the partitions.
  const out = new Array(nnode).fill(null).map(() => []);
  
  // For each edge, add the descendant (or its partition if internal)
  // to the partition of the parent.
  for (let i = 0; i < parent.length; i++) {
    const j = parent[i] - nTips - 1; // 0-indexed position in out
    if (children[i] > nTips) {
      // Child is an internal node: concatenate its partition.
      const y = out[children[i] - nTips - 1];
      out[j] = out[j].concat(y);
    } else {
      // Child is a tip: add its number.
      out[j].push(children[i]);
    }
  }
  // Sort each partition.
  for (let i = 0; i < nnode; i++) {
    out[i].sort((a, b) => a - b);
  }
  return out;
}

/**
 * Internal function to compute the proportional partitions from a list of trees.
 *
 * @param {Array<object>} trees - An array of phylogenetic tree objects.
 *                                Each tree is expected to have an `edge` property.
 * @param {number} nTips - The number of tip nodes.
 * @returns {object} An object with properties:
 *                   - partitions: Array of clade arrays.
 *                   - number: Array of counts corresponding to each clade.
 */
function propPart2(trees, nTips) {
  const nbtree = trees.length;
  // Use the first tree to initialize partitions.
  const M = trees[0];
  const E = M.edge; // edge matrix of the first tree
  let ans = bipartition2(E, nTips); // ans is an array of arrays (clades)
  // Initialize count array: one for each partition.
  let no = ans.map(() => 1);
  if (no.length > 0) {
    no[0] = nbtree; // first partition count becomes the number of trees
  }
  // Process the remaining trees.
  for (let k = 1; k < nbtree; k++) {
    const tmpTree = trees[k];
    const tmpE = tmpTree.edge;
    const bp = bipartition2(tmpE, nTips);
    // For each partition (skip index 0) in bp, compare with those in ans.
    for (let i = 1; i < bp.length; i++) {
      let found = false;
      // Search for an equivalent partition in ans (starting from index 1).
      for (let j = 1; j < ans.length; j++) {
        if (arrayEquals(bp[i], ans[j])) {
          no[j]++;
          found = true;
          break;
        }
      }
      // If not found, add it to the list.
      if (!found) {
        ans.push(bp[i]);
        no.push(1);
      }
    }
  }
  // Return an object with attributes.
  return {
    partitions: ans,
    number: no,
    class: "prop.part"
  };
}

/**
 * Computes the proportional partitions (prop.part) for one or more phylogenetic trees.
 *
 * The function reorders the tree(s) in "postorder", computes the clades
 * (bipartitions) using `propPart2`, and attaches the tip labels as the "labels" attribute.
 *
 * @param {object|Array<object>} tre - A phylogenetic tree object or an array of such objects.
 *                                     Each tree is expected to have properties:
 *                                     - edge: 2D array of edges ([parent, child])
 *                                     - tipLabel: array of tip labels.
 * @returns {object} An object with the following properties:
 *                   - partitions: Array of clade arrays.
 *                   - number: Array of counts.
 *                   - labels: The tip labels of the tree.
 *                   - class: "prop.part"
 */
function propPart(tre) {
  // Check if tre is an array of trees.
  let trees = Array.isArray(tre) ? tre.slice() : [tre];
  // Reorder each tree in postorder.
  trees = trees.map(t => reorder(t, "postorder"));
  const nTips = trees[0].tipLabel.length;
  const clades = propPart2(trees, nTips);
  // Attach the tip labels.
  clades.labels = trees[0].tipLabel;
  return clades;
}

export { propPart };
