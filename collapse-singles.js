/**
 * @file collapse-singles.js
 * @module collapseSingles
 *
 * This module provides two functions:
 *   - hasSingles(tree): Tests if a tree (or array of trees) has any single-child (singleton) nodes.
 *   - collapseSingles(tree, rootEdge = false): Collapses single-child nodes in a phylogenetic tree,
 *     mimicking the behavior of the collapse.singles function from the R ape package.
 *
 * The algorithm uses a series of steps:
 *   1. Reorder the tree (assumed to be available via reorder()).
 *   2. Tabulate parent frequencies from the edge matrix.
 *   3. If no internal node (nodes numbered > number of tips) appears only once, return the tree unchanged.
 *   4. If branch lengths are available, update them accordingly.
 *   5. Remove basal (root) edges that are singletons.
 *   6. Collapse any remaining singleton edges by reassigning the descendant of the parent edge.
 *   7. Renumber nodes to conform to the phylo format.
 *
 * Note: This implementation assumes 1-indexed node numbers (as in R's ape).
 */

import { reorder } from './reorder.js';
import { tabulate } from './common-utilities.js';

/**
 * Determines whether a phylogenetic tree (or array of trees) has any single-child nodes.
 *
 * @param {object|Array<object>} tree - A phylogenetic tree object or an array of such objects.
 *                                      Each tree is expected to have an `edge` property
 *                                      (a 2D array where each row is [parent, child]).
 * @returns {boolean|Array<boolean>} True if at least one singleton is found; if multiple trees are provided,
 *                                   returns an array of booleans.
 */
export function hasSingles(tree) {
  const fun = (x) => {
    const e1 = x.edge.map(row => row[0]);
    const tab = tabulate(e1);
    // If any count equals 1, return true.
    for (const count of tab) {
      if (count === 1) return true;
    }
    return false;
  };
  if (!Array.isArray(tree)) return fun(tree);
  return tree.map(fun);
}

/**
 * Collapses single-child (singleton) nodes in a phylogenetic tree.
 *
 * @param {object} tree - A phylogenetic tree object in "phylo" format.
 *                        Expected properties:
 *                          - tipLabel: array of tip labels.
 *                          - edge: 2D array (each row: [parent, child], 1-indexed).
 *                          - edgeLength: (optional) array of branch lengths.
 *                          - nodeLabel: (optional) array of node labels.
 *                          - Nnode: number of internal nodes.
 * @param {boolean} [rootEdge=false] - If TRUE and branch lengths exist, accumulates branch length for the root.
 * @returns {object} The tree with singleton nodes collapsed.
 */
export function collapseSingles(tree, rootEdge = false) {
  const n = tree.tipLabel.length;
  if (n === 0) return tree;

  // Reorder the tree.
  tree = reorder(tree);
  // Extract parent and child columns from the edge matrix.
  let e1 = tree.edge.map(row => row[0]);
  let e2 = tree.edge.map(row => row[1]);

  // Compute frequency of each parent.
  let tab = tabulate(e1);
  // Check if all internal nodes (indices > n) have counts > 1.
  let allInternalMultiple = true;
  for (let i = n + 1; i < tab.length; i++) {
    if (tab[i] <= 1) {
      allInternalMultiple = false;
      break;
    }
  }
  if (allInternalMultiple) return tree;

  // Determine if branch lengths are available.
  let wbl = false;
  let el = [];
  if (tree.edgeLength == null) {
    rootEdge = false;
  } else {
    wbl = true;
    el = tree.edgeLength.slice();
  }

  let ROOTEDGE = 0;
  if (rootEdge) ROOTEDGE = 0;

  // Start with the root node (by convention, n + 1).
  let ROOT = n + 1;
  // Remove basal singletons: while the count for ROOT is exactly 1,
  // find the unique edge with parent == ROOT, update ROOT to its child, and remove that edge.
  while (true) {
    tab = tabulate(e1);
    if (tab[ROOT] === 1) {
      // Find index where e1 equals ROOT.
      const iIndices = [];
      for (let i = 0; i < e1.length; i++) {
        if (e1[i] === ROOT) iIndices.push(i);
      }
      if (iIndices.length === 0) break;
      const iIndex = iIndices[0];
      ROOT = e2[iIndex];
      e1.splice(iIndex, 1);
      e2.splice(iIndex, 1);
      if (wbl) {
        if (rootEdge) ROOTEDGE += el[iIndex];
        el.splice(iIndex, 1);
      }
    } else {
      break;
    }
  }

  // Identify singleton internal nodes.
  const tabE1 = tabulate(e1);
  const singles = [];
  for (let i = 0; i < tabE1.length; i++) {
    if (i > n && tabE1[i] === 1) {
      singles.push(i);
    }
  }
  if (singles.length > 0) {
    // For each singleton value, find its first occurrence in e1.
    const ii = singles.map(s => e1.indexOf(s)).sort((a, b) => b - a);
    // For each index in ii, find the first index in e2 matching e1[ii].
    const jj = ii.map(idx => e2.indexOf(e1[idx]));
    for (let k = 0; k < ii.length; k++) {
      e2[jj[k]] = e2[ii[k]];
      if (wbl) {
        el[jj[k]] = el[jj[k]] + el[ii[k]];
      }
    }
    // Remove the edges at indices in ii.
    for (const idx of ii) {
      e1.splice(idx, 1);
      e2.splice(idx, 1);
      if (wbl) el.splice(idx, 1);
    }
  }

  // Update number of internal nodes.
  const Nnode = e1.length - n + 1;

  // Determine unique internal nodes from the (updated) parent column.
  const oldnodes = [...new Set(e1)];
  if (tree.nodeLabel) {
    // Update node labels: use the (old node - n) index.
    tree.nodeLabel = oldnodes.map(val => tree.nodeLabel[val - n - 1]);
  }

  // Create a new numbering for internal nodes.
  const maxOld = Math.max(...oldnodes);
  const newNb = new Array(maxOld).fill(0);
  newNb[ROOT - 1] = n + 1;

  // For all nodes in e2 that are internal (i.e. > n), assign new numbers.
  const internalNodes = [...new Set(e2.filter(val => val > n))].sort((a, b) => a - b);
  const newNumbers = [];
  for (let num = n + 2; num <= n + Nnode; num++) {
    newNumbers.push(num);
  }
  for (let i = 0; i < internalNodes.length; i++) {
    newNb[internalNodes[i] - 1] = newNumbers[i];
  }
  // Update e2 values using the new numbering.
  e2 = e2.map(val => (val > n ? newNb[val - 1] : val));
  // Update e1 values as well.
  e1 = e1.map(val => newNb[val - 1]);

  // Rebuild the edge matrix.
  tree.edge = [];
  for (let i = 0; i < e1.length; i++) {
    tree.edge.push([e1[i], e2[i]]);
  }
  tree.Nnode = Nnode;
  if (wbl) {
    if (rootEdge) tree.rootEdge = ROOTEDGE;
    tree.edgeLength = el;
  }
  return tree;
}
