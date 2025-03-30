/**
 * @file newick-writer.js
 * @module newickWriter
 *
 * This module provides a function to convert a phylo object (as produced by the
 * Newick parser and other tree manipulations) back to a Newick format string.
 *
 * The phylo object is assumed to have the following properties:
 *   - edge:      2D array of edges, where each row is [parent, child] (nodes are 1-indexed)
 *   - edgeLength: (optional) array of branch lengths corresponding to each edge
 *   - tipLabel:  Array of tip (leaf) labels (for nodes 1..Ntip)
 *   - Nnode:     Number of internal nodes (internal nodes are numbered Ntip+1 to Ntip+Nnode)
 *   - nodeLabel: (optional) Array of internal node labels (indexed by node - Ntip - 1)
 *   - rootEdge:  (optional) Branch length for the root edge (ignored in Newick output)
 *
 * The returned Newick string will include branch lengths (if provided) and internal node labels.
 */

/**
 * Converts a phylo object back into a Newick string.
 *
 * @param {object} phy - The phylo object.
 * @returns {string} The Newick format string representing the tree.
 */
export function writeNewick(phy) {
    const Ntip = phy.tipLabel.length;
    const totalNodes = Ntip + phy.Nnode;
  
    // Build a children map: maps a parent node to an array of { child, length } objects.
    const childrenMap = {};
    for (let i = 0; i < phy.edge.length; i++) {
      const [p, c] = phy.edge[i];
      const len = phy.edgeLength && phy.edgeLength[i] != null ? phy.edgeLength[i] : null;
      if (!childrenMap[p]) childrenMap[p] = [];
      childrenMap[p].push({ child: c, length: len });
    }
  
    // Identify the root: a node that never appears as a child.
    const childSet = new Set(phy.edge.map(row => row[1]));
    let root = null;
    for (let i = 1; i <= totalNodes; i++) {
      if (!childSet.has(i)) {
        root = i;
        break;
      }
    }
    if (root === null) {
      throw new Error("Could not identify the root of the tree.");
    }
  
    /**
     * Recursively builds the Newick string for a given node.
     *
     * @param {number} node - The current node number.
     * @returns {string} Newick string for the subtree rooted at this node.
     */
    function recurse(node) {
      let subtreeStr = "";
      // Check if node is an internal node (has children)
      if (childrenMap[node] && childrenMap[node].length > 0) {
        const childStrs = childrenMap[node].map(childObj => {
          // Recursively process each child and append branch length if available.
          const childNewick = recurse(childObj.child);
          return childNewick + (childObj.length != null ? ":" + childObj.length : "");
        });
        subtreeStr = "(" + childStrs.join(",") + ")";
        // Append internal node label if available.
        if (phy.nodeLabel && (node - Ntip - 1) < phy.nodeLabel.length && phy.nodeLabel[node - Ntip - 1] !== "") {
          subtreeStr += phy.nodeLabel[node - Ntip - 1];
        }
      } else {
        // Leaf node.
        subtreeStr = phy.tipLabel[node - 1];
      }
      return subtreeStr;
    }
  
    // Build the Newick string from the root.
    // (The root edge (if any) is ignored in the output.)
    const newickStr = recurse(root) + ";";
    return newickStr;
  }
  