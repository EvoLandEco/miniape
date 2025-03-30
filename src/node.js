/**
 * @file node.js
 * @module node
 *
 * This module implements functions for computing node depths and heights for a
 * phylogenetic tree, porting the corresponding C routines from R's ape package.
 * The functions assume that edge arrays (edge1 and edge2) are 1-indexed.
 */

/**
 * Computes node depths using edge lengths.
 *
 * This function performs a preorder tree traversal (starting from the bottom of the edge matrix)
 * and updates the node depths in array `xx` using the edge lengths.
 *
 * @param {Array<number>} edge1 - Array of parent node numbers (1-indexed).
 * @param {Array<number>} edge2 - Array of child node numbers (1-indexed).
 * @param {Array<number>} edge_length - Array of branch lengths corresponding to each edge.
 * @param {Array<number>} xx - Array of node depths. Must be preallocated.
 */
function nodeDepthEdgeLength(edge1, edge2, edge_length, xx) {
  const nedge = edge1.length;
  for (let i = nedge - 1; i >= 0; i--) {
    // Convert 1-indexed node numbers to 0-indexed array positions.
    xx[edge2[i] - 1] = xx[edge1[i] - 1] + edge_length[i];
  }
}

/**
 * Computes node depths.
 *
 * Method 1: Node depths are proportional to the number of tips.
 * Method 2: Node depths are evenly spaced.
 *
 * @param {number} ntip - The number of tip nodes.
 * @param {Array<number>} e1 - Array of parent node numbers (1-indexed).
 * @param {Array<number>} e2 - Array of child node numbers (1-indexed).
 * @param {Array<number>} xx - Array of node depths. It will be modified in place.
 * @param {number} method - Calculation method: 1 for proportional to tips, 2 for evenly spaced.
 */
function nodeDepth(ntip, e1, e2, xx, method) {
  const nedge = e1.length;
  
  // Initialize the depths for all tip nodes.
  for (let i = 0; i < ntip; i++) {
    xx[i] = 1;
  }
  
  if (method === 1) {
    // Sum descendant depths into the ancestor.
    for (let i = 0; i < nedge; i++) {
      xx[e1[i] - 1] += xx[e2[i] - 1];
    }
  } else { // method === 2: evenly spaced
    for (let i = 0; i < nedge; i++) {
      // If a value > 0 is already assigned and is at least one more than the descendant's depth, skip.
      if (xx[e1[i] - 1] && xx[e1[i] - 1] >= xx[e2[i] - 1] + 1) continue;
      xx[e1[i] - 1] = xx[e2[i] - 1] + 1;
    }
  }
}

/**
 * Computes node heights.
 *
 * The function assumes that the coordinates (heights) for the tip nodes have already been computed.
 * It processes the edge matrix in pruningwise order.
 *
 * @param {Array<number>} edge1 - Array of parent node numbers (1-indexed).
 * @param {Array<number>} edge2 - Array of child node numbers (1-indexed).
 * @param {Array<number>} yy - Array of node heights. It will be modified in place.
 */
function nodeHeight(edge1, edge2, yy) {
  const nedge = edge1.length;
  let S = 0;
  let n = 0;
  let i;
  
  // Process all edges except the last one.
  for (i = 0; i < nedge - 1; i++) {
    S += yy[edge2[i] - 1];
    n++;
    // When the next edge has a different parent, average the accumulated heights.
    if (edge1[i + 1] !== edge1[i]) {
      yy[edge1[i] - 1] = S / n;
      S = 0;
      n = 0;
    }
  }
  // Process the last edge.
  S += yy[edge2[i] - 1];
  n++;
  yy[edge1[i] - 1] = S / n;
}

/**
 * Computes node heights in a cladewise manner.
 *
 * First, node depths are computed using method 1 (proportional to number of tips),
 * and then node heights are computed as a weighted average based on the node depths.
 *
 * @param {number} ntip - The number of tip nodes.
 * @param {Array<number>} edge1 - Array of parent node numbers (1-indexed).
 * @param {Array<number>} edge2 - Array of child node numbers (1-indexed).
 * @param {Array<number>} xx - Array for node depths. Will be computed by this function.
 * @param {Array<number>} yy - Array of node heights. It will be modified in place.
 */
function nodeHeightClado(ntip, edge1, edge2, xx, yy) {
  // Compute node depths with method 1.
  nodeDepth(ntip, edge1, e2, xx, 1);
  
  const nedge = edge1.length;
  let S = 0;
  let n = 0;
  let i, j;
  
  for (i = 0; i < nedge - 1; i++) {
    j = edge2[i] - 1;
    S += yy[j] * xx[j];
    n += xx[j];
    if (edge1[i + 1] !== edge1[i]) {
      yy[edge1[i] - 1] = S / n;
      S = 0;
      n = 0;
    }
  }
  // Process the last edge.
  j = edge2[i] - 1;
  S += yy[j] * xx[j];
  n += xx[j];
  yy[edge1[i] - 1] = S / n;
}

export {
  nodeDepthEdgeLength,
  nodeDepth,
  nodeHeight,
  nodeHeightClado
};
