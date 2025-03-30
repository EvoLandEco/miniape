/**
 * Port of ape's reorderRcpp from C++ into JavaScript.
 *
 * @param {Array<Array<number>>} orig - The original edge matrix (0-based),
 *        where each row is [parent, child].
 * @param {number} nTips - The number of tips (1..nTips).
 * @param {number} root - The root node number (usually nTips + 1).
 * @param {number} order - 1 for cladewise, 2 for postorder.
 * @returns {Array<number>} A 0-based integer array giving the new edge order.
 */
function reorderRcpp(orig, nTips, root, order) {
    // Extract edge vectors e1 and e2 from the two-dimensional array "orig"
    const n = orig.length;
    const e1 = new Array(n);
    const e2 = new Array(n);
    for (let i = 0; i < n; i++) {
      e1[i] = orig[i][0];
      e2[i] = orig[i][1];
    }
  
    // Determine the maximum value in e1 and the number of internal nodes.
    const m = Math.max(...e1);
    const nnode = m - nTips;

    // Initialize arrays
    const L = new Array(n).fill(0);
    const neworder = new Array(n).fill(0);
    const pos = new Array(nnode).fill(0);
    const xi = new Array(nnode).fill(0);
    const xj = new Array(nnode).fill(0);
  
    // Build the xj array: count how many times each internal node appears in e1.
    for (let i = 0; i < n; i++) {
      const idx = e1[i] - nTips - 1; // convert node label to 0-index for internal nodes
      xj[idx]++;
    }
  
    // Build the xi array: cumulative sum of xj (xi[0] is assumed to be 0).
    for (let i = 1; i < nnode; i++) {
      xi[i] = xi[i - 1] + xj[i - 1];
    }
  
    // Populate the L array using positions from xi and pos.
    for (let i = 0; i < n; i++) {
      const k = e1[i] - nTips - 1;
      const j = pos[k];
      L[xi[k] + j] = i;
      pos[k]++;
    }
  
    // "iii" will be the running index into neworder.
    let iii;
  
    // Recursive function corresponding to foo_reorderRcpp in C++
    function foo_reorderRcpp(node) {
      const i = node - nTips - 1; // convert node to index in the arrays
      for (let j = 0; j < xj[i]; j++) {
        const k = L[xi[i] + j];
        neworder[iii] = k;  // store k in neworder (as in C++ code)
        iii++;
        // If the child node is an internal node, recurse.
        if (e2[k] > nTips) {
          foo_reorderRcpp(e2[k]);
        }
      }
    }
  
    // Recursive function corresponding to bar_reorderRcpp in C++
    function bar_reorderRcpp(node) {
      const i = node - nTips - 1;
      // Process in reverse order.
      for (let j = xj[i] - 1; j >= 0; j--) {
        neworder[iii] = L[xi[i] + j] + 1; // Note: ot sure if +1 is needed here
        iii--;
      }
      // Recurse on children that are internal nodes.
      for (let j = 0; j < xj[i]; j++) {
        const k = e2[L[xi[i] + j]];
        if (k > nTips) {
          bar_reorderRcpp(k);
        }
      }
    }
  
    // Choose which ordering to use based on the 'order' parameter.
    if (order === 1) {
      iii = 0;
      foo_reorderRcpp(root);
    } else if (order === 2) {
      iii = n - 1;
      bar_reorderRcpp(root);
    } else {
      throw new Error("order must be 1 or 2");
    }
  
    return neworder;
  }  

/**
 * Internal function to reorder a phylogenetic tree object.
 *
 * @param {object} tree - The phylogenetic tree object. Expected properties:
 *   - edge: Array of edges (each edge is an array [parent, child])
 *   - edgeLength: Array of branch lengths (optional)
 *   - tipLabel: Array of tip labels
 *   - Nnode: Number of internal nodes
 *   - order: (optional) current order attribute
 * @param {string} order - The output order ("cladewise" or "postorder").
 * @param {boolean} indexOnly - If true, only return the new order indices.
 * @param {number} nbTip - Number of tip nodes.
 * @param {number} io - Order code (1 for "cladewise", 2 for "postorder").
 * @returns {object|Array<number>} The reordered tree object or new order indices if indexOnly is true.
 */
function _reorder_ape(tree, order, indexOnly, nbTip, io) {
    const nbEdge = tree.edge.length;

    // If the tree is already ordered, return early.
    if (tree.order != null && tree.order === order) {
        if (indexOnly) {
            return Array.from({ length: nbEdge }, (_, i) => i);
        }
        return tree;
    }

    // If there is only one internal node, no reordering is needed.
    if (tree.Nnode === 1) {
        if (indexOnly) {
            return Array.from({ length: nbEdge }, (_, i) => i);
        }
        return tree;
    }

    // Compute new order indices using the helper function.
    const neworder = reorderRcpp(tree.edge, nbTip, nbTip + 1, io);

    if (indexOnly) {
        return neworder;
    }

    // Reorder the edge matrix using the new order.
    tree.edge = neworder.map(idx => tree.edge[idx]);

    // Reorder branch lengths if available.
    if (tree.edgeLength != null) {
        tree.edgeLength = neworder.map(idx => tree.edgeLength[idx]);
    }

    // Set the order attribute on the tree.
    tree.order = order;

    return tree;
}

/**
 * Reorders a phylogenetic tree object.
 *
 * This function mimics the behavior of R's reorder.phylo function,
 * but assumes a 0-indexed tree.
 *
 * @param {object} tree - The phylogenetic tree object with expected properties:
 *   - edge: Array of edges ([parent, child])
 *   - edgeLength: Array of branch lengths (optional)
 *   - tipLabel: Array of tip labels
 *   - Nnode: Number of internal nodes
 *   - order: (optional) current order attribute
 * @param {string} [order="cladewise"] - The output order ("cladewise" or "postorder").
 * @param {boolean} [indexOnly=false] - If true, returns only the new order indices.
 * @returns {object|Array<number>} The reordered tree object or the new order indices if indexOnly is true.
 * @throws {Error} If the specified order is not recognized.
 */
function reorder(tree, order = "cladewise", indexOnly = false) {
    const ORDER = ["cladewise", "postorder"];
    const ioIndex = ORDER.indexOf(order);
    if (ioIndex === -1) {
        throw new Error("ambiguous order");
    }
    // Convert to 1-based order code in the original code; here, 1 for "cladewise", 2 for "postorder".
    const io = ioIndex + 1;

    const n = tree.tipLabel.length;
    if (n < 2) {
        return tree;
    } else {
        return _reorder_ape(tree, ORDER[ioIndex], indexOnly, n, io);
    }
}

export { reorder };
