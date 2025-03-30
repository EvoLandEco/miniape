/**
 * @file drop-tip.js
 * @module dropTip
 *
 * This module exports the function `dropTip` which removes one or more tips
 * from a phylogenetic tree (in "phylo" format), mimicking the behavior of R's
 * drop.tip.phylo function from the ape package.
 *
 * Note: This implementation relies on the following utilities:
 *   - reorder(phy, order?)        (from reorder.js)
 *   - isRooted(phy)               (from root.js)
 *   - collapseSingles(phy)        (from collapse-singles.js)
 *   - nodeDepth(ntip, e1, e2, xx, method) (from node.js)
 *   - rootPhylo(phy, tip)              (from root.js)
 *
 */

/* Imports from other modules */
import { reorder } from "./reorder.js";
import { isRooted } from "./root.js";
import { collapseSingles } from "./collapse-singles.js";
import { nodeDepth } from "./node.js";
import { rank } from "./common-utilities.js";
import { rootPhylo } from "./root.js";
import { tabulate } from "./common-utilities.js";

/**
 * Removes tips from a phylogenetic tree.
 *
 * This function mimics the behavior of R's drop.tip.phylo.
 *
 * @param {object} phy - A phylogenetic tree object. Expected properties:
 *   - tipLabel: Array of tip labels.
 *   - edge: 2D array of edges ([parent, child]) with 1-indexed node numbers.
 *   - edgeLength: (optional) Array of branch lengths.
 *   - nodeLabel: (optional) Array of node labels.
 *   - Nnode: Number of internal nodes.
 * @param {number|string|Array<number|string>} tip - Tip(s) to drop. If tip labels are provided,
 *   they are converted to tip indices.
 * @param {boolean} [trimInternal=true] - Whether to trim internal edges that lose descendants.
 * @param {boolean} [subtree=false] - If TRUE, keep the subtree subtending dropped tips.
 * @param {number} [rootEdge=0] - Used when adjusting the root edge.
 * @param {boolean} [rooted=isRooted(phy)] - Whether the tree is rooted.
 * @param {boolean} [collapseSinglesFlag=true] - Whether to collapse single-child nodes after dropping.
 * @returns {object|null} The pruned tree, or null if all tips are dropped.
 * @throws {Error} If a required utility (e.g. root) is missing.
 */
export function dropTip(
  phy,
  tip,
  trimInternal = true,
  subtree = false,
  rootEdge = 0,
  rooted = isRooted(phy),
  collapseSinglesFlag = true
) {
  const Ntip = phy.tipLabel.length;

  // If tip is a string, convert it to an array.
  if (typeof tip === "string") tip = [tip];
  // If tip is provided as tip labels (array of strings), convert to indices (1-indexed).
  if (Array.isArray(tip) && typeof tip[0] === "string") {
    tip = tip
      .map((t) => {
        const idx = phy.tipLabel.indexOf(t);
        return idx !== -1 ? idx + 1 : null;
      })
      .filter((x) => x !== null);
  }

  // Remove tip numbers that are out of range.
  const outOfRange = tip.filter((t) => t > Ntip);
  if (outOfRange.length > 0) {
    console.warn(
      "some tip numbers were larger than the number of tips: they were ignored"
    );
    tip = tip.filter((t) => t <= Ntip);
  }

  if (tip.length === 0) return phy;

  // If dropping all tips.
  if (tip.length === Ntip) {
    if (phy.Nnode < 3 || trimInternal) {
      console.warn("drop all tips of the tree: returning null");
      return null;
    }
  }

  const wbl = phy.edgeLength != null;

  // Special case: if dropping all but one tip and trimInternal is TRUE.
  if (tip.length === Ntip - 1 && trimInternal) {
    // Identify the single remaining tip.
    const allTips = Array.from({ length: Ntip }, (_, i) => i + 1);
    const remaining = allTips.filter((x) => !tip.includes(x));
    const iEdge = phy.edge.findIndex((row) => row[1] === remaining[0]);
    if (iEdge === -1) {
      console.warn("Edge for the remaining tip not found");
      return phy;
    }
    const res = {
      edge: [[2, 1]], // mimics matrix(2:1, 1, 2)
      tipLabel: [phy.tipLabel[phy.edge[iEdge][1] - 1]],
      Nnode: 1,
    };
    if (wbl) res.edgeLength = [phy.edgeLength[iEdge]];
    if (phy.nodeLabel)
      res.nodeLabel = [phy.nodeLabel[phy.edge[iEdge][0] - Ntip - 1]];
    return res;
  }

  // If tree is unrooted.
  if (!rooted) {
    phy.rootEdge = null;
    if (subtree) {
      // Root the tree using the first remaining tip.
      const allTips = Array.from({ length: Ntip }, (_, i) => i + 1);
      const remaining = allTips.filter((x) => !tip.includes(x));
      phy = rootPhylo(phy, remaining[0]);
      rootEdge = 0;
    }
  }

  phy = reorder(phy);

  let NEWROOT = Ntip + 1,
    ROOT = Ntip + 1;
  const Nnode = phy.Nnode;
  const Nedge = phy.edge.length;
  let Nvec = null;
  if (subtree) {
    trimInternal = true;
    const tr = reorder(phy, "postorder");
    // Allocate array for node depths.
    Nvec = new Array(Ntip + Nnode).fill(0);
    // Compute node depths (method 1: proportional to the number of tips).
    const tr_e1 = tr.edge.map((row) => row[0]);
    const tr_e2 = tr.edge.map((row) => row[1]);
    nodeDepth(Ntip, tr_e1, tr_e2, Nvec, 1);
  }

  // Create local copies of edge columns.
  const edge1 = phy.edge.map((row) => row[0]);
  const edge2 = phy.edge.map((row) => row[1]);
  const keep = new Array(Nedge).fill(true);

  // Delete terminal edges corresponding to the tips to drop.
  tip.forEach((t) => {
    const idx = edge2.indexOf(t);
    if (idx !== -1) keep[idx] = false;
  });

  if (trimInternal) {
    const ints = edge2.map((x) => x > Ntip);
    // Repeat: remove internal edges that lose all descendants.
    while (true) {
      const allowed = new Set();
      for (let i = 0; i < Nedge; i++) {
        if (keep[i]) allowed.add(edge1[i]);
      }
      const sel = [];
      for (let i = 0; i < Nedge; i++) {
        if (keep[i] && ints[i] && !allowed.has(edge2[i])) {
          sel.push(i);
        }
      }
      if (sel.length === 0) break;
      sel.forEach((i) => {
        keep[i] = false;
      });
    }
    if (subtree) {
      // Keep the subtending edge(s).
      const setKeep = new Set(edge1.filter((v, i) => keep[i]));
      const setNotKeep = new Set(edge1.filter((v, i) => !keep[i]));
      for (let i = 0; i < Nedge; i++) {
        if (setKeep.has(edge1[i]) && setNotKeep.has(edge1[i])) {
          keep[i] = true;
        }
      }
    }
    if (rootEdge && wbl) {
      const keptEdge1 = edge1.filter((v, i) => keep[i]);
      const degree = tabulate(keptEdge1);
      if (degree[ROOT] === 1) {
        let j = [];
        while (true) {
          const indices = [];
          for (let i = 0; i < Nedge; i++) {
            if (keep[i] && edge1[i] === NEWROOT) {
              indices.push(i);
            }
          }
          if (indices.length === 0) break;
          j = indices.concat(j);
          NEWROOT = edge2[indices[0]];
          if ((degree[NEWROOT] || 0) > 1) break;
        }
        j.forEach((idx) => {
          keep[idx] = false;
        });
        j = j.slice(0, rootEdge);
        let NewRootEdge = 0;
        j.forEach((idx) => {
          NewRootEdge += phy.edgeLength ? phy.edgeLength[idx] : 0;
        });
        if (j.length < rootEdge && phy.rootEdge != null)
          NewRootEdge += phy.rootEdge;
        phy.rootEdge = NewRootEdge;
      }
    }
  }

  // Drop the edges not kept.
  const newEdge = [];
  const newEdgeLength = wbl ? [] : null;
  for (let i = 0; i < Nedge; i++) {
    if (keep[i]) {
      newEdge.push(phy.edge[i]);
      if (wbl) newEdgeLength.push(phy.edgeLength[i]);
    }
  }
  phy.edge = newEdge;
  if (wbl) phy.edgeLength = newEdgeLength;

  // Identify new terminal edges: those whose descendant is not a parent.
  const parentSet = new Set(phy.edge.map((row) => row[0]));
  const TERMS = phy.edge.map((row) => !parentSet.has(row[1]));
  // Get the old node numbers of those edges (column 2).
  const oldNoOfNewTips = [];
  for (let i = 0; i < phy.edge.length; i++) {
    if (TERMS[i]) {
      oldNoOfNewTips.push(phy.edge[i][1]);
    }
  }

  // If subtree is TRUE, mark tips dropped but kept due to subtree = true.
  if (subtree) {
    const indicesToRemove = [];
    tip.forEach((t, index) => {
      if (oldNoOfNewTips.includes(t)) {
        phy.tipLabel[t - 1] = "[1_tip]";
        indicesToRemove.push(index);
      }
    });
    tip = tip.filter((_, index) => !indicesToRemove.includes(index));
  }

  const newNtip = oldNoOfNewTips.length;

  // Replace terminal edge numbers with their rank.
  const termEdges = phy.edge
    .filter((row, i) => TERMS[i])
    .map((row) => row[1]);
  const termRanks = rank(termEdges);
  let tIdx = 0;
  for (let i = 0; i < phy.edge.length; i++) {
    if (TERMS[i]) {
      phy.edge[i][1] = termRanks[tIdx];
      tIdx++;
    }
  }

  // Remove tip labels corresponding to dropped tips.
  if (tip.length > 0) {
    tip.sort((a, b) => a - b);
    for (let i = tip.length - 1; i >= 0; i--) {
      phy.tipLabel.splice(tip[i] - 1, 1);
    }
  }

  // Make new tip labels if necessary.
  if (subtree || !trimInternal) {
    const node2tip = oldNoOfNewTips.filter((x) => x > Ntip);
    let newTipLabel;
    if (node2tip.length === 0) {
      newTipLabel = [];
    } else if (subtree && !phy.nodeLabel) {
      newTipLabel = node2tip.map(
        (x) => "[" + (Nvec ? Nvec[x - 1] : "?") + "_tips]"
      );
    } else {
      newTipLabel = !phy.nodeLabel
        ? new Array(node2tip.length).fill("NA")
        : node2tip.map((x) => phy.nodeLabel[x - Ntip - 1]);
    }
    phy.tipLabel = phy.tipLabel.concat(newTipLabel);
  }

  // Update the number of internal nodes.
  phy.Nnode = phy.edge.length - newNtip + 1;

  // newNtip is the new tip count (after dropping)
  const newNb = new Array(Ntip + Nnode).fill(0);

  // In ape code: newNb[NEWROOT] = n + 1
  // Here zero-based indexing is used: 
  newNb[NEWROOT - 1] = newNtip + 1;

  // Identify the edges whose child is an internal node (child > newNtip)
  const sndcol = [];
  for (let i = 0; i < phy.edge.length; i++) {
    if (phy.edge[i][1] > newNtip) {
      sndcol.push(i);
    }
  }

  // Gather the set of child nodes that are internal
  let nodesToRenumber = [];
  for (const i of sndcol) {
    nodesToRenumber.push(phy.edge[i][1]);
  }
  nodesToRenumber = Array.from(new Set(nodesToRenumber)).sort((a, b) => a - b);

  // In ape: newNb[sort(...) ] = (n+2):(n+phy$Nnode)
  nodesToRenumber.forEach((nodeVal, i) => {
    // nodeVal is the old node number, 1-based
    // newNtip+2 + i is the new node number
    newNb[nodeVal - 1] = newNtip + 2 + i;
  });

  // First renumber the child column
  for (let i = 0; i < sndcol.length; i++) {
    const idx = sndcol[i];
    const oldChild = phy.edge[idx][1];
    phy.edge[idx][1] = newNb[oldChild - 1];
  }

  // Then renumber all parents
  for (let i = 0; i < phy.edge.length; i++) {
    const oldParent = phy.edge[i][0];
    phy.edge[i][0] = newNb[oldParent - 1];
  }

  // Convert any undefined or fractional values to integers (shouldn't occur, but just in case)
  phy.edge = phy.edge.map(row => row.map(x => (x ? Math.floor(x) : 0)));

  // In ape code:
  //   if (!is.null(phy$node.label)) phy$node.label <- phy$node.label[which(newNb > 0) - oldNtip]

  // Update nodeLabel if it exists
  if (phy.nodeLabel) {
    // Assume oldNtip = the original number of tips before dropping. 
    // Originally ape's logic: "which(newNb>0) - Ntip"
    const indices = [];
    for (let i = 0; i < newNb.length; i++) {
      // newNb[i] > 0 indicates the old node i+1 still exists
      // and i >= Ntip means it's an internal node
      if (newNb[i] > 0 && i >= Ntip) {
        indices.push(i - Ntip);
      }
    }
    phy.nodeLabel = indices.map(i => phy.nodeLabel[i]);
  }

  if (collapseSinglesFlag) {
    if (typeof collapseSingles === "function") {
      phy = collapseSingles(phy);
    } else {
      throw new Error("Utility function 'collapseSingles' is missing.");
    }
  }

  return phy;
}

