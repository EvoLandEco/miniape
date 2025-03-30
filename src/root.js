/**
 * @file root.js
 * @module root
 *
 * This module provides functionality for handling rooting of phylogenetic trees.
 * It exports:
 *   - isRooted: Checks whether a tree is rooted.
 *   - rootPhylo: Roots a phylogenetic tree using an outgroup or an explicit node.
 *   - unrootPhylo: Unroots a phylogenetic tree.
 *
 * Assumptions:
 *   - Each tree is a "phylo" object with:
 *       edge: 2D array of [parent, child] (1-indexed)
 *       edgeLength: (optional) array of branch lengths
 *       tipLabel: array of tip labels
 *       Nnode: number of internal nodes
 *   - The reorder function is available from reorder.js
 *   - The propPart function is available from dist-topo.js
 */

import { reorder } from "./reorder.js";
import { propPart } from "./dist-topo.js";

/**
 * Checks if a phylogenetic tree is rooted.
 *
 * @param {object} phy - A phylogenetic tree object.
 * @param {number} ntips - Number of tip nodes in the tree.
 * @returns {boolean} True if the tree is rooted, false otherwise.
 */
export function isRooted(phy, ntips) {
  if (phy.rootEdge != null) return true;
  const parentCounts = {};
  phy.edge.forEach(edge => {
    const parent = edge[0];
    parentCounts[parent] = (parentCounts[parent] || 0) + 1;
  });
  const count = parentCounts[ntips + 1] || 0;
  return count <= 2;
}

function arrayEquals(a, b) {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) return false;
  }
  return true;
}

/**
 * Roots a phylogenetic tree.
 *
 * @param {object} phy - A phylogenetic tree object.
 * @param {number|string|Array<number|string>} outgroup - Tip(s) used as outgroup.
 * @param {number|null} [node=null] - If set, use this node as the new root.
 * @param {boolean} [resolveRoot=false] - Whether to resolve the root if ambiguous.
 * @param {boolean} [edgelabel=false] - Whether to transfer edge labels.
 * @returns {object} The re-rooted phylogenetic tree.
 * @throws {Error} If invalid parameters are given or utilities are missing.
 */
export function rootPhylo(
  phy,
  outgroup,
  node = null,
  resolveRoot = false,
  edgelabel = false
) {
  if (!phy || !phy.tipLabel || !phy.edge) {
    throw new Error("object not of class 'phylo'");
  }
  phy = reorder(phy);
  const n = phy.tipLabel.length;
  const ROOT = n + 1;

  if (node === null && Array.isArray(outgroup) && outgroup.length > 1 && resolveRoot) {
    phy = unrootPhylo(phy);
  }
  const e1 = phy.edge.map(row => row[0]);
  const e2 = phy.edge.map(row => row[1]);
  const wbl = phy.edgeLength != null;

  let newroot;
  let MRCA_outgroup;

  if (node !== null) {
    if (node <= n) {
      throw new Error("incorrect node#: should be greater than the number of taxa");
    }
    outgroup = null;
    newroot = node;
  } else {
    if (Array.isArray(outgroup)) {
      if (typeof outgroup[0] === "number") {
        if (outgroup.some(t => t > n)) {
          throw new Error("incorrect taxa#: should not be greater than the number of taxa");
        }
      } else if (typeof outgroup[0] === "string") {
        outgroup = outgroup.map(label => {
          const idx = phy.tipLabel.indexOf(label);
          if (idx === -1) throw new Error("specified outgroup not in labels of the tree");
          return idx + 1;
        });
      }
    } else {
      outgroup = [outgroup];
    }
    if (outgroup.length === n) return phy;
    outgroup = outgroup.sort((a, b) => a - b);

    if (outgroup.length > 1) {
      if (typeof propPart !== "function") {
        throw new Error("Utility function 'propPart' is missing.");
      }
      const pp = propPart(phy).partitions;
      const ingroup = [];
      for (let i = 1; i <= n; i++) {
        if (!outgroup.includes(i)) ingroup.push(i);
      }
      newroot = 0;
      for (let i = 1; i < phy.Nnode; i++) {
        const part = pp[i];
        if (arrayEquals(part, ingroup)) {
          const idx = e2.findIndex(val => val === i + n);
          newroot = e1[idx];
          break;
        }
        if (arrayEquals(part, outgroup)) {
          newroot = i + n;
          break;
        }
      }
      if (!newroot) throw new Error("the specified outgroup is not monophyletic");
      MRCA_outgroup = phy.Nnode + n;
    } else {
      const idx = e2.findIndex(val => val === outgroup[0]);
      newroot = e1[idx];
    }
  }
  const N = phy.edge.length;
  const oldNnode = phy.Nnode;
  const degree = {};
  e1.forEach(val => {
    degree[val] = (degree[val] || 0) + 1;
  });
  const Nclade = degree[ROOT] || 0;
  const fuseRoot = Nclade === 2;

  if (newroot === ROOT) {
    if (!resolveRoot) return phy;
    if (outgroup && outgroup.length > 1) outgroup = [MRCA_outgroup];
    if (node !== null) {
      throw new Error("ambiguous resolution of the root node: specify an outgroup");
    }
    const k = [];
    e1.forEach((val, i) => {
      if (val === ROOT) k.push(i);
    });
    if (k.length > 2) {
      const i_index = e2.findIndex(val => val === outgroup[0]);
      const j = k.filter(idx => idx !== i_index);
      const newnod = oldNnode + n + 1;
      j.forEach(idx => {
        phy.edge[idx][0] = newnod;
      });
      phy.edge = [[ROOT, newnod]].concat(phy.edge);
      if (wbl) {
        phy.edgeLength = [0].concat(phy.edgeLength);
      }
      phy.Nnode = phy.Nnode + 1;
    }
  } else {
    phy.rootEdge = null;
    const INV = new Array(N).fill(false);
    const w = [];
    for (let i = 0; i < N; i++) {
      if (e2[i] === newroot) w.push(i);
    }
    const anc = w.map(idx => e1[idx]);
    let i_idx = w.length > 0 ? w[0] : -1;
    let nod = anc.length > 0 ? anc[0] : null;
    if (nod !== ROOT) {
      w.forEach(idx => {
        INV[idx] = true;
      });
      i_idx = w[0] - 1;
      while (i_idx >= 0) {
        if (e2[i_idx] === nod) {
          if (e1[i_idx] === ROOT) break;
          INV[i_idx] = true;
          nod = e1[i_idx];
        }
        i_idx--;
      }
    }
    if (!fuseRoot && i_idx >= 0) {
      INV[i_idx] = true;
    }
    if (fuseRoot) {
      const k = [];
      for (let i = 0; i < N; i++) {
        if (e1[i] === ROOT) k.push(i);
      }
      let k_idx;
      if (k.length >= 2 && k[1] > w[0]) {
        k_idx = k[1];
      } else {
        k_idx = k[0];
      }
      phy.edge[k_idx][0] = phy.edge[i_idx][1];
      if (wbl) {
        phy.edgeLength[k_idx] =
          phy.edgeLength[k_idx] + phy.edgeLength[i_idx];
      }
    }
    if (fuseRoot) phy.Nnode = oldNnode - 1;
    if (edgelabel && phy.nodeLabel) {
      for (let i = 0; i < N; i++) {
        if (INV[i]) {
          const idx1 = e1[i] - n - 1;
          const idx2 = e2[i] - n - 1;
          if (phy.nodeLabel[idx2] !== undefined) {
            phy.nodeLabel[idx1] = phy.nodeLabel[idx2];
          }
        }
      }
      const newrootIdx = newroot - n - 1;
      if (phy.nodeLabel[newrootIdx] !== undefined) {
        phy.nodeLabel[newrootIdx] = "";
      }
    }
    for (let i = 0; i < N; i++) {
      if (INV[i]) {
        const temp = phy.edge[i][0];
        phy.edge[i][0] = phy.edge[i][1];
        phy.edge[i][1] = temp;
      }
    }
    if (fuseRoot && i_idx >= 0) {
      phy.edge.splice(i_idx, 1);
      if (wbl) phy.edgeLength.splice(i_idx, 1);
    }
    if (resolveRoot) {
      const newnod = oldNnode + n + 1;
      if (outgroup.length === 1) {
        const wh = [];
        for (let i = 0; i < N; i++) {
          if (phy.edge[i][1] === outgroup[0]) wh.push(i);
        }
        const k = [];
        for (let i = 0; i < N; i++) {
          if (phy.edge[i][0] === newroot) k.push(i);
        }
        k.forEach(idx => {
          if (!wh.includes(idx)) {
            phy.edge[idx][0] = newnod;
          }
        });
        const o = [];
        for (let i = 0; i < N; i++) {
          if (!wh.includes(i)) o.push(i);
        }
        o.push(...wh);
        phy.edge = [[newroot, newnod]].concat(o.map(i => phy.edge[i]));
        if (wbl) {
          phy.edgeLength = [0].concat(o.map(i => phy.edgeLength[i]));
        }
      } else {
        const wh = [];
        for (let i = 0; i < N; i++) {
          if (phy.edge[i][0] === newroot) wh.push(i);
        }
        wh.slice(1).forEach(idx => {
          phy.edge[idx][0] = newnod;
        });
        const s1 = phy.edge.slice(0, wh[1]);
        const s2 = phy.edge.slice(wh[1]);
        phy.edge = s1.concat([[newroot, newnod]], s2);
        if (wbl) {
          const s1len = phy.edgeLength.slice(0, wh[1]);
          const s2len = phy.edgeLength.slice(wh[1]);
          phy.edgeLength = s1len.concat([0], s2len);
        }
      }
      phy.Nnode = phy.Nnode + 1;
    }
  }
  const totalNodes = n + phy.Nnode;
  const newNb = new Array(totalNodes).fill(0);
  newNb[newroot - 1] = n + 1;
  const sndcolIndices = [];
  for (let i = 0; i < phy.edge.length; i++) {
    if (phy.edge[i][1] > n) {
      sndcolIndices.push(i);
    }
  }
  const nodesToRenumber = Array.from(
    new Set(sndcolIndices.map(i => phy.edge[i][1]))
  ).sort((a, b) => a - b);
  let num = n + 2;
  nodesToRenumber.forEach(nodeVal => {
    newNb[nodeVal - 1] = num;
    num++;
  });

  // Renumber child column first
  for (const i of sndcolIndices) {
    const oldChild = phy.edge[i][1];
    phy.edge[i][1] = newNb[oldChild - 1];
  }
  // Then renumber parents
  for (let i = 0; i < phy.edge.length; i++) {
    const oldParent = phy.edge[i][0];
    phy.edge[i][0] = newNb[oldParent - 1];
  }

  if (phy.nodeLabel) {
    const newNbSub = newNb.slice(n);
    if (newNbSub.length && fuseRoot) {
      newNbSub.shift();
      phy.nodeLabel.shift();
    }
    const entries = [];
    for (let i = 0; i < newNbSub.length; i++) {
      entries.push({ val: newNbSub[i], idx: i });
    }
    entries.sort((a, b) => a.val - b.val);
    phy.nodeLabel = entries.map(e => phy.nodeLabel[e.idx]);
    if (resolveRoot) {
      for (let i = 0; i < phy.nodeLabel.length; i++) {
        if (phy.nodeLabel[i] == null) {
          phy.nodeLabel[i] = phy.nodeLabel[0];
        }
      }
      phy.nodeLabel[0] = "Root";
    }
  }
  if (phy.order !== undefined) {
    delete phy.order;
  }
  phy = reorder(phy);
  return phy;
}

/**
 * Unroots a phylogenetic tree.
 *
 * @param {object} phy - A phylogenetic tree object in "phylo" format.
 * @param {number} n - The number of tips in the tree.
 * @param {boolean} [collapseSinglesFlag=false] - If true, collapse single-child nodes first.
 * @param {boolean} [keepRootEdge=false] - If true, retains the root edge by adding a terminal edge.
 * @returns {object} The unrooted tree.
 * @throws {Error} If the tree has too few edges or if all nodes are singleton.
 */
export function unrootPhylo(phy, n, collapseSinglesFlag = false, keepRootEdge = false) {
  if (collapseSinglesFlag) {
    if (typeof collapseSingles === "function") {
      phy = collapseSingles(phy);
    } else {
      throw new Error("Utility function 'collapseSingles' is missing.");
    }
  }
  let N = phy.edge.length;
  if (N < 3) {
    throw new Error("cannot unroot a tree with less than three edges.");
  }
  const totalNodes = n + phy.Nnode;
  const dgr = new Array(totalNodes).fill(0);
  for (let i = 0; i < phy.edge.length; i++) {
    const parent = phy.edge[i][0];
    const child = phy.edge[i][1];
    if (parent >= 1 && parent <= totalNodes) dgr[parent - 1]++;
    if (child >= 1 && child <= totalNodes) dgr[child - 1]++;
  }
  if (dgr.every(count => count < 3)) {
    throw new Error("cannot unroot a tree where all nodes are singleton");
  }
  if (phy.rootEdge == null) {
    keepRootEdge = false;
  } else {
    if (!keepRootEdge) phy.rootEdge = null;
  }
  if (!isRooted(phy, n)) return phy;

  const wbl = phy.edgeLength != null;
  let ROOT = n + 1;
  let basal = dgr[ROOT - 1] === 1;
  if (keepRootEdge) basal = false;

  if (keepRootEdge) {
    for (let i = 0; i < phy.edge.length; i++) {
      for (let j = 0; j < 2; j++) {
        if (phy.edge[i][j] > n) {
          phy.edge[i][j] = phy.edge[i][j] + 1;
        }
      }
    }
    ROOT++;
    phy.edge.push([ROOT, n + 1]);
    if (wbl) {
      phy.edgeLength.push(phy.rootEdge);
    }
    phy.rootEdge = null;
    phy.tipLabel.push("[ROOT]");
    n++;
    N++;
  } else {
    if (basal) {
      const i = phy.edge.findIndex(row => row[0] === ROOT);
      if (i === -1) {
        throw new Error("Basal edge not found.");
      }
      const NEWROOT = phy.edge[i][1];
      ROOT++;
      for (let r = 0; r < phy.edge.length; r++) {
        for (let c = 0; c < 2; c++) {
          if (phy.edge[r][c] === NEWROOT) {
            phy.edge[r][c] = ROOT;
          }
        }
      }
      n++;
      phy.edge[i][0] = ROOT;
      phy.edge[i][1] = n;
      let newlab;
      if (phy.nodeLabel && phy.nodeLabel.length > 0) {
        newlab = phy.nodeLabel.shift();
      } else {
        newlab = "[ROOT]";
      }
      phy.tipLabel.push(newlab);
      phy.Nnode--;
    }
  }
  let ophy = phy.order;
  if (!ophy || (ophy !== "cladewise" && ophy !== "postorder")) {
    phy = reorder(phy, "cladewise");
    ophy = phy.order;
  }
  let EDGEROOT = [];
  let NEWROOT;
  if (ophy !== "cladewise") {
    NEWROOT = phy.edge[N - 2][0];
    EDGEROOT = [N - 1, N - 2];
    if (phy.edge[EDGEROOT[0]][1] !== NEWROOT) {
      EDGEROOT.reverse();
    }
  } else {
    const edgerootCandidates = [];
    for (let i = 0; i < phy.edge.length; i++) {
      if (phy.edge[i][0] === ROOT) {
        edgerootCandidates.push(i);
      }
    }
    EDGEROOT = edgerootCandidates.slice();
    if (phy.edge[EDGEROOT[0]][1] <= n) {
      EDGEROOT.reverse();
    }
    NEWROOT = phy.edge[EDGEROOT[0]][1];
  }
  let idxToRemove = EDGEROOT[0];
  let idxOther = EDGEROOT[1];
  if (idxToRemove < idxOther) {
    idxOther--;
  }
  phy.edge.splice(EDGEROOT[0], 1);
  if (wbl) {
    phy.edgeLength[idxOther] =
      phy.edgeLength[idxOther] + phy.edgeLength[EDGEROOT[0]];
    phy.edgeLength.splice(EDGEROOT[0], 1);
  }
  phy.Nnode--;

  for (let i = 0; i < phy.edge.length; i++) {
    for (let j = 0; j < 2; j++) {
      if (phy.edge[i][j] === NEWROOT) {
        phy.edge[i][j] = ROOT;
      }
    }
  }
  for (let i = 0; i < phy.edge.length; i++) {
    for (let j = 0; j < 2; j++) {
      if (phy.edge[i][j] > NEWROOT) {
        phy.edge[i][j] = phy.edge[i][j] - 1;
      }
    }
  }
  if (phy.nodeLabel) {
    if (NEWROOT === n + 2) {
      phy.nodeLabel.shift();
    } else {
      const lbs = phy.nodeLabel.slice();
      const idxToRemoveLabel = NEWROOT - n - 1;
      const tmp = lbs[idxToRemoveLabel];
      lbs.splice(0, 1);
      lbs.splice(idxToRemoveLabel - 1, 1);
      phy.nodeLabel = [tmp].concat(lbs);
    }
  }
  return phy;
}
