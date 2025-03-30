/**
 * @file phylo.js
 * @module phylo
 * 
 * This module defines the Phylo class, a JavaScript representation of a phylogenetic tree,
 * mimicking the "phylo" object from R's ape package.
 * 
 * Components:
 *  - edge: A 2D array representing tree edges, where each sub-array is [parent, child].
 *          **Note:** The tree is assumed to be 0-indexed: tip nodes are numbered 0 to (nTips-1),
 *          and internal nodes are numbered nTips to (nTips + Nnode - 1).
 *  - edgeLength: An array of branch lengths corresponding to each edge.
 *  - tipLabel: An array of tip labels (leaf names).
 *  - nodeLabel: An array of node labels (internal node names).
 *  - Nnode: The number of internal nodes.
 *  - rootEdge: (Optional) The branch length of the root edge.
 */

class Phylo {
    /**
     * Creates an instance of Phylo.
     *
     * @param {Array<Array<number>>} edge - A 2D array where each sub-array represents an edge as [parent, child].
     *                                      **Assumption:** The tree is 0-indexed: tip nodes are numbered 0 to (nTips-1),
     *                                      and internal nodes are numbered nTips to (nTips+Nnode-1).
     * @param {Array<number>} edgeLength - An array of branch lengths corresponding to each edge.
     * @param {Array<string>} tipLabel - An array of tip labels (leaf names).
     * @param {Array<string>} nodeLabel - An array of node labels (internal node names).
     * @param {number} Nnode - The number of internal nodes.
     * @param {number|null} [rootEdge=null] - Optional branch length from the root.
     */
    constructor(edge, edgeLength, tipLabel, nodeLabel, Nnode, rootEdge = null) {
      this._edge = edge;
      this._edgeLength = edgeLength;
      this._tipLabel = tipLabel;
      this._nodeLabel = nodeLabel;
      this._Nnode = Nnode;
      this._rootEdge = rootEdge;
    }
  
    /**
     * Gets the edge matrix.
     *
     * @return {Array<Array<number>>} The edge matrix.
     */
    get edge() {
      return this._edge;
    }
  
    /**
     * Sets the edge matrix.
     *
     * @param {Array<Array<number>>} newEdge - The new edge matrix.
     */
    set edge(newEdge) {
      this._edge = newEdge;
    }
  
    /**
     * Gets the branch lengths.
     *
     * @return {Array<number>} The branch lengths.
     */
    get edgeLength() {
      return this._edgeLength;
    }
  
    /**
     * Sets the branch lengths.
     *
     * @param {Array<number>} newEdgeLength - The new branch lengths.
     */
    set edgeLength(newEdgeLength) {
      this._edgeLength = newEdgeLength;
    }
  
    /**
     * Gets the tip labels.
     *
     * @return {Array<string>} The tip labels.
     */
    get tipLabel() {
      return this._tipLabel;
    }
  
    /**
     * Sets the tip labels.
     *
     * @param {Array<string>} newTipLabel - The new tip labels.
     */
    set tipLabel(newTipLabel) {
      this._tipLabel = newTipLabel;
    }
  
    /**
     * Gets the node labels.
     *
     * @return {Array<string>} The node labels.
     */
    get nodeLabel() {
      return this._nodeLabel;
    }
  
    /**
     * Sets the node labels.
     *
     * @param {Array<string>} newNodeLabel - The new node labels.
     */
    set nodeLabel(newNodeLabel) {
      this._nodeLabel = newNodeLabel;
    }
  
    /**
     * Gets the number of internal nodes.
     *
     * @return {number} The number of internal nodes.
     */
    get Nnode() {
      return this._Nnode;
    }
  
    /**
     * Sets the number of internal nodes.
     *
     * @param {number} newNnode - The new number of internal nodes.
     */
    set Nnode(newNnode) {
      this._Nnode = newNnode;
    }
  
    /**
     * Gets the branch length from the root.
     *
     * @return {number|null} The root edge length or null if not set.
     */
    get rootEdge() {
      return this._rootEdge;
    }
  
    /**
     * Sets the branch length from the root.
     *
     * @param {number|null} newRootEdge - The new root edge length.
     */
    set rootEdge(newRootEdge) {
      this._rootEdge = newRootEdge;
    }
  
    // TODO: Implement more methods here
  }
  
  export default Phylo;
  