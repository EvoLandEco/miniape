import Phylo from "./phylo.js";

/**
 * Parses a Newick string into a nested tree structure.
 *
 * Each node is an object with optional properties:
 *   - name:   node label.
 *   - length: branch length from its parent.
 *   - children: an array of child nodes.
 *
 * @param {string} newick - A preprocessed Newick string.
 * @returns {Object} The nested tree.
 */
export function parseNested(newick) {
  // Tokenize the string by splitting on parentheses, commas, colons, and semicolons.
  const tokens = newick.split(/([\(\),:;])/).filter(s => s !== "");
  let index = 0;

  function parseSubtree() {
    const node = { children: [] };
    if (tokens[index] === "(") {
      index++; // skip '('
      // Parse first child
      node.children.push(parseSubtree());
      // Parse additional children, if any.
      while (tokens[index] === ",") {
        index++; // skip comma
        node.children.push(parseSubtree());
      }
      if (tokens[index] !== ")") {
        throw new Error("Missing closing parenthesis in Newick string.");
      }
      index++; // skip ')'
      // A label may follow.
      if (tokens[index] && tokens[index] !== ":" && tokens[index] !== "," && tokens[index] !== ")") {
        node.name = tokens[index];
        index++;
      }
    } else {
      // Leaf node: the token is the label.
      node.name = tokens[index];
      index++;
    }
    // A branch length may follow a colon.
    if (tokens[index] === ":") {
      index++; // skip ':'
      node.length = parseFloat(tokens[index]);
      index++;
    }
    return node;
  }

  return parseSubtree();
}

/**
 * Converts a nested tree structure into a phylo object.
 *
 * Numbering scheme:
 *   - Tips (leaves) are numbered from 1 to nTips.
 *   - The root is always numbered nTips + 1.
 *   - All other internal nodes (non-root) are numbered consecutively starting at nTips + 2.
 *
 * This function builds the edge matrix, branch lengths, tip labels, and node (internal) labels.
 *
 * @param {Object} tree - The nested tree structure.
 * @returns {Phylo} A phylo object.
 */
export function convertToPhylo(tree) {
  const edges = [];
  const edgeLengths = [];
  const tipLabels = [];
  const nodeLabels = [];  // For internal nodes; nodeLabels[0] will be the root's label.
  let rootEdge = null;

  // First pass: count the number of tips.
  function countTips(node) {
    if (node.children && node.children.length > 0) {
      return node.children.reduce((sum, child) => sum + countTips(child), 0);
    } else {
      return 1;
    }
  }
  const totalTips = countTips(tree);

  let tipIndex = 1;
  // For non-root internal nodes, numbering will start at totalTips + 2 because the root gets totalTips + 1.
  let nextInternal = totalTips + 2;

  // Recursive function to assign node numbers.
  // The "isRoot" flag distinguishes the root from other internal nodes.
  function assignNumbers(node, isRoot) {
    if (node.children && node.children.length > 0) {
      if (isRoot) {
        // Note: always assign the root number as totalTips + 1.
        node.number = totalTips + 1;
        nodeLabels.push(node.name || "");  // Store root label as first internal node label.
        // Process children (non-root internal nodes will be numbered later).
        node.children.forEach(child => assignNumbers(child, false));
      } else {
        // For non-root internal nodes, assign numbers in postorder.
        node.children.forEach(child => assignNumbers(child, false));
        node.number = nextInternal++;
        nodeLabels.push(node.name || "");
      }
    } else {
      // Leaf node: assign tip numbers.
      node.number = tipIndex++;
      tipLabels.push(node.name || "");
    }
  }
  assignNumbers(tree, true);

  // Build the edge matrix.
  // For every node (except the root) record an edge from its parent to itself.
  function buildEdges(node, parentNum) {
    if (parentNum !== null) {
      edges.push([parentNum, node.number]);
      edgeLengths.push(node.length != null ? node.length : undefined);
    } else {
      // For the root, if a branch length exists, store it as rootEdge.
      if (node.length != null) {
        rootEdge = node.length;
      }
    }
    if (node.children && node.children.length > 0) {
      node.children.forEach(child => buildEdges(child, node.number));
    }
  }
  buildEdges(tree, null);

  // The number of internal nodes is the total count of nodes with children.
  const Nnode = nodeLabels.length;
  const finalEdgeLengths = edgeLengths.some(el => el !== undefined) ? edgeLengths : null;
  return new Phylo(edges, finalEdgeLengths, tipLabels, nodeLabels, Nnode, rootEdge);
}

/**
 * Parses a Newick string and returns a phylo object.
 *
 * Preprocessing steps include:
 *   - Removing comments enclosed in square brackets.
 *   - Handling single-quoted labels (by temporarily replacing them with placeholders).
 *   - Stripping extraneous whitespace and underscores.
 *   - Removing a trailing semicolon.
 *
 * @param {string} newickStr - A string in Newick format.
 * @returns {Phylo} A phylo object representing the parsed tree.
 */
export function parseNewick(newickStr) {
  if (typeof newickStr !== "string") {
    throw new Error("Newick input must be a string.");
  }

  // Remove comments enclosed in square brackets.
  newickStr = newickStr.replace(/\[[^\]]*\]/g, '');

  // Process single quotes: replace quoted labels with placeholders.
  // Store mapping from placeholder to original label.
  const singleQuoteMapping = {};
  newickStr = newickStr.replace(/'([^']*)'/g, (match, p1) => {
    const placeholder = "IMPROBABLEPREFIX" + Object.keys(singleQuoteMapping).length + "IMPROBABLESUFFIX";
    singleQuoteMapping[placeholder] = p1;
    return placeholder;
  });

  // Remove possible leading and trailing underscores.
  newickStr = newickStr.replace(/^_+|_+$/g, '');
  newickStr = newickStr.trim();

  // Remove trailing semicolon, if present.
  if (newickStr.endsWith(";")) {
    newickStr = newickStr.slice(0, -1);
  }

  // Remove spaces and TABs.
  newickStr = newickStr.replace(/[ \t]/g, '');

  // Parse the nested tree structure.
  const nestedTree = parseNested(newickStr);

  // Replace placeholders in node labels with the original quoted text.
  function replacePlaceholders(node) {
    if (node.name) {
      const regex = /IMPROBABLEPREFIX(\d+)IMPROBABLESUFFIX/;
      const match = node.name.match(regex);
      if (match) {
        const placeholder = "IMPROBABLEPREFIX" + match[1] + "IMPROBABLESUFFIX";
        if (singleQuoteMapping.hasOwnProperty(placeholder)) {
          node.name = singleQuoteMapping[placeholder];
        }
      }
    }
    if (node.children && node.children.length > 0) {
      node.children.forEach(child => replacePlaceholders(child));
    }
  }
  replacePlaceholders(nestedTree);

  return convertToPhylo(nestedTree);
}
