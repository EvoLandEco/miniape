/**
 * Tabulates the counts of each integer in an array.
 *
 * @param {Array<number>} arr - An array of positive integers.
 * @returns {Array<number>} An array of counts where index i gives the frequency of i in arr.
 */
export function tabulate(arr) {
    const maxVal = Math.max(...arr);
    const counts = new Array(maxVal + 1).fill(0);
    for (const val of arr) {
        counts[val] = (counts[val] || 0) + 1;
    }
    return counts;
}

/**
* Computes the rank of each number in an array, mimicking R's base rank function.
*
* @param {Array<number>} arr - Array of numbers (with missing values as NaN or null).
* @param {Object} [options={}] - Options object.
* @param {boolean} [options.naLast=true] - If true, missing values are ranked last; if false, ranked first.
* @param {string} [options.tiesMethod="average"] - Method for handling ties:
*         "average", "first", "last", "random", "max", or "min".
* @returns {Array<number>} Array of ranks.
*/
export function rank(arr, { naLast = true, tiesMethod = "average" } = {}) {
    const n = arr.length;
    const nonMissing = [];
    const missingIndices = [];

    // Separate non-missing values and record indices for missing values.
    for (let i = 0; i < n; i++) {
        const x = arr[i];
        if (x === null || Number.isNaN(x)) {
            missingIndices.push(i);
        } else {
            nonMissing.push({ val: x, idx: i });
        }
    }

    // Sort nonMissing by value; if equal, by original index.
    nonMissing.sort((a, b) => {
        if (a.val === b.val) return a.idx - b.idx;
        return a.val - b.val;
    });

    const ranks = new Array(n);
    let currentRank = 1;
    let i = 0;

    // Process non-missing values groupwise (ties)
    while (i < nonMissing.length) {
        let j = i;
        while (j < nonMissing.length && nonMissing[j].val === nonMissing[i].val) {
            j++;
        }
        const groupSize = j - i;

        // Compute ranks for the tie group based on tiesMethod.
        if (tiesMethod === "first") {
            for (let k = i; k < j; k++) {
                ranks[nonMissing[k].idx] = currentRank++;
            }
        } else if (tiesMethod === "last") {
            for (let k = i; k < j; k++) {
                // Reverse order: first element gets the highest rank in group.
                ranks[nonMissing[k].idx] = currentRank + groupSize - 1 - (k - i);
            }
            currentRank += groupSize;
        } else if (tiesMethod === "min") {
            for (let k = i; k < j; k++) {
                ranks[nonMissing[k].idx] = currentRank;
            }
            currentRank += groupSize;
        } else if (tiesMethod === "max") {
            for (let k = i; k < j; k++) {
                ranks[nonMissing[k].idx] = currentRank + groupSize - 1;
            }
            currentRank += groupSize;
        } else if (tiesMethod === "average") {
            const avg = currentRank + (groupSize - 1) / 2;
            for (let k = i; k < j; k++) {
                ranks[nonMissing[k].idx] = avg;
            }
            currentRank += groupSize;
        } else if (tiesMethod === "random") {
            // Create an array of possible ranks for the group.
            const groupRanks = [];
            for (let k = 0; k < groupSize; k++) {
                groupRanks.push(currentRank + k);
            }
            // Shuffle the groupRanks array.
            for (let k = groupRanks.length - 1; k > 0; k--) {
                const randIdx = Math.floor(Math.random() * (k + 1));
                [groupRanks[k], groupRanks[randIdx]] = [groupRanks[randIdx], groupRanks[k]];
            }
            for (let k = i; k < j; k++) {
                ranks[nonMissing[k].idx] = groupRanks[k - i];
            }
            currentRank += groupSize;
        } else {
            throw new Error("Unknown tiesMethod: " + tiesMethod);
        }
        i = j;
    }

    // Handle missing values.
    if (naLast) {
        // Missing values get ranks after all non-missing values.
        for (let i = 0; i < missingIndices.length; i++) {
            ranks[missingIndices[i]] = currentRank++;
        }
    } else {
        // Missing values get ranks at the beginning.
        const numMissing = missingIndices.length;
        for (let i = 0; i < missingIndices.length; i++) {
            ranks[missingIndices[i]] = i + 1;
        }
        // Shift non-missing ranks by the number of missing values.
        for (let i = 0; i < n; i++) {
            if (ranks[i] !== undefined) {
                ranks[i] += numMissing;
            }
        }
    }

    return ranks;
}
