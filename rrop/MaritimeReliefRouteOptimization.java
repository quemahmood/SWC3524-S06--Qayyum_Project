import java.util.*;

public class MaritimeReliefRouteOptimization {

    // Distance Matrix (Adjacency Matrix)
    static int[][] distanceMatrix = {
            { 0, 15, 25, 35 },
            { 15, 0, 30, 28 },
            { 25, 30, 0, 20 },
            { 35, 28, 20, 0 }
    };

    // Location names
    static String[] locations = { "Port A", "Port B", "Relief Center C", "Relief Center D" };

    // --- Global variables for Backtracking and Divide & Conquer to store the
    // result ---
    // These are used because the main methods must return a string,
    // but the recursive helpers find the minimum cost.

    // For Backtracking
    static int globalMinCost = Integer.MAX_VALUE;
    static String globalBestPath = "";

    // For Divide & Conquer
    static int dcMinCost = Integer.MAX_VALUE;
    static String dcBestPath = "";

    /**
     * Solves the TSP using a simple Greedy Algorithm.
     * Starts at Port A, and at each step, moves to the nearest unvisited location.
     */
    public static String greedyTSP(int[][] dist) {
        int n = dist.length;
        boolean[] visited = new boolean[n];
        int current = 0; // Start at Port A
        int totalDistance = 0;
        StringBuilder path = new StringBuilder(locations[0]);

        visited[0] = true;
        int visitedCount = 1;

        while (visitedCount < n) {
            int nearest = -1;
            int minDistance = Integer.MAX_VALUE;

            for (int next = 0; next < n; next++) {
                if (!visited[next] && dist[current][next] < minDistance) {
                    minDistance = dist[current][next];
                    nearest = next;
                }
            }

            if (nearest != -1) {
                visited[nearest] = true;
                visitedCount++;
                totalDistance += minDistance;
                current = nearest;
                path.append(" -> ").append(locations[current]);
            } else {
                // Should not happen in a complete graph
                break;
            }
        }

        // Return to start
        totalDistance += dist[current][0];
        path.append(" -> ").append(locations[0]);

        return "Greedy TSP Route: " + path.toString() + " | Total Distance: " + totalDistance + " nm";
    }

    /**
     * Solves the TSP using Dynamic Programming (Held-Karp algorithm).
     */
    public static String dynamicProgrammingTSP(int[][] dist) {
        int n = dist.length;
        int VISITED_ALL = (1 << n) - 1;
        int[][] memo = new int[n][1 << n];
        String[][] paths = new String[n][1 << n]; // To reconstruct the path

        // Initialize memo table with -1 (unvisited)
        for (int[] row : memo) {
            Arrays.fill(row, -1);
        }

        int totalDistance = dynamicProgrammingTSPHelper(0, 1, dist, memo, VISITED_ALL, paths);

        // Reconstruct the path string from the `paths` array
        String path = locations[0] + " -> " + paths[0][1];
        return "Dynamic Programming TSP Route: " + path + " | Total Distance: " + totalDistance + " nm";
    }

    /**
     * Helper for Dynamic Programming TSP.
     * * @param pos       Current location index
     * @param mask      Bitmask of visited locations
     * @param dist      Distance matrix
     * @param memo      Memoization table for costs
     * @param VISITED_ALL Target bitmask (all 1s)
     * @param paths     Table to store the path taken
     * @return Minimum cost from this state (pos, mask)
     */
    private static int dynamicProgrammingTSPHelper(int pos, int mask, int[][] dist,
            int[][] memo, int VISITED_ALL, String[][] paths) {

        // Base case: If all cities are visited, return the cost to go back to the start
        if (mask == VISITED_ALL) {
            paths[pos][mask] = locations[0]; // Path from last city is to the start
            return dist[pos][0];
        }

        // Memoization: If this state is already computed, return it
        if (memo[pos][mask] != -1) {
            return memo[pos][mask];
        }

        int minCost = Integer.MAX_VALUE;
        int bestNext = -1;

        // Try visiting an unvisited city `next`
        for (int next = 0; next < dist.length; next++) {
            // Check if `next` city is not visited (bit is 0)
            if ((mask & (1 << next)) == 0) {
                // Visit `next` city
                int newMask = mask | (1 << next);
                int cost = dist[pos][next] + dynamicProgrammingTSPHelper(next, newMask, dist, memo, VISITED_ALL, paths);

                if (cost < minCost) {
                    minCost = cost;
                    bestNext = next;
                }
            }
        }

        // Store the path and cost
        if (bestNext != -1) {
            paths[pos][mask] = locations[bestNext] + " -> " + paths[bestNext][mask | (1 << bestNext)];
        }
        memo[pos][mask] = minCost;
        return minCost;
    }

    /**
     * SolVes the TSP using a recursive Backtracking (brute-force) algorithm.
     */
    public static String backtrackingTSP(int[][] dist) {
        int n = dist.length;
        boolean[] visited = new boolean[n];
        // Reset global variables for this run
        globalMinCost = Integer.MAX_VALUE;
        globalBestPath = "";

        // Start at Port A (index 0)
        visited[0] = true;
        StringBuilder path = new StringBuilder(locations[0]);
        tspBacktracking(0, dist, visited, n, 1, 0, path);

        return globalBestPath; // This global variable is set by the helper
    }

    /**
     * Helper for Backtracking TSP.
     * Finds the minimum cost path and updates global variables.
     */
    private static int tspBacktracking(int pos, int[][] dist, boolean[] visited,
            int n, int count, int cost, StringBuilder path) {

        // Base case: If all cities have been visited
        if (count == n) {
            int totalCost = cost + dist[pos][0]; // Add cost to return home
            if (totalCost < globalMinCost) {
                globalMinCost = totalCost;
                // Save the complete path string
                globalBestPath = "Backtracking TSP Route: " + path.toString() + " -> " + locations[0]
                        + " | Total Distance: " + globalMinCost + " nm";
            }
            return dist[pos][0]; // Cost to return home
        }

        int minCostFromHere = Integer.MAX_VALUE;

        // Recursive step: Try all unvisited neighbors
        for (int next = 0; next < n; next++) {
            if (!visited[next]) {
                visited[next] = true;
                path.append(" -> ").append(locations[next]);

                // Recursive call
                int costFromNext = tspBacktracking(next, dist, visited, n, count + 1, cost + dist[pos][next], path);
                minCostFromHere = Math.min(minCostFromHere, dist[pos][next] + costFromNext);

                // Backtrack
                visited[next] = false;
                path.setLength(path.lastIndexOf(" -> "));
            }
        }
        return minCostFromHere;
    }

    /**
     * Solves the TSP using a Divide and Conquer approach (which is functionally
     * similar to backtracking for this problem).
     */
    public static String divideAndConquerTSP(int[][] dist) {
        int n = dist.length;
        boolean[] visited = new boolean[n];
        // Reset global variables for this run
        dcMinCost = Integer.MAX_VALUE;
        dcBestPath = "";

        // Start at Port A (index 0)
        visited[0] = true;
        StringBuilder path = new StringBuilder(locations[0]);
        divideAndConquerHelper(0, visited, 0, dist, n, path);

        return dcBestPath; // This global variable is set by the helper
    }

    /**
     * Helper for Divide and Conquer TSP.
     * Finds the minimum cost path and updates global variables.
     */
    private static int divideAndConquerHelper(int pos, boolean[] visited,
            int currentCost, int[][] dist, int n, StringBuilder path) {

        // Base case: If all cities have been visited
        if (allVisited(visited)) {
            int totalCost = currentCost + dist[pos][0]; // Add cost to return home
            if (totalCost < dcMinCost) {
                dcMinCost = totalCost;
                // Save the complete path string
                dcBestPath = "Divide & Conquer TSP Route: " + path.toString() + " -> " + locations[0]
                        + " | Total Distance: " + dcMinCost + " nm";
            }
            return dist[pos][0]; // Cost to return home
        }

        int minCostFromHere = Integer.MAX_VALUE;

        // Recursive step: Try all unvisited neighbors
        for (int next = 0; next < n; next++) {
            if (!visited[next]) {
                visited[next] = true;
                path.append(" -> ").append(locations[next]);

                // Recursive call
                int costFromNext = divideAndConquerHelper(next, visited, currentCost + dist[pos][next], dist, n, path);
                minCostFromHere = Math.min(minCostFromHere, dist[pos][next] + costFromNext);

                // Backtrack
                visited[next] = false;
                path.setLength(path.lastIndexOf(" -> "));
            }
        }
        return minCostFromHere;
    }

    /**
     * Helper method to check if all locations have been visited.
     */
    private static boolean allVisited(boolean[] visited) {
        for (boolean v : visited) {
            if (!v) {
                return false;
            }
        }
        return true;
    }

    /**
     * Sorts an array in-place using Insertion Sort.
     * Returns a string representation (as required by the stub).
     */
    public static String insertionSort(int[] arr) {
        int n = arr.length;
        for (int i = 1; i < n; ++i) {
            int key = arr[i];
            int j = i - 1;

            // Move elements of arr[0..i-1], that are greater than key,
            // to one position ahead of their current position
            while (j >= 0 && arr[j] > key) {
                arr[j + 1] = arr[j];
                j = j - 1;
            }
            arr[j + 1] = key;
        }
        // The main method expects this to sort in-place, but the stub requires a
        // String return.
        // We return the string, but the main method relies on the in-place sort.
        return Arrays.toString(arr);
    }

    /**
     * Searches for a target in a sorted array using Binary Search.
     * Returns the index as a String, or "-1" if not found.
     */
    public static String binarySearch(int[] arr, int target) {
        int low = 0;
        int high = arr.length - 1;
        while (low <= high) {
            int mid = low + (high - low) / 2;
            if (arr[mid] == target) {
                return String.valueOf(mid); // Found
            } else if (arr[mid] < target) {
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
        return "-1"; // Not found
    }

    /**
     * Min-Heap implementation using PriorityQueue.
     * This class is defined as required by the main method.
     */
    static class MinHeap {
        // TO BE IMPLEMENTED: Heap Algorithm logic
        // Use Java's built-in PriorityQueue as a min-heap
        private PriorityQueue<Integer> heap = new PriorityQueue<>();

        /**
         * Inserts a value into the min-heap.
         */
        public void insert(int value) {
            heap.add(value);
        }

        /**
         * Removes and returns the minimum value from the heap.
         * Returns -1 or throws exception if heap is empty (here, we'll let it
         * return null, which poll() does, and main will get a NullPointerException
         * if empty. Or, we can check.)
         */
        public int extractMin() {
            if (heap.isEmpty()) {
                // Or throw new NoSuchElementException("Heap is empty");
                return -1;
            }
            return heap.poll();
        }
    }

    /**
     * Splay Tree implementation.
     */
    static class SplayTree {
        private Node root;

        // Node class for the Splay Tree
        private static class Node {
            int key;
            Node left, right;

            Node(int key) {
                this.key = key;
                this.left = null;
                this.right = null;
            }
        }

        // --- Splay Tree Core Operations ---

        // Right rotate
        private Node rightRotate(Node x) {
            Node y = x.left;
            x.left = y.right;
            y.right = x;
            return y;
        }

        // Left rotate
        private Node leftRotate(Node x) {
            Node y = x.right;
            x.right = y.left;
            y.left = x;
            return y;
        }

        /**
         * Splay operation: brings the node with the given key to the root.
         * If key is not in the tree, brings the last accessed node to the root.
         */
        private Node splay(Node root, int key) {
            if (root == null || root.key == key) {
                return root;
            }

            if (key < root.key) {
                // Key is in the left subtree
                if (root.left == null)
                    return root;

                if (key < root.left.key) {
                    // Zig-Zig (Left-Left)
                    root.left.left = splay(root.left.left, key);
                    root = rightRotate(root);
                } else if (key > root.left.key) {
                    // Zig-Zag (Left-Right)
                    root.left.right = splay(root.left.right, key);
                    if (root.left.right != null) {
                        root.left = leftRotate(root.left);
                    }
                }

                return (root.left == null) ? root : rightRotate(root);

            } else {
                // Key is in the right subtree
                if (root.right == null)
                    return root;

                if (key < root.right.key) {
                    // Zag-Zig (Right-Left)
                    root.right.left = splay(root.right.left, key);
                    if (root.right.left != null) {
                        root.right = rightRotate(root.right);
                    }
                } else if (key > root.right.key) {
                    // Zag-Zag (Right-Right)
                    root.right.right = splay(root.right.right, key);
                    root = leftRotate(root);
                }

                return (root.right == null) ? root : leftRotate(root);
            }
        }

        /**
         * Inserts a key into the Splay Tree.
         */
        public void insert(int key) {
            if (root == null) {
                root = new Node(key);
                return;
            }

            // Bring the closest leaf or the key itself to the root
            root = splay(root, key);

            // If key is already present
            if (root.key == key)
                return;

            Node newNode = new Node(key);

            if (key < root.key) {
                newNode.right = root;
                newNode.left = root.left;
                root.left = null;
            } else {
                newNode.left = root;
                newNode.right = root.right;
                root.right = null;
            }
            root = newNode;
        }

        /**
         * Searches for a key in the Splay Tree.
         * Spleys the node to the root.
         * * @return true if found, false otherwise.
         */
        public boolean search(int key) {
            root = splay(root, key);
            return (root != null && root.key == key);
        }

    } // end of SplayTree class

    // Driver method (as provided in the prompt)
    public static void main(String[] args) {
        System.out.println(greedyTSP(distanceMatrix));
        System.out.println(dynamicProgrammingTSP(distanceMatrix));
        System.out.println(backtrackingTSP(distanceMatrix));
        System.out.println(divideAndConquerTSP(distanceMatrix));

        // Sorting and Searching
        int[] arr = { 8, 3, 5, 1, 9, 2 };
        insertionSort(arr); // Sorts in-place
        System.out.println("Sorted Array: " + Arrays.toString(arr));
        System.out.println("Binary Search (5 found at index): " + binarySearch(arr, 5));

        // Min-Heap Test
        MinHeap heap = new MinHeap();
        heap.insert(10);
        heap.insert(3);
        heap.insert(15);
        System.out.println("Min-Heap Extract Min: " + heap.extractMin());

        // Splay Tree Test
        SplayTree tree = new SplayTree();
        tree.insert(20);
        tree.insert(10);
        tree.insert(30);
        System.out.println("Splay Tree Search (10 found): " + tree.search(10));
    }
}