/*
 * (C) Copyright 2012 Dr. Stefan Funke
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */

package de.tourenplaner.chconstruction.graph;

import fmi.graph.metaio.MetaData;

import java.util.Random;

public class RAMGraph extends SGraph {

    private MetaData meta;

    private int nofNodes;
    private int nofEdges;

    // DATA for VERTICES
    private float[] lats;
    private float[] lons;
    private int[] altID;// EXT; alternative  ID of node if it was derived from another network; for contracted
    // network in particular the original ID
    private int[] height;

    private int[] OSMID;

    private int[] level;            // hold level in contraction hierarchy

    private int[] outEdgeOffset;    // outEdges[i] denotes the start of i's outgoing edges in the edgeList
    private int[] inEdgeOffset;    // inEdges[i] denotes the start of i's incoming edges in the edgeList

    // EdgeLists (not containing real edge data, but only IDs of real edges in edgeList
    // we maintain two lists
    private int[] inEdgeList;        // list of edge IDs ordered according to TARGET
    private int[] outEdgeList;        // list of edge IDs ordered according to SOURCE

    // DATA for EDGES
    // potential tweaks: ordering edges according to source/target/whatever
    private int[] edgeSource;        // source of i-th edge
    private int[] edgeTarget;        // target of i-th edge
    private int[] edgeWeight;        // cost of i-th edge
    private int[] edgeLength;      // euclidian distance of i-th edge
    private int[] edgeAltitudeDifference;      // only positive altitude of i-th edge
    private int[] edgeSkippedA;    // ID of the first skipped edge (in the previous network) in case of shortcuts
    private int[] edgeSkippedB;    // ID of the second skipped edge

    // variables for online creation of a graph
    private int edgesAdded;
    private int nodesAdded;

    // usefull data about graph
    private int maxLevel;

    // random generator
    Random generator = new Random();

    /**
     * Create a new graph with nofNodes nodes and nofEdges edges
     * as well as associated meta data meta (can be null if not available)
     * @param nofNodes
     * @param nofEdges
     * @param meta
     */
    public RAMGraph(int nofNodes, int nofEdges, MetaData meta) {
        this.meta = meta;
        this.nofNodes = nofNodes;
        this.nofEdges = nofEdges;
        nodesAdded = edgesAdded = 0;
        maxLevel = -1;
        setupMemory();
        generator.setSeed(42);
    }

    /**
     * Create a partial copy of orignal graph, use .clone() for
     * real copy. Note that we can only do partial copies "from the start"
     * up to newNofNodes, newNofEdges
     * @param _original
     * @param newNofEdges
     * @param newNofNodes
     */
    public RAMGraph(RAMGraph _original, int newNofNodes, int newNofEdges) {
        meta = _original.meta;
        nofNodes = newNofNodes;
        nofEdges = newNofEdges;
        maxLevel = 0;
        setupMemory();
        for (int i = 0; i < nofNodes; i++) {
            lats[i] = _original.lats[i];
            lons[i] = _original.lons[i];
            altID[i] = _original.altID[i];
            height[i] = _original.height[i];
            OSMID[i] = _original.OSMID[i];
            level[i] = _original.level[i];
            if (level[i] > maxLevel)
                maxLevel = level[i];
            outEdgeOffset[i] = _original.outEdgeOffset[i];
            inEdgeOffset[i] = _original.inEdgeOffset[i];
        }
        outEdgeOffset[nofNodes] = inEdgeOffset[nofNodes] = nofEdges;
        for (int i = 0; i < nofEdges; i++) {
            inEdgeList[i] = _original.inEdgeList[i];
            outEdgeList[i] = _original.outEdgeList[i];
            edgeSource[i] = _original.edgeSource[i];
            edgeTarget[i] = _original.edgeTarget[i];
            edgeWeight[i] = _original.edgeWeight[i];
            edgeLength[i] = _original.edgeLength[i];
            edgeAltitudeDifference[i] = _original.edgeAltitudeDifference[i];
            edgeSkippedA[i] = _original.edgeSkippedA[i];
            edgeSkippedB[i] = _original.edgeSkippedB[i];
        }

    }

    @Override
    public  MetaData getMetaData(){
        return meta;
    }


    void setupMemory()    // sets up memory for a graph with nofNodes and nofEdges
    {
        lats = new float[nofNodes];
        lons = new float[nofNodes];
        altID = new int[nofNodes];
        height = new int[nofNodes];
        OSMID = new int[nofNodes];
        level = new int[nofNodes];
        outEdgeOffset = new int[nofNodes + 1];
        inEdgeOffset = new int[nofNodes + 1];

        inEdgeList = new int[nofEdges + 1];
        outEdgeList = new int[nofEdges + 1];

        edgeSource = new int[nofEdges];
        edgeTarget = new int[nofEdges];
        edgeWeight = new int[nofEdges];
        edgeLength = new int[nofEdges];
        edgeAltitudeDifference = new int[nofEdges];
        edgeSkippedA = new int[nofEdges];
        edgeSkippedB = new int[nofEdges];
    }

    public RAMGraph compressGraph()    // creates new graph object with exactly added nodes/edges
    {
        RAMGraph resultGraph = new RAMGraph(nodesAdded, edgesAdded, this.meta);

        for (int i = 0; i < nodesAdded; i++)
            resultGraph.addNode(getLat(i), getLon(i), getAltNodeID(i), getHeight(i), getOSMID(i), level[i]);
        for (int j = 0; j < edgesAdded; j++) {
            int curSrc = getSource(j), curTrg = getTarget(j), curWeight = getWeight(j), curLength = getEuclidianLength(j), curAltDiff = getAltitudeDifference(j),
                    curA = getSkippedA(j), curB = getSkippedB(j);
            resultGraph.addEdge(curSrc, curTrg, curWeight, curLength, curAltDiff, curA, curB);
        }
        resultGraph.setupOffsets();

        return resultGraph;
    }


    public void setupOffsets()    // computes edgeOffsets for given vertices and edges
    {
        System.err.println("SetupOffsets for " + nofNodes + "/" + nofEdges + " and added " + nodesAdded + "/" + edgesAdded);
        int[] outCount = new int[nofNodes];
        int[] inCount = new int[nofNodes];

        // count how many out/in-edges there are for each node
        for (int j = 0; j < nofEdges; j++) {
            //System.out.println("edge "+j+" target= "+getTarget[j]+ ", skippedA="+getSkippedA[j]);
            outCount[edgeSource[j]]++;
            inCount[edgeTarget[j]]++;
        }

        // set the offsets for the edge lists based on edge count
        int outEdgeSum = 0;
        int inEdgeSum = 0;
        for (int i = 0; i < nofNodes; i++) {
            outEdgeOffset[i] = outEdgeSum;
            inEdgeOffset[i] = inEdgeSum;
            outEdgeSum += outCount[i];
            inEdgeSum += inCount[i];
        }
        outEdgeOffset[nofNodes] = inEdgeOffset[nofNodes] = nofEdges;
        assert (inEdgeSum == outEdgeSum);
        assert (inEdgeSum == nofEdges);

        // now put edgeIDs into out/inEdgeLists
        for (int j = 0; j < nofEdges; j++) {
            int curSrc = edgeSource[j];
            int curTrg = edgeTarget[j];
            int curInOffset = inEdgeOffset[curTrg + 1] - inCount[curTrg];
            int curOutOffset = outEdgeOffset[curSrc + 1] - outCount[curSrc];
            inEdgeList[curInOffset] = j;
            outEdgeList[curOutOffset] = j;
            inCount[curTrg]--;
            outCount[curSrc]--;
        }
        for (int i = 0; i < nofNodes; i++) {
            assert (inCount[i] == 0);
            assert (outCount[i] == 0);
        }

    }

    void swapEdges(int e, int f) {// swaps edges in EdgeArray at positions e and f
        int src, trg, weight, length, altDiff, skipA, skipB;
        src = edgeSource[e];
        trg = edgeTarget[e];
        weight = edgeWeight[e];
        length = edgeLength[e];
        altDiff = edgeAltitudeDifference[e];
        skipA = edgeSkippedA[e];
        skipB = edgeSkippedB[e];

        edgeSource[e] = edgeSource[f];
        edgeTarget[e] = edgeTarget[f];
        edgeWeight[e] = edgeWeight[f];
        edgeLength[e] = edgeLength[f];
        edgeAltitudeDifference[e] = edgeAltitudeDifference[f];
        edgeSkippedA[e] = edgeSkippedA[f];
        edgeSkippedB[e] = edgeSkippedB[f];

        edgeSource[f] = src;
        edgeTarget[f] = trg;
        edgeWeight[f] = weight;
        edgeLength[f] = length;
        edgeAltitudeDifference[f] = altDiff;
        edgeSkippedA[f] = skipA;
        edgeSkippedB[f] = skipB;
    }

    void quickSortEdgeArray(int start, int end) {
        // sorts Edge Array according to <source, target> including boundaries (!)
        if (start >= end)
            return;
        //int pivot=(start+end)/2;	// of course suboptimal pivot :-)
        int pivot = start + generator.nextInt(end - start);
        swapEdges(pivot, end);        // move pivot to the end
        int storage = start;            // all elements up to [storage-1] are to be
        // smaller than pivot
        for (int j = start; j < end; j++) {
            if ((edgeSource[j] < edgeSource[end])
                    || ((edgeSource[j] == edgeSource[end])
                    && (edgeTarget[j] < edgeTarget[end]))) {
                swapEdges(j, storage);
                storage++;
            }
        }
        swapEdges(storage, end);    // now have pivot at storage
        //if (end-start>1024)
        //	System.err.println(start+"--"+end+" with pivot "+storage+" and "+pivot);
        quickSortEdgeArray(start, storage - 1);
        quickSortEdgeArray(storage + 1, end);
    }

    public RAMGraph pruneGraphSelfloops()
    // get rid of superfluous edges
    {
        int selfLoops = 0;

        boolean[] survivorEdge = new boolean[nofEdges()];
        for (int j = 0; j < nofEdges(); j++)
            survivorEdge[j] = true;
        for (int i = 0; i < nofNodes(); i++)
            for (int j = 0; j < nofOutEdges(i); j++) {
                int curEdge = outEdgeID(i, j);
                int curTarget = getTarget(curEdge);
                if (curTarget == i) {
                    survivorEdge[curEdge] = false;
                    selfLoops++;
                }
            }
        int newNofEdges = 0;
        for (int j = 0; j < nofEdges(); j++)
            if (survivorEdge[j])
                newNofEdges++;
        RAMGraph resultGraph = new RAMGraph(nofNodes(), newNofEdges, this.meta);

        for (int i = 0; i < nofNodes(); i++)
            resultGraph.addNode(getLat(i), getLon(i), getAltNodeID(i), getHeight(i), getOSMID(i), level[i]);
        for (int j = 0; j < nofEdges(); j++) {
            if (survivorEdge[j]) {
                int curSrc = getSource(j), curTrg = getTarget(j), curWeight = getWeight(j), curLength = getEuclidianLength(j), curAltDiff = getAltitudeDifference(j),
                        curA = getSkippedA(j), curB = getSkippedB(j);
                resultGraph.addEdge(curSrc, curTrg, curWeight, curLength, curAltDiff, curA, curB);
            }
        }
        resultGraph.setupOffsets();
        System.err.println("pruneGraph: " + nofEdges() + "/" + newNofEdges + " with " + selfLoops + " selfLoops");
        return resultGraph;
    }

    public RAMGraph pruneGraph()
    // get rid of superfluous edges
    {
        int selfLoops = 0;

        boolean[] survivorEdge = new boolean[nofEdges()];
        for (int j = 0; j < nofEdges(); j++)
            survivorEdge[j] = true;
        for (int i = 0; i < nofNodes(); i++)
            for (int j = 0; j < nofOutEdges(i); j++) {
                int curEdge = outEdgeID(i, j);
                int curTarget = getTarget(curEdge);
                int curWeight = getWeight(curEdge);
                if (curTarget == i) {
                    survivorEdge[curEdge] = false;
                    selfLoops++;
                }
                for (int jj = j + 1; jj < nofOutEdges(i); jj++) {
                    int cur2Edge = outEdgeID(i, jj);
                    int cur2Target = getTarget(cur2Edge);
                    int cur2Weight = getWeight(cur2Edge);
                    if ((cur2Target == curTarget) && (cur2Weight <= curWeight))    // let the last survive!
                        survivorEdge[curEdge] = false;
                }
            }
        int newNofEdges = 0;
        for (int j = 0; j < nofEdges(); j++)
            if (survivorEdge[j])
                newNofEdges++;
        RAMGraph resultGraph = new RAMGraph(nofNodes(), newNofEdges, this.meta);

        for (int i = 0; i < nofNodes(); i++)
            resultGraph.addNode(getLat(i), getLon(i), getAltNodeID(i), getHeight(i), getOSMID(i), level[i]);
        for (int j = 0; j < nofEdges(); j++) {
            if (survivorEdge[j]) {
                int curSrc = getSource(j), curTrg = getTarget(j), curWeight = getWeight(j), curLength = getEuclidianLength(j), curAltDiff = getAltitudeDifference(j),
                        curA = getSkippedA(j), curB = getSkippedB(j);
                resultGraph.addEdge(curSrc, curTrg, curWeight, curLength, curAltDiff, curA, curB);
            }
        }
        resultGraph.setupOffsets();
        System.err.println("pruneGraph: " + nofEdges() + "/" + newNofEdges + " with " + selfLoops + " selfLoops");
        return resultGraph;
    }

    private void ensureMaxLevel(){
        if (maxLevel == -1) {
            for (int i = 0; i < nodesAdded; i++) {
                if ((maxLevel < level[i]) && (level[i] != Integer.MAX_VALUE)) maxLevel = level[i];
            }
        }
    }

    public RAMGraph rearrangeGraph() {
        // rearrange graph according to levels of the nodes (small levels first)
        // does not rearrange within the nodes of one getLevel
        RAMGraph resultGraph = new RAMGraph(nodesAdded, edgesAdded, this.meta);

        System.err.println("We have a graph with " + nodesAdded + " nodes and " + edgesAdded + " edges");
        int[] old2new = new int[nodesAdded];
        ensureMaxLevel();
        int countIMAX = 0;
        for (int i = 0; i < nodesAdded; i++) {
            if (level[i] == Integer.MAX_VALUE) {
                level[i] = maxLevel + 1;
                countIMAX++;
            }
        }
        System.err.println("We had " + countIMAX + " high levels");
        maxLevel++;

        System.err.println("We have maxlevel=" + maxLevel);
        int[] levelCount = new int[maxLevel + 2];
        int[] levelOffset = new int[maxLevel + 2];

        for (int l = 0; l <= maxLevel + 1; l++)
            levelCount[l] = 0;

        for (int i = 0; i < nodesAdded; i++)
            levelCount[level[i]]++;

        levelOffset[0] = 0;
        for (int l = 1; l <= maxLevel + 1; l++)
            levelOffset[l] = levelOffset[l - 1] + levelCount[l - 1];
        assert (levelOffset[maxLevel + 1] == nodesAdded);

        for (int i = 0; i < nodesAdded; i++) {// iterate through all levels and put nodes in respective order
            assert (levelCount[level[i]] > 0);
            old2new[i] = levelOffset[level[i] + 1] - levelCount[level[i]];
            resultGraph.addNodeAt(lats[i], lons[i], altID[i], height[i], OSMID[i], level[i], old2new[i]);
            levelCount[level[i]]--;
        }
        for (int i = 0; i < nodesAdded; i++)
            assert (levelCount[level[i]] == 0);

        for (int j = 0; j < edgesAdded; j++) {

            int curSrc = getSource(j), curTrg = getTarget(j), curWeight = getWeight(j), curLength = getEuclidianLength(j), curAltDiff = getAltitudeDifference(j),
                    curA = getSkippedA(j), curB = getSkippedB(j);
            resultGraph.addEdge(old2new[curSrc], old2new[curTrg], curWeight, curLength, curAltDiff, curA, curB);

        }
        System.err.println("Before sorting");
        resultGraph.quickSortEdgeArray(0, edgesAdded - 1);
        System.err.println("After sorting");

        resultGraph.setupOffsets();

        System.err.println("Graph rearrangement complete!");
        return resultGraph;
    }

    /**
     * Makes sure that edgeWeight[i] > 0 for all edges
     * by setting all lower weights to 1.
     * Safe weights are assumed when setting up CH Shortcuts!
     */
    public void safeWeights(){
        for (int i = 0; i < nofEdges; i++){
            if (edgeWeight[i] <= 0)
                edgeWeight[i] = 1;
        }
    }

    public void sanityCheck() {
        int minWeight = Integer.MAX_VALUE;
        int maxWeight = 0;
        long inSum = 0, inDegSum = 0;
        for (int i = 0; i < nofNodes(); i++) {
            inDegSum += nofInEdges(i);
            for (int j = 0; j < nofInEdges(i); j++) {
                int curEdge = inEdgeID(i, j);
                inSum += getWeight(curEdge);
                if (minWeight > getWeight(curEdge))
                    minWeight = getWeight(curEdge);
                if (maxWeight < getWeight(curEdge))
                    maxWeight = getWeight(curEdge);
            }
        }
        System.err.println(nofNodes + "/" + nofEdges + ": Sum of inEdges=" + inSum + " with degree sum=" + inDegSum);

        long outSum = 0, outDegSum = 0;
        for (int i = 0; i < nofNodes(); i++) {
            outDegSum += nofOutEdges(i);
            for (int j = 0; j < nofOutEdges(i); j++) {
                int curEdge = outEdgeID(i, j);
                outSum += getWeight(curEdge);
                if (minWeight > getWeight(curEdge))
                    minWeight = getWeight(curEdge);
                if (maxWeight < getWeight(curEdge))
                    maxWeight = getWeight(curEdge);
            }
        }
        System.err.println(nofNodes + "/" + nofEdges + ": Sum of outEdges=" + outSum + " with degree sum=" + outDegSum);

        System.err.println("MaxWeight: " + maxWeight + " and MinWeight: " + minWeight + " AvgWeight:"
                + (inSum + outSum) / nofEdges());
    }

    public void printLevelStats(){
        ensureMaxLevel();

        int[] levelNodeHist = new int[maxLevel+1];

        for (int i = 0; i < nofNodes(); i++) {
            levelNodeHist[getLevel(i)]++;
        }
        int sumHigher = nofNodes();
        System.err.println("Level Stats Full:");
        System.err.println("Level, # nodes, # higher nodes (including), percent of total");
        for (int i = 0; i <= maxLevel; i++){
            System.err.println(i + ", " + levelNodeHist[i] + ", "+ sumHigher+", "+((double)levelNodeHist[i])/((double)nofNodes())*100.0);
            sumHigher -= levelNodeHist[i];
        }

        int low = 0;
        sumHigher = nofNodes();
        int high, sumNodes;
        System.err.println("Level Stats Summary:");
        System.err.println("Levels, # nodes, # higher nodes (including), percent of total");
        while (low <= maxLevel) {
            sumNodes = 0;
            high = Math.min(low + 10, maxLevel + 1);
            for (int i = low; i < high; i++){
               sumNodes += levelNodeHist[i];
            }
            System.err.println("["+low + "-" + high + "), " + sumNodes + ", " + sumHigher + ", " + ((double) sumNodes) / ((double) nofNodes()) * 100.0);
            sumHigher -= sumNodes;
            low = high;
        }

    }


    public void setCHShortCuts() {
        //only makes real sense for CH computations
        System.err.println("Setting CH shortcuts");
        int count_shortcuts = 0;
        for (int j = 0; j < nofEdges(); j++) {
            int edgeSrc = getSource(j);
            int edgeTrg = getTarget(j);
            int edgeCst = getWeight(j);

            for (int sj = 0; sj < nofOutEdges(edgeSrc); sj++) {
                int cur_out_edge = outEdgeID(edgeSrc, sj);
                int cur_out_edge_target = getTarget(cur_out_edge);
                int cur_out_edge_cost = getWeight(cur_out_edge);
                for (int tj = 0; tj < nofInEdges(edgeTrg); tj++) {
                    int cur_in_edge = inEdgeID(edgeTrg, tj);
                    int cur_in_edge_source = getSource(cur_in_edge);
                    int cur_in_edge_cost = getWeight(cur_in_edge);

                    if ((cur_out_edge_target == cur_in_edge_source)
                            && (edgeCst == cur_out_edge_cost + cur_in_edge_cost)
                            && (cur_out_edge_cost != 0)
                            && (cur_in_edge_cost != 0)) {
                        edgeSkippedA[j] = cur_out_edge;
                        edgeSkippedB[j] = cur_in_edge;
                        count_shortcuts++;
                    }
                }
            }
        }
        System.err.println("We have found " + count_shortcuts + " shortcuts");
    }


    public void addNode(float _lat, float _lon, int _altID, int _height, int _OSMID) {
        assert (nodesAdded < nofNodes);
        lats[nodesAdded] = _lat;
        lons[nodesAdded] = _lon;
        altID[nodesAdded] = _altID;
        height[nodesAdded] = _height;
        OSMID[nodesAdded] = _OSMID;
        nodesAdded++;
    }

    public int addNode(float _lat, float _lon, int _altID, int _height, int _OSMID, int _level) {
        assert (nodesAdded < nofNodes);
        lats[nodesAdded] = _lat;
        lons[nodesAdded] = _lon;
        altID[nodesAdded] = _altID;
        height[nodesAdded] = _height;
        OSMID[nodesAdded] = _OSMID;
        level[nodesAdded] = _level;
        nodesAdded++;
        return (nodesAdded - 1);    // return ID of added node
    }

    public void addNodeAt(float _lat, float _lon, int _altID, int _height, int _OSMID, int _level, int _pos) {    // add node at specified position; voids nodesAdded variable!!!!
        nodesAdded = nofNodes;
        assert (_pos < nofNodes);
        lats[_pos] = _lat;
        lons[_pos] = _lon;
        altID[_pos] = _altID;
        height[_pos] = _height;
        OSMID[_pos] = _OSMID;
        level[_pos] = _level;
    }

    public void addEdge(int _src, int _trg, int _weight, int _length, int _alt_diff) {
        edgeSource[edgesAdded] = _src;
        edgeTarget[edgesAdded] = _trg;
        edgeWeight[edgesAdded] = _weight;
        edgeLength[edgesAdded] = _length;
        edgeAltitudeDifference[edgesAdded] = _alt_diff;
        edgesAdded++;
    }

    public void addEdge(int _src, int _trg, int _weight, int _length, int _alt_diff, int _skipA, int _skipB) {
        edgeSource[edgesAdded] = _src;
        edgeTarget[edgesAdded] = _trg;
        edgeWeight[edgesAdded] = _weight;
        edgeLength[edgesAdded] = _length;
        edgeAltitudeDifference[edgesAdded] = _alt_diff;
        edgeSkippedA[edgesAdded] = _skipA;
        edgeSkippedB[edgesAdded] = _skipB;
        edgesAdded++;
    }


    @Override
    public int nofEdges() {
        return nofEdges;
    }

    @Override
    public int nofNodes() {
        return nofNodes;
    }

    // overwritten methods
    public float getLat(int nodeID) {
        return lats[nodeID];
    }

    public float getLon(int nodeID) {
        return lons[nodeID];
    }

    public int getAltNodeID(int nodeID) {
        return altID[nodeID];
    }

    public int getHeight(int nodeID) {
        return height[nodeID];
    }

    public int getOSMID(int nodeID) {
        return OSMID[nodeID];
    }

    public int getLevel(int nodeID) {
        return level[nodeID];
    }

    public void setLevel(int nodeID, int newLevel) {
        level[nodeID] = newLevel;
    }

    public int nofOutEdges(int nodeID) {
        return (outEdgeOffset[nodeID + 1] - outEdgeOffset[nodeID]);
    }

    public int nofInEdges(int nodeID) {
        return (inEdgeOffset[nodeID + 1] - inEdgeOffset[nodeID]);
    }

    public int outEdgeID(int nodeID, int edgePos)    // returns edge ID of edgePos-th outEdge of nodeID
    {
        return outEdgeList[outEdgeOffset[nodeID] + edgePos];
    }

    public int inEdgeID(int nodeID, int edgePos)    // returns edge ID of edgePos-th inEdge of nodeID
    {
        return inEdgeList[inEdgeOffset[nodeID] + edgePos];
    }

    public int getSource(int edgeID) {
        return edgeSource[edgeID];
    }

    public int getTarget(int edgeID) {
        return edgeTarget[edgeID];
    }

    public int getWeight(int edgeID) {
        return edgeWeight[edgeID];
    }

    public int getEuclidianLength(int edgeID) {
        return edgeLength[edgeID];
    }

    public int getAltitudeDifference(int edgeID) {
        return edgeAltitudeDifference[edgeID];
    }

    public int getSkippedA(int edgeID) {
        return edgeSkippedA[edgeID];
    }

    public int getSkippedB(int edgeID) {
        return edgeSkippedB[edgeID];
    }
}
