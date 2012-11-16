/*
 * (C) Copyright 2012 Dr. Stefan Funke
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */

package de.tourenplaner.chconstruction;

import java.util.Random;

public class RAMGraph extends SGraph {


    // DATA for VERTICES
    float[] xCoord;
    float[] yCoord;
    int[] altID;// EXT; alternative  ID of node if it was derived from another network; for contracted
    // network in particular the original ID
    int[] height;

    int[] OSMID;

    int[] level;            // hold level in contraction hierarchy

    int[] outEdgeOffset;    // outEdges[i] denotes the start of i's outgoing edges in the edgeList
    int[] inEdgeOffset;    // inEdges[i] denotes the start of i's incoming edges in the edgeList

    // EdgeLists (not containing real edge data, but only IDs of real edges in edgeList
    // we maintain two lists
    int[] inEdgeList;        // list of edge IDs ordered according to TARGET
    int[] outEdgeList;        // list of edge IDs ordered according to SOURCE

    // DATA for EDGES
    // potential tweaks: ordering edges according to source/target/whatever
    int[] edgeSource;        // source of i-th edge
    int[] edgeTarget;        // target of i-th edge
    int[] edgeWeight;        // cost of i-th edge
    int[] edgeLength;      // euclidian distance of i-th edge
    int[] edgeAltitudeDifference;      // only positive altitude of i-th edge
    int[] edgeSkippedA;    // ID of the first skipped edge (in the previous network) in case of shortcuts
    int[] edgeSkippedB;    // ID of the second skipped edge

    // variables for online creation of a graph
    int edgesAdded;
    int nodesAdded;

    // random generator
    Random generator = new Random();


    RAMGraph(int nofNodes, int nofEdges) {
        this.nofNodes = nofNodes;
        this.nofEdges = nofEdges;
        nodesAdded = edgesAdded = 0;
        setupMemory();
    }

    RAMGraph(RAMGraph _original) {
        nofNodes = _original.nofNodes;
        nofEdges = _original.nofEdges;
        setupMemory();
        for (int i = 0; i < nofNodes; i++) {
            xCoord[i] = _original.xCoord[i];
            yCoord[i] = _original.yCoord[i];
            altID[i] = _original.altID[i];
            height[i] = _original.height[i];
            OSMID[i] = _original.OSMID[i];
            level[i] = _original.level[i];
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


    void setupMemory()    // sets up memory for a graph with nofNodes and nofEdges
    {
        xCoord = new float[nofNodes];
        yCoord = new float[nofNodes];
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

    RAMGraph compressGraph()    // creates new graph object with exactly added nodes/edges
    {
        RAMGraph resultGraph = new RAMGraph(nodesAdded, edgesAdded);

        for (int i = 0; i < nodesAdded; i++)
            resultGraph.addNode(xCoord(i), yCoord(i), altNodeID(i), height(i), OSMID(i), level[i]);
        for (int j = 0; j < edgesAdded; j++) {
            int curSrc = edgeSource(j), curTrg = edgeTarget(j), curWeight = edgeWeight(j), curLength = edgeLength(j), curAltDiff = edgeAltitudeDifference(j),
                    curA = edgeSkippedA(j), curB = edgeSkippedB(j);
            resultGraph.addEdge(curSrc, curTrg, curWeight, curLength, curAltDiff, curA, curB);
        }
        resultGraph.setupOffsets();

        return resultGraph;
    }


    void setupOffsets()    // computes edgeOffsets for given vertices and edges
    {
        System.err.println("SetupOffsets for " + nofNodes + "/" + nofEdges + " and added " + nodesAdded + "/" + edgesAdded);
        int[] outCount = new int[nofNodes];
        int[] inCount = new int[nofNodes];

        // count how many out/in-edges there are for each node
        for (int j = 0; j < nofEdges; j++) {
            //System.out.println("edge "+j+" target= "+edgeTarget[j]+ ", skippedA="+edgeSkippedA[j]);
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

    RAMGraph pruneGraphSelfloops()
    // get rid of superfluous edges
    {
        int selfLoops = 0;

        boolean[] survivorEdge = new boolean[nofEdges()];
        for (int j = 0; j < nofEdges(); j++)
            survivorEdge[j] = true;
        for (int i = 0; i < nofNodes(); i++)
            for (int j = 0; j < nofOutEdges(i); j++) {
                int curEdge = outEdgeID(i, j);
                int curTarget = edgeTarget(curEdge);
                if (curTarget == i) {
                    survivorEdge[curEdge] = false;
                    selfLoops++;
                }
            }
        int newNofEdges = 0;
        for (int j = 0; j < nofEdges(); j++)
            if (survivorEdge[j])
                newNofEdges++;
        RAMGraph resultGraph = new RAMGraph(nofNodes(), newNofEdges);

        for (int i = 0; i < nofNodes(); i++)
            resultGraph.addNode(xCoord(i), yCoord(i), altNodeID(i), height(i), OSMID(i), level[i]);
        for (int j = 0; j < nofEdges(); j++) {
            if (survivorEdge[j]) {
                int curSrc = edgeSource(j), curTrg = edgeTarget(j), curWeight = edgeWeight(j), curLength = edgeLength(j), curAltDiff = edgeAltitudeDifference(j),
                        curA = edgeSkippedA(j), curB = edgeSkippedB(j);
                resultGraph.addEdge(curSrc, curTrg, curWeight, curLength, curAltDiff, curA, curB);
            }
        }
        resultGraph.setupOffsets();
        System.err.println("pruneGraph: " + nofEdges() + "/" + newNofEdges + " with " + selfLoops + " selfLoops");
        return resultGraph;
    }

    RAMGraph pruneGraph()
    // get rid of superfluous edges
    {
        int selfLoops = 0;

        boolean[] survivorEdge = new boolean[nofEdges()];
        for (int j = 0; j < nofEdges(); j++)
            survivorEdge[j] = true;
        for (int i = 0; i < nofNodes(); i++)
            for (int j = 0; j < nofOutEdges(i); j++) {
                int curEdge = outEdgeID(i, j);
                int curTarget = edgeTarget(curEdge);
                int curWeight = edgeWeight(curEdge);
                if (curTarget == i) {
                    survivorEdge[curEdge] = false;
                    selfLoops++;
                }
                for (int jj = j + 1; jj < nofOutEdges(i); jj++) {
                    int cur2Edge = outEdgeID(i, jj);
                    int cur2Target = edgeTarget(cur2Edge);
                    int cur2Weight = edgeWeight(cur2Edge);
                    if ((cur2Target == curTarget) && (cur2Weight <= curWeight))    // let the last survive!
                        survivorEdge[curEdge] = false;
                }
            }
        int newNofEdges = 0;
        for (int j = 0; j < nofEdges(); j++)
            if (survivorEdge[j])
                newNofEdges++;
        RAMGraph resultGraph = new RAMGraph(nofNodes(), newNofEdges);

        for (int i = 0; i < nofNodes(); i++)
            resultGraph.addNode(xCoord(i), yCoord(i), altNodeID(i), height(i), OSMID(i), level[i]);
        for (int j = 0; j < nofEdges(); j++) {
            if (survivorEdge[j]) {
                int curSrc = edgeSource(j), curTrg = edgeTarget(j), curWeight = edgeWeight(j), curLength = edgeLength(j), curAltDiff = edgeAltitudeDifference(j),
                        curA = edgeSkippedA(j), curB = edgeSkippedB(j);
                resultGraph.addEdge(curSrc, curTrg, curWeight, curLength, curAltDiff, curA, curB);
            }
        }
        resultGraph.setupOffsets();
        System.err.println("pruneGraph: " + nofEdges() + "/" + newNofEdges + " with " + selfLoops + " selfLoops");
        return resultGraph;
    }

    RAMGraph rearrangeGraph() {
        // rearrange graph according to levels of the nodes (small levels first)
        // does not rearrange within the nodes of one level
        RAMGraph resultGraph = new RAMGraph(nodesAdded, edgesAdded);

        System.err.println("We have a graph with " + nodesAdded + " nodes and " + edgesAdded + " edges");
        int[] old2new = new int[nodesAdded];
        int maxLevel = 0;
        for (int i = 0; i < nodesAdded; i++) {
            if ((maxLevel < level[i]) && (level[i] != Integer.MAX_VALUE))
                maxLevel = level[i];
        }
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
            resultGraph.addNodeAt(xCoord[i], yCoord[i], altID[i], height[i], OSMID[i], level[i], old2new[i]);
            levelCount[level[i]]--;
        }
        for (int i = 0; i < nodesAdded; i++)
            assert (levelCount[level[i]] == 0);

        for (int j = 0; j < edgesAdded; j++) {

            int curSrc = edgeSource(j), curTrg = edgeTarget(j), curWeight = edgeWeight(j), curLength = edgeLength(j), curAltDiff = edgeAltitudeDifference(j),
                    curA = edgeSkippedA(j), curB = edgeSkippedB(j);
            resultGraph.addEdge(old2new[curSrc], old2new[curTrg], curWeight, curLength, curAltDiff, curA, curB);

        }
        System.err.println("Before sorting");
        resultGraph.quickSortEdgeArray(0, edgesAdded - 1);
        System.err.println("After sorting");

        resultGraph.setupOffsets();

        System.err.println("Graph rearrangement complete!");
        return resultGraph;
    }

    void sanityCheck() {
        int minWeight = Integer.MAX_VALUE;
        int maxWeight = 0;
        long inSum = 0, inDegSum = 0;
        for (int i = 0; i < nofNodes(); i++) {
            inDegSum += nofInEdges(i);
            for (int j = 0; j < nofInEdges(i); j++) {
                int curEdge = inEdgeID(i, j);
                inSum += edgeWeight(curEdge);
                if (minWeight > edgeWeight(curEdge))
                    minWeight = edgeWeight(curEdge);
                if (maxWeight < edgeWeight(curEdge))
                    maxWeight = edgeWeight(curEdge);
            }
        }
        System.err.println(nofNodes + "/" + nofEdges + ": Sum of inEdges=" + inSum + " with degree sum=" + inDegSum);

        long outSum = 0, outDegSum = 0;
        for (int i = 0; i < nofNodes(); i++) {
            outDegSum += nofOutEdges(i);
            for (int j = 0; j < nofOutEdges(i); j++) {
                int curEdge = outEdgeID(i, j);
                outSum += edgeWeight(curEdge);
                if (minWeight > edgeWeight(curEdge))
                    minWeight = edgeWeight(curEdge);
                if (maxWeight < edgeWeight(curEdge))
                    maxWeight = edgeWeight(curEdge);
            }
        }
        System.err.println(nofNodes + "/" + nofEdges + ": Sum of outEdges=" + outSum + " with degree sum=" + outDegSum);

        System.err.println("MaxWeight: " + maxWeight + " and MinWeight: " + minWeight + " AvgWeight:"
                + (inSum + outSum) / nofEdges());

        // statistics about levels

        int[] levelCount = new int[9999];
        int maxLevel = 0;
        for (int i = 0; i < 9999; i++)
            levelCount[i] = 0;
        for (int i = 0; i < nofNodes(); i++) {
            levelCount[level(i)]++;
            if (level(i) > maxLevel)
                maxLevel = level(i);
        }
        int levelSum = 0;
        for (int i = maxLevel; i >= 0; i--) {
            levelSum += levelCount[i];
            System.err.println(i + ": " + levelSum);
        }
    }


    void setCHShortCuts() {
        //only makes real sense for CH computations
        System.err.println("Setting CH shortcuts");
        int count_shortcuts = 0;
        for (int j = 0; j < nofEdges(); j++) {
            int edgeSrc = edgeSource(j);
            int edgeTrg = edgeTarget(j);
            int edgeCst = edgeWeight(j);

            for (int sj = 0; sj < nofOutEdges(edgeSrc); sj++) {
                int cur_out_edge = outEdgeID(edgeSrc, sj);
                int cur_out_edge_target = edgeTarget(cur_out_edge);
                int cur_out_edge_cost = edgeWeight(cur_out_edge);
                for (int tj = 0; tj < nofInEdges(edgeTrg); tj++) {
                    int cur_in_edge = inEdgeID(edgeTrg, tj);
                    int cur_in_edge_source = edgeSource(cur_in_edge);
                    int cur_in_edge_cost = edgeWeight(cur_in_edge);

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


    void addNode(float _x, float _y, int _altID, int _height, int _OSMID) {
        assert (nodesAdded < nofNodes);
        xCoord[nodesAdded] = _x;
        yCoord[nodesAdded] = _y;
        altID[nodesAdded] = _altID;
        height[nodesAdded] = _height;
        OSMID[nodesAdded] = _OSMID;
        nodesAdded++;
    }

    int addNode(float _x, float _y, int _altID, int _height, int _OSMID, int _level) {
        assert (nodesAdded < nofNodes);
        xCoord[nodesAdded] = _x;
        yCoord[nodesAdded] = _y;
        altID[nodesAdded] = _altID;
        height[nodesAdded] = _height;
        OSMID[nodesAdded] = _OSMID;
        level[nodesAdded] = _level;
        nodesAdded++;
        return (nodesAdded - 1);    // return ID of added node
    }

    void addNodeAt(float _x, float _y, int _altID, int _height, int _OSMID, int _level, int _pos) {    // add node at specified position; voids nodesAdded variable!!!!
        nodesAdded = nofNodes;
        assert (_pos < nofNodes);
        xCoord[_pos] = _x;
        yCoord[_pos] = _y;
        altID[_pos] = _altID;
        height[_pos] = _height;
        OSMID[_pos] = _OSMID;
        level[_pos] = _level;
    }

    void addEdge(int _src, int _trg, int _weight, int _length, int _alt_diff) {
        edgeSource[edgesAdded] = _src;
        edgeTarget[edgesAdded] = _trg;
        edgeWeight[edgesAdded] = _weight;
        edgeLength[edgesAdded] = _length;
        edgeAltitudeDifference[edgesAdded] = _alt_diff;
        edgesAdded++;
    }

    void addEdge(int _src, int _trg, int _weight, int _length, int _alt_diff, int _skipA, int _skipB) {
        edgeSource[edgesAdded] = _src;
        edgeTarget[edgesAdded] = _trg;
        edgeWeight[edgesAdded] = _weight;
        edgeLength[edgesAdded] = _length;
        edgeAltitudeDifference[edgesAdded] = _alt_diff;
        edgeSkippedA[edgesAdded] = _skipA;
        edgeSkippedB[edgesAdded] = _skipB;
        edgesAdded++;
    }


    // overwritten methods
    float xCoord(int nodeID) {
        return xCoord[nodeID];
    }

    float yCoord(int nodeID) {
        return yCoord[nodeID];
    }

    int altNodeID(int nodeID) {
        return altID[nodeID];
    }

    int height(int nodeID) {
        return height[nodeID];
    }

    int OSMID(int nodeID) {
        return OSMID[nodeID];
    }

    int level(int nodeID) {
        return level[nodeID];
    }

    int nofOutEdges(int nodeID) {
        return (outEdgeOffset[nodeID + 1] - outEdgeOffset[nodeID]);
    }

    int nofInEdges(int nodeID) {
        return (inEdgeOffset[nodeID + 1] - inEdgeOffset[nodeID]);
    }

    int outEdgeID(int nodeID, int edgePos)    // returns edge ID of edgePos-th outEdge of nodeID
    {
        return outEdgeList[outEdgeOffset[nodeID] + edgePos];
    }

    int inEdgeID(int nodeID, int edgePos)    // returns edge ID of edgePos-th inEdge of nodeID
    {
        return inEdgeList[inEdgeOffset[nodeID] + edgePos];
    }


    int edgeSource(int edgeID) {
        return edgeSource[edgeID];
    }

    int edgeTarget(int edgeID) {
        return edgeTarget[edgeID];
    }

    int edgeWeight(int edgeID) {
        return edgeWeight[edgeID];
    }

    int edgeLength(int edgeID) {
        return edgeLength[edgeID];
    }

    int edgeAltitudeDifference(int edgeID) {
        return edgeAltitudeDifference[edgeID];
    }

    int edgeSkippedA(int edgeID) {
        return edgeSkippedA[edgeID];
    }

    int edgeSkippedB(int edgeID) {
        return edgeSkippedB[edgeID];
    }
}
