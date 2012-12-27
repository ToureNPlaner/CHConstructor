/*
 * (C) Copyright 2012 Dr. Stefan Funke
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */

package de.tourenplaner.chconstruction;

import java.util.PriorityQueue;

/* computes contraction hierarchy for the given graph
 * KEEPS the given order and the altID
 * IDs of the original graph
 * ONLY WORKS with RAMGraph (!)
 * 
 * tempGraph stores current residual graph
 * myCHGraph stores graph with currently added Shortcuts 
 * 
 * levels are set such that myCHGraph allows CH-SP-Queries at any time (after setting up offsets)
 * 
 * skipped edges are set after the very end only !!!
 */

public class CHConstructor extends Constructor {

    RAMGraph myGraph;        // original graph

    RAMGraph tempGraph; // store the residual graph

    CHConstructor(RAMGraph _myGraph) {
        myGraph = _myGraph;
        myCHGraph = new RAMGraph(myGraph.nofNodes(), 3 * myGraph.nofEdges());  /// KONSTANTE GRUSEL!!!!!!!
        tempGraph = new RAMGraph(myGraph.nofNodes(), myGraph.nofEdges());
        for (int i = 0; i < myGraph.nofNodes(); i++)    // first add the original graph
        {
            myCHGraph.addNode(myGraph.getLat(i), _myGraph.getLon(i), myGraph.getAltNodeID(i), myGraph.getHeight(i), myGraph.getOSMID(i), Integer.MAX_VALUE);
            tempGraph.addNode(myGraph.getLat(i), _myGraph.getLon(i), i, myGraph.getHeight(i), myGraph.getOSMID(i), Integer.MAX_VALUE); // here alt denotes ID
        }
        for (int j = 0; j < myGraph.nofEdges(); j++) {
            int orgSrc = myGraph.getSource(j), orgTrg = myGraph.getTarget(j), orgWeight = myGraph.getWeight(j), orgLength = myGraph.getEuclidianLength(j), orgHeight = myGraph.getAltitudeDifference(j);
            myCHGraph.addEdge(orgSrc, orgTrg, orgWeight, orgLength, orgHeight, -1, -1);
            tempGraph.addEdge(orgSrc, orgTrg, orgWeight, orgLength, orgHeight, -1, -1);
        }
        tempGraph.setupOffsets();

        // tempGraph.sanityCheck();
    }

    int contractNode(int curNode, Dijkstra myDijkstra, int[] srcSC, int[] trgSC, int[] wgtSC, int[] lgthSC, int[] hgtSC, int boundSC)    // return number of computed SCs
    {                                                                                                            // stops if boundSC is reached
        int nofSC = 0;
        for (int i = 0; i < tempGraph.nofInEdges(curNode); i++) {
            int curSrcEdge = tempGraph.inEdgeID(curNode, i);
            int curSrc = tempGraph.getSource(curSrcEdge);
            for (int j = 0; j < tempGraph.nofOutEdges(curNode); j++) {
                int curTrgEdge = tempGraph.outEdgeID(curNode, j);
                int curTrg = tempGraph.getTarget(curTrgEdge);
                int weightSC = tempGraph.getWeight(curSrcEdge) + tempGraph.getWeight(curTrgEdge);
                int lengthSC = tempGraph.getEuclidianLength(curSrcEdge) + tempGraph.getEuclidianLength(curTrgEdge);
                int heightSC = tempGraph.getAltitudeDifference(curSrcEdge) + tempGraph.getAltitudeDifference(curTrgEdge);
                myDijkstra.runDijkstra(curSrc, curTrg);
                //if (d==weightSC) // better: check if pred[curTrg]==curNode and pred[curNode]==curSrc
                if (((myDijkstra.pred(curTrg) == curNode) && (myDijkstra.pred(curNode) == curSrc))) {
                    srcSC[nofSC] = curSrc;
                    trgSC[nofSC] = curTrg;
                    wgtSC[nofSC] = weightSC;
                    lgthSC[nofSC] = lengthSC;
                    hgtSC[nofSC] = heightSC;
                    nofSC++;
                }
                if (nofSC == boundSC) return boundSC;
            }
        }
        return nofSC;
    }

    int contractLevel(int newLevel)    // contracts an independent set of the current tempGraph
    {
        PriorityQueue<PQElement> degreePQ = new PriorityQueue<PQElement>();
        boolean[] stillIndep = new boolean[tempGraph.nofNodes()];
        boolean[] contracted = new boolean[tempGraph.nofNodes()];
        int[] candidates = new int[tempGraph.nofNodes()];
        int[] candSCoffset = new int[tempGraph.nofNodes() + 1];        // offsets into list of shortcuts
        int nofCandidates = 0;

        // create priority queue with nodes acc. to their degrees (MULTIPLIED)
        for (int i = 0; i < tempGraph.nofNodes(); i++) {
            stillIndep[i] = true;
            contracted[i] = false;
            int degree = tempGraph.nofInEdges(i) * tempGraph.nofOutEdges(i);
            degreePQ.add(new PQElement(degree, i));
        }

        PQElement curEl;

        // now pick independent set as candidates; greedy, starting with small degrees
        int degSum = 0;
        while (!degreePQ.isEmpty()) {
            curEl = degreePQ.remove();
            int curNode = curEl.value;
            if (stillIndep[curNode]) {
                degSum += curEl.key;
                candidates[nofCandidates] = curNode;
                nofCandidates++;
                for (int j = 0; j < tempGraph.nofInEdges(curNode); j++) {
                    int edgeID = tempGraph.inEdgeID(curNode, j);
                    int src = tempGraph.getSource(edgeID);
                    stillIndep[src] = false;
                }
                for (int j = 0; j < tempGraph.nofOutEdges(curNode); j++) {
                    int edgeID = tempGraph.outEdgeID(curNode, j);
                    int trg = tempGraph.getTarget(edgeID);
                    stillIndep[trg] = false;
                }
            }
        }
        System.err.println("We have an IS of size " + nofCandidates);
        //PQElement curEl=degreePQ.peek();
        int boundSC = degSum / nofCandidates + 1;
        // we know that we can find a node which produces at most boundSC shortcuts !!!
        if (boundSC < 6)
            boundSC = 6;

        System.err.println("boundSC=" + boundSC);


        int allSCUB = tempGraph.nofEdges();
        if (allSCUB < boundSC * nofCandidates)
            allSCUB = boundSC * nofCandidates;

        // allocate memory for all SCs
        int[] srcSCall = new int[allSCUB];
        int[] trgSCall = new int[allSCUB];
        int[] wgtSCall = new int[allSCUB];
        int[] lgthSCall = new int[allSCUB];
        int[] altDiffSCall = new int[allSCUB];


        // instantiate Dijkstra
        Dijkstra myDijkstra = new Dijkstra(tempGraph);

        int tentNofSC;
        candSCoffset[0] = 0;

        // now we try to contract the nodes of the independent set
        PriorityQueue<PQElement> contractionPQ;

        int sumED;
        int validED;

        boundSC = (boundSC + 4) / 2;
        if (newLevel < 3)
            boundSC = 2;
        int candBound = nofCandidates;    // how many candidates to look at
        int validBound = 5 * candBound / 6;        // how many of the considered candidates to fully evaluate (at least)

        do {
            contractionPQ = new PriorityQueue<PQElement>();
            boundSC = boundSC * 2;
            System.err.print("\n Current boundSC: " + boundSC + "  ");
            validED = 0;
            sumED = 0;
            tentNofSC = 0;

            // temporary memory for single contraction
            int[] srcSC = new int[boundSC];
            int[] trgSC = new int[boundSC];
            int[] wgtSC = new int[boundSC];
            int[] lgthSC = new int[boundSC];
            int[] altDiffSC = new int[boundSC];


            for (int i = 0; i < candBound; i++) {
                int nofSC = contractNode(candidates[i], myDijkstra, srcSC, trgSC, wgtSC, lgthSC, altDiffSC, boundSC);
                int edgeDiff = nofSC - tempGraph.nofInEdges(candidates[i]) - tempGraph.nofOutEdges(candidates[i]);
                if (nofSC < boundSC) {
                    sumED += edgeDiff;
                    validED++;
                }

                for (int j = 0; j < nofSC; j++) {
                    srcSCall[tentNofSC] = srcSC[j];
                    trgSCall[tentNofSC] = trgSC[j];
                    wgtSCall[tentNofSC] = wgtSC[j];
                    lgthSCall[tentNofSC] = lgthSC[j];
                    altDiffSCall[tentNofSC] =  altDiffSC[j];
                    tentNofSC++;
                }
                contractionPQ.add(new PQElement(edgeDiff, i));
                candSCoffset[i + 1] = tentNofSC;

                if ((i % (nofCandidates / 10 + 1)) == 0) {
                    System.err.print((10 * i / (nofCandidates / 10 + 1) + "% "));
                    System.err.print("(" + nofSC + "/" + edgeDiff + ") ");

                }
            }
        } while (validED < validBound);

        int newNofNodes = 0;
        int newNofEdges;
        int realContract = 0;
        int totalNofSC = 0;

        // allocate memory for final SCs
        int[] srcSCfinal = new int[allSCUB];
        int[] trgSCfinal = new int[allSCUB];
        int[] wgtSCfinal = new int[allSCUB];
        int[] lgthSCfinal = new int[allSCUB];
        int[] altDiffSCfinal = new int[allSCUB];

        int avgED = sumED / validED + 1;
        System.err.println("\n AvgED=" + avgED + " validED=" + validED);


        while (!contractionPQ.isEmpty()) {
            PQElement curCand = contractionPQ.remove();
            int curNode = curCand.value;
            int curED = curCand.key;
            int curNofSC = candSCoffset[curNode + 1] - candSCoffset[curNode];

            if ((curNofSC < boundSC) &&                                // we contract if ED<=0 but at least 3/4 of valid candidates
                    ((curED <= 0) || (curED < avgED) || (realContract <= validBound)
                    )
                    ) {
                realContract++;
                contracted[candidates[curNode]] = true;
                // now copy its SCs in final SC list
                for (int i = candSCoffset[curNode]; i < candSCoffset[curNode + 1]; i++) {
                    srcSCfinal[totalNofSC] = srcSCall[i];
                    trgSCfinal[totalNofSC] = trgSCall[i];
                    wgtSCfinal[totalNofSC] = wgtSCall[i];
                    lgthSCfinal[totalNofSC] = lgthSCall[i];
                    altDiffSCfinal[totalNofSC] = altDiffSCall[i];
                    totalNofSC++;
                }
            }
        }


        System.err.println("\n Will contract " + realContract + " nodes and creating " + totalNofSC + " SCs");

        // count surviving nodes and edges
        int nofDeletedEdges = 0;
        newNofEdges = totalNofSC;
        for (int i = 0; i < tempGraph.nofNodes(); i++)
            if (!contracted[i])
                newNofNodes++;
        assert (realContract == tempGraph.nofNodes() - newNofNodes);
        for (int j = 0; j < tempGraph.nofEdges(); j++)
            if ((!contracted[tempGraph.getSource(j)]) && (!contracted[tempGraph.getTarget(j)]))
                newNofEdges++;
            else
                nofDeletedEdges++;
        assert ((newNofEdges + nofDeletedEdges) == tempGraph.nofEdges() + totalNofSC);
        System.err.println("\n New Graph has " + newNofNodes + " nodes and " + newNofEdges + " edges having deleted " + nofDeletedEdges);

        // * assign all contracted nodes the new level in myCHgraph
        // * add all created shortcuts to myCHgraph
        // * construct new tempGraph consisting of all surviving nodes and edges and shortcuts


        RAMGraph newTempGraph = new RAMGraph(newNofNodes, newNofEdges);
        int[] old2new = new int[tempGraph.nofNodes()];

        for (int i = 0; i < tempGraph.nofNodes(); i++)
            if (!contracted[i]) {
                old2new[i] = newTempGraph.addNode(tempGraph.getLat(i), tempGraph.getLon(i), tempGraph.getAltNodeID(i), tempGraph.getHeight(i), tempGraph.getOSMID(i), 0);
            } else {
                old2new[i] = -1;
                assert (myCHGraph.getLevel(tempGraph.getAltNodeID(i)) == Integer.MAX_VALUE);
                myCHGraph.setLevel(tempGraph.getAltNodeID(i), newLevel);
            }

        // copy surviving edges to newTempGraph
        for (int j = 0; j < tempGraph.nofEdges(); j++) {
            int curSrc = tempGraph.getSource(j), curTrg = tempGraph.getTarget(j), curWgt = tempGraph.getWeight(j), curEdgeLength = tempGraph.getEuclidianLength(j), curEdgeHeight = tempGraph.getAltitudeDifference(j);
            if ((!contracted[curSrc]) && (!contracted[curTrg]))    // copy edge to newTempGraph
            {
                newTempGraph.addEdge(old2new[curSrc], old2new[curTrg], curWgt, curEdgeLength, curEdgeHeight, -2, -2);
            }
        }

        // now add SC edges to newTempGRaph as well as myCHGraph

        for (int j = 0; j < totalNofSC; j++) {
            newTempGraph.addEdge(old2new[srcSCfinal[j]], old2new[trgSCfinal[j]], wgtSCfinal[j], lgthSCfinal[j], altDiffSCfinal[j], -2, -2);
            myCHGraph.addEdge(tempGraph.getAltNodeID(srcSCfinal[j]), tempGraph.getAltNodeID(trgSCfinal[j]), wgtSCfinal[j], lgthSCfinal[j], altDiffSCfinal[j], -2, -2);
        }

        newTempGraph.setupOffsets();

        tempGraph = newTempGraph.pruneGraph();
        // now create new tempGraph which consists of old graph (with edges between uncontracted nodes)
        // and all shortcuts created
        //
        return tempGraph.nofNodes();

    }


}

