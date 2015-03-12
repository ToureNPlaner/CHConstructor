/*
 * (C) Copyright 2012 Peter Vollmer
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */

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
import java.util.concurrent.*;

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

public class CSPCHConstructor extends Constructor{

    RAMGraph myGraph;        // original graph

    RAMGraph tempGraph; // store the residual graph

    CSPCHConstructor(RAMGraph _myGraph) {
        myGraph = _myGraph;
        myCHGraph = new RAMGraph(myGraph.nofNodes(), 3 * myGraph.nofEdges());  /// KONSTANTE GRUSEL!!!!!!!
        tempGraph = new RAMGraph(myGraph.nofNodes(), myGraph.nofEdges());
        for (int i = 0; i < myGraph.nofNodes(); i++)    // first add the original graph
        {
            myCHGraph.addNode(myGraph.xCoord(i), _myGraph.yCoord(i), myGraph.altNodeID(i), myGraph.height(i), myGraph.OSMID(i), Integer.MAX_VALUE);
            tempGraph.addNode(myGraph.xCoord(i), _myGraph.yCoord(i), i, myGraph.height(i), myGraph.OSMID(i), Integer.MAX_VALUE); // here alt denotes ID
        }
        for (int j = 0; j < myGraph.nofEdges(); j++) {
            int orgSrc = myGraph.edgeSource(j), orgTrg = myGraph.edgeTarget(j), orgWeight = myGraph.edgeWeight(j), orgLength = myGraph.edgeLength(j), orgHeight = myGraph.edgeAltitudeDifference(j);
            myCHGraph.addEdge(orgSrc, orgTrg, orgWeight, orgLength, orgHeight, -1, -1);
            tempGraph.addEdge(orgSrc, orgTrg, orgWeight, orgLength, orgHeight, -1, -1);
        }
        tempGraph.setupOffsets();

        // tempGraph.sanityCheck();
    }

    private class NodeContracter implements Runnable{

        CSPDijkstra cspDijkstra;
        SGraph graph;
        public int curNode;
        public int[] srcSC;
        public int[] trgSC;
        public int[] wgtSC;
        public int[] lgthSC;
        public int[] altSC;
        public int boundSC;
        Semaphore sem;
        int nofSC;

        NodeContracter (SGraph graph, Semaphore sem)
        {
            this.graph = graph;
            this.cspDijkstra = new CSPDijkstra(graph);
            this.sem = sem;
        }

        boolean addShortcut (CSPDijkstra myDijkstra, int curSrc, int curTrg, int weightSC, int altitudeSC) {
            int maxLambda = 2048;
            int upperBound = maxLambda;
            int lowerBound = 0;
            int midLambda;
            while (true) {

                if (upperBound-lowerBound <= 1){
                    return true;
                }
                midLambda = (upperBound + lowerBound)/2;
                int dijkDist=myDijkstra.runDijkstra(curSrc,curTrg,midLambda,maxLambda);
                backTrack(myDijkstra,curSrc,curTrg);
                int distSC=(weightSC-altitudeSC)*midLambda + maxLambda*altitudeSC;
                assert(dijkDist<=distSC);

                if (weightSC == weightSum && altitudeSC == resourceSum){
                    //System.err.println("gleich");
                    return true;
                }
                if(weightSC >= weightSum && altitudeSC >= resourceSum){
                    //System.err.println("dominanter pfad");
                    assert(dijkDist<distSC);
                    return false;
                }
                //lambda of the intersection point of the both functions.
                int lambdaIntersect = maxLambda*(resourceSum-altitudeSC)/((weightSC-weightSum)-(altitudeSC-resourceSum));
                if(lambdaIntersect < lowerBound || lambdaIntersect > upperBound){
                   return false;
                }
                if (weightSC-altitudeSC < weightSum-resourceSum){
                    lowerBound = midLambda;
                    //System.err.println("linke mitte");
                } else {
                    upperBound = midLambda;
                    //System.err.println("rechte mitte");
                }

            }
        }

        int weightSum;
        int resourceSum;
        private void backTrack(CSPDijkstra myDijkstra, int curSrc, int curTrg) {
            weightSum = 0;
            resourceSum = 0;
            int curNode = curTrg;
            int curEdge;
            while (curNode != curSrc){
                curEdge = myDijkstra.pred(curNode);
                weightSum += myDijkstra.myGraph.edgeWeight(curEdge);
                resourceSum += myDijkstra.myGraph.edgeAltitudeDifference(curEdge);
                curNode = myDijkstra.myGraph.edgeSource(curEdge);
            }
        }

        @Override
        public void run() {
            nofSC = 0;
            OUTER: for (int i = 0; i < graph.nofInEdges(curNode); i++) {
                int curSrcEdge = graph.inEdgeID(curNode, i);
                int curSrc = graph.edgeSource(curSrcEdge);
                for (int j = 0; j < graph.nofOutEdges(curNode); j++) {
                    int curTrgEdge = graph.outEdgeID(curNode, j);
                    int curTrg = graph.edgeTarget(curTrgEdge);
                    int weightSC = graph.edgeWeight(curSrcEdge) + graph.edgeWeight(curTrgEdge);
                    int lengthSC = graph.edgeLength(curSrcEdge) + graph.edgeLength(curTrgEdge);
                    int altitudeSC = graph.edgeAltitudeDifference(curSrcEdge) + graph.edgeAltitudeDifference(curTrgEdge);
                    boolean adden = addShortcut(cspDijkstra,curSrc,curTrg,weightSC,altitudeSC);
                    //if (d==weightSC) // better: check if pred[curTrg]==curNode and pred[curNode]==curSrc
                    if (adden){ //&& myDijkstra.pred(curTrg) == curNode) {
                        srcSC[nofSC] = curSrc;
                        trgSC[nofSC] = curTrg;
                        wgtSC[nofSC] = weightSC;
                        lgthSC[nofSC] = lengthSC;
                        altSC[nofSC] = altitudeSC;
                        nofSC++;
                    }
                    if (nofSC == boundSC){
                        break OUTER;
                    }
                }
            }

            //sem.release();
            sem.release(1);
            return;
        }

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
                    int src = tempGraph.edgeSource(edgeID);
                    stillIndep[src] = false;
                }
                for (int j = 0; j < tempGraph.nofOutEdges(curNode); j++) {
                    int edgeID = tempGraph.outEdgeID(curNode, j);
                    int trg = tempGraph.edgeTarget(edgeID);
                    stillIndep[trg] = false;
                }
            }
        }
        System.err.println("We have an IS of size " + nofCandidates);
        for(int i=0; i<nofCandidates; i++)
        {
            assert(stillIndep[candidates[i]]==true);
        }
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
        int nofTask = 28;
        int nofThreads = 14;
        ExecutorService executorService = Executors.newFixedThreadPool(nofThreads);
        NodeContracter[] nodeContracters = new NodeContracter[nofTask];
        Semaphore sem = new Semaphore(0);
        //CSPDijkstra myDijkstra = new CSPDijkstra(tempGraph);
        for (int i = 0 ; i < nodeContracters.length ; ++i){
            nodeContracters[i] = new NodeContracter(tempGraph,sem);
        }

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


            // temporary for each thread
            for (int i = 0 ; i < nodeContracters.length ; ++i){
                nodeContracters[i].srcSC = new int[boundSC];
                nodeContracters[i].trgSC = new int[boundSC];
                nodeContracters[i].wgtSC = new int[boundSC];
                nodeContracters[i].lgthSC = new int[boundSC];
                nodeContracters[i].altSC = new int[boundSC];
                nodeContracters[i].boundSC = boundSC;

            }


            for (int i = 0; i < candBound; i+=nofTask) {
                int end = Math.min(nofTask,candBound-i);
                for (int j = 0 ; j < end; ++j){
                    nodeContracters[j].curNode = candidates[i+j];
                    executorService.submit(nodeContracters[j]);
                }

                try {

                    sem.acquire(end);

                } catch (InterruptedException e) {
                    e.printStackTrace();
                }


                for (int j = 0 ; j < end ; ++j){
                    NodeContracter curNodeContracter = nodeContracters[j];

                    int edgeDiff = curNodeContracter.nofSC - tempGraph.nofInEdges(candidates[i+j]) - tempGraph.nofOutEdges(candidates[i+j]);
                    if (curNodeContracter.nofSC < boundSC) {
                        sumED += edgeDiff;
                        validED++;
                    }

                    for (int k = 0; k < curNodeContracter.nofSC; k++) {
                        srcSCall[tentNofSC] = curNodeContracter.srcSC[k];
                        trgSCall[tentNofSC] = curNodeContracter.trgSC[k];
                        wgtSCall[tentNofSC] = curNodeContracter.wgtSC[k];
                        lgthSCall[tentNofSC] = curNodeContracter.lgthSC[k];
                        altDiffSCall[tentNofSC] =  curNodeContracter.altSC[k];
                        tentNofSC++;
                    }
                    contractionPQ.add(new PQElement(edgeDiff, i+j));
                    candSCoffset[i+j + 1] = tentNofSC;

                    if ((i % (nofCandidates / 10 + 1)) == 0) {
                        System.err.print((10 * i / (nofCandidates / 10 + 1) + "% "));
                        System.err.print("(" + curNodeContracter.nofSC + "/" + edgeDiff + ") ");

                    }
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
                assert(stillIndep[candidates[curNode]]==true);
                contracted[candidates[curNode]] = true;
                // now copy its SCs in final SC list
                for (int i = candSCoffset[curNode]; i < candSCoffset[curNode + 1]; i++) {
                    assert(stillIndep[trgSCall[i]]==false);
                    assert(stillIndep[srcSCall[i]]==false);
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
            if ((!contracted[tempGraph.edgeSource(j)]) && (!contracted[tempGraph.edgeTarget[j]]))
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
                old2new[i] = newTempGraph.addNode(tempGraph.xCoord(i), tempGraph.yCoord(i), tempGraph.altNodeID(i), tempGraph.height(i), tempGraph.OSMID(i), 0);
            } else {
                old2new[i] = -1;
                assert (myCHGraph.level[tempGraph.altNodeID(i)] == Integer.MAX_VALUE);
                myCHGraph.level[tempGraph.altNodeID(i)] = newLevel;
            }


        /// / copy surviving edges to newTempGraph
        for (int j = 0; j < tempGraph.nofEdges(); j++) {
            int curSrc = tempGraph.edgeSource(j), curTrg = tempGraph.edgeTarget(j), curWgt = tempGraph.edgeWeight(j), curEdgeLength = tempGraph.edgeLength(j), curEdgeHeight = tempGraph.edgeAltitudeDifference(j);

            if ((!contracted[curSrc]) && (!contracted[curTrg]))    // copy edge to newTempGraph
            {
                assert(old2new[curSrc]!=-1);
                assert(old2new[curTrg]!=-1);
                newTempGraph.addEdge(old2new[curSrc], old2new[curTrg], curWgt, curEdgeLength, curEdgeHeight, -2, -2);
            }
        }

        // now add SC edges to newTempGRaph as well as myCHGraph
        System.err.println("No of added SCs: " + totalNofSC);
        for (int j = 0; j < totalNofSC; j++) {
            newTempGraph.addEdge(old2new[srcSCfinal[j]], old2new[trgSCfinal[j]], wgtSCfinal[j], lgthSCfinal[j], altDiffSCfinal[j], -2, -2);
            myCHGraph.addEdge(tempGraph.altNodeID(srcSCfinal[j]), tempGraph.altNodeID(trgSCfinal[j]), wgtSCfinal[j], lgthSCfinal[j], altDiffSCfinal[j], -2, -2);
        }

        newTempGraph.setupOffsets();

        tempGraph=newTempGraph.pruneGraphSelfloops();
        //tempGraph = newTempGraph.pruneGraph();
        // now create new tempGraph which consists of old graph (with edges between uncontracted nodes)
        // and all shortcuts created
        //
        return tempGraph.nofNodes();

    }


}

