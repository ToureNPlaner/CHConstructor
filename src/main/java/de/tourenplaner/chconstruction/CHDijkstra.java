/*
 * (C) Copyright 2012 Dr. Stefan Funke
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */

package de.tourenplaner.chconstruction;

import de.tourenplaner.chconstruction.graph.SGraph;

public class CHDijkstra extends BDDijkstra {

    CHDijkstra(SGraph _myGraph) {
        super(_myGraph);
    }

    int runDijkstra(int src, int trg)    // returns distance from src to trg
    {

        // clean up previously touched nodes
        for (int i = 0; i < nofTouchedNodes; i++) {
            distFwd[touchedNodes[i]] = distBwd[touchedNodes[i]] = Integer.MAX_VALUE;
            settledFwd[touchedNodes[i]] = settledBwd[touchedNodes[i]] = false;
        }
        nofTouchedNodes = 0;
        myQueue.clear();
        // start with src
        labelFwd(src, 0);
        labelBwd(trg, 0);

        int bestDist = Integer.MAX_VALUE;

        int edgeCount = 0;
        // otherwise we have to process pq until settling trg
        boolean phaseFinished = false;
        while ((!myQueue.isEmpty()) && !phaseFinished) {
            //System.err.print(".");
            BDPQElement cur = myQueue.remove();
            int cur_dist = cur.key;
            int cur_node = cur.value;
            int cur_side = cur.queue;

            if (cur_dist>bestDist)
            	phaseFinished=true;

            if (cur_side == 0)    // we are in forward search
            {
                if (cur_dist == distFwd[cur_node]) {
                    settledFwd[cur_node] = true;

                    boolean stalled = false;

                    // check for stalling (if there is a node tmp_node (ABOVE) and an edge (tmp_node,cur_node)
                    // which sum to a smaller distance (!)
                    //for(int j=0; j<myGraph.nofInEdges(cur_node); j++)
                    for (int j = myGraph.nofInEdges(cur_node) - 1; j >= 0; j--) {
                        int tmp_edge = myGraph.inEdgeID(cur_node, j);
                        int tmp_wgt = myGraph.getWeight(tmp_edge);
                        int tmp_node = myGraph.getSource(tmp_edge);
                        if (distFwd[cur_node] - tmp_wgt > distFwd[tmp_node]) {
                            stalled = true;
                            break;
                        }
                        if (myGraph.getLevel(tmp_node) < myGraph.getLevel(cur_node))
                            break;
                    }

                    if ((settledBwd[cur_node]) && (distFwd[cur_node] + distBwd[cur_node] < bestDist))
                        bestDist = distFwd[cur_node] + distBwd[cur_node];

                    if (!stalled)
                        for (int i = myGraph.nofOutEdges(cur_node) - 1; i >= 0; i--) {
                            int cur_edge = myGraph.outEdgeID(cur_node, i);
                            int cur_trg = myGraph.getTarget(cur_edge);
                            int cur_weight = myGraph.getWeight(cur_edge);
                            if (myGraph.getLevel(cur_trg) >= myGraph.getLevel(cur_node))
                                edgeCount++;
                            else
                                break;
                            if (distFwd[cur_trg] > cur_dist + cur_weight) { // Levels have already been checked
                                labelFwd(cur_trg, cur_dist + cur_weight);
                            }
                        }
                }
            } else    // we are in backward search
            {
                if (cur_dist == distBwd[cur_node]) {
                    settledBwd[cur_node] = true;
                    boolean stalled = false;

                    // check for stalling: if there is a node ABOVE cur_node ...
                    //for(int j=0; j<myGraph.nofOutEdges(cur_node); j++)
                    for (int j = myGraph.nofOutEdges(cur_node) - 1; j >= 0; j--) {
                        int tmp_edge = myGraph.outEdgeID(cur_node, j);
                        int tmp_wgt = myGraph.getWeight(tmp_edge);
                        int tmp_node = myGraph.getTarget(tmp_edge);
                        if (distBwd[cur_node] - tmp_wgt > distBwd[tmp_node]) {
                            stalled = true;
                            break;
                        }
                        if (myGraph.getLevel(cur_node) > myGraph.getLevel(tmp_node))
                            break;
                    }


                    if ((settledFwd[cur_node]) && (distFwd[cur_node] + distBwd[cur_node] < bestDist))
                        bestDist = distFwd[cur_node] + distBwd[cur_node];

                    if (!stalled)
                        for (int i = myGraph.nofInEdges(cur_node) - 1; i >= 0; i--) {
                            int cur_edge = myGraph.inEdgeID(cur_node, i);
                            int cur_trg = myGraph.getSource(cur_edge);
                            int cur_weight = myGraph.getWeight(cur_edge);
                            if (myGraph.getLevel(cur_trg) >= myGraph.getLevel(cur_node))
                                edgeCount++;
                            else
                                break;
                            if (distBwd[cur_trg] > cur_dist + cur_weight) { // Levels have already been checked
                                labelBwd(cur_trg, cur_dist + cur_weight);
                            }
                        }
                }
            }
        }
        System.err.println("CH-BD has touched " + nofTouchedNodes + " and looked at " + edgeCount + " edges");

        if (1 == 0) {
            // ONLY FOR DEBUGGING
            // now check how many nodes bear correct distances
            int bestTmpDist = Integer.MAX_VALUE;
            Dijkstra myTmpDijkstraF = new Dijkstra(myGraph);
            Dijkstra myTmpDijkstraB = new BwDijkstra(myGraph);
            int fwdOK = 0, fwdTotal = 0, bwdOK = 0, bwdTotal = 0, bothOK = 0;

            for (int i = 0; i < nofTouchedNodes; i++) {
                //if ((i%100)==0)
                //	System.err.println(i+" ");
                int curNode = touchedNodes[i];
                boolean fwdUse = false, bwdUse = false;

                int fwdDist = distFwd[curNode];
                if (fwdDist != Integer.MAX_VALUE) {
                    fwdTotal++;
                    myTmpDijkstraF.runDijkstra(src, curNode);
                    assert (myTmpDijkstraF.dist[curNode] != Integer.MAX_VALUE);
                    if (myTmpDijkstraF.dist[curNode] == fwdDist) {
                        if (myGraph.getLevel(curNode) > 80)
                            System.err.print(curNode + "(" + myGraph.getLevel(curNode) + ") ");
                        fwdOK++;
                        fwdUse = true;
                    }
                }
                int bwdDist = distBwd[curNode];
                if (bwdDist != Integer.MAX_VALUE) {
                    bwdTotal++;
                    myTmpDijkstraB.runDijkstra(trg, curNode);
                    assert (myTmpDijkstraB.dist[curNode] != Integer.MAX_VALUE);
                    if (myTmpDijkstraB.dist[curNode] == bwdDist) {
                        bwdOK++;
                        bwdUse = true;
                    }
                }
                if ((bwdDist != Integer.MAX_VALUE) && (fwdDist != Integer.MAX_VALUE)) {
                    if (myTmpDijkstraF.dist[curNode] + myTmpDijkstraB.dist[curNode] < bestTmpDist)
                        bestTmpDist = myTmpDijkstraF.dist[curNode] + myTmpDijkstraB.dist[curNode];
                }
                if (bwdUse || fwdUse)
                    bothOK++;
            }
            System.err.println("\n Forward search: " + fwdOK + "/" + fwdTotal);
            System.err.println("Backward search: " + bwdOK + "/" + bwdTotal);
            System.err.println("Best Distance:" + bestTmpDist + " and total nodes: " + bothOK);
        }
        return bestDist;

    }

    void checkCHreqs() {
        for (int i = 0; i < myGraph.nofNodes(); i++) {
            for (int j = 0; j < myGraph.nofOutEdges(i) - 1; j++) {
                int cur_edge = myGraph.outEdgeID(i, j);
                int next_edge = myGraph.outEdgeID(i, j + 1);
                int trg_node = myGraph.getTarget(cur_edge);
                int next_trg_node = myGraph.getTarget(next_edge);
                assert (myGraph.getLevel(trg_node) <= myGraph.getLevel(next_trg_node));
            }
            for (int j = 0; j < myGraph.nofInEdges(i) - 1; j++) {
                int cur_edge = myGraph.inEdgeID(i, j);
                int next_edge = myGraph.inEdgeID(i, j + 1);
                int src_node = myGraph.getSource(cur_edge);
                int next_src_node = myGraph.getSource(next_edge);
                assert (myGraph.getLevel(src_node) <= myGraph.getLevel(next_src_node));
            }
        }
        System.err.println("CH Reqs ok!");
    }
}
