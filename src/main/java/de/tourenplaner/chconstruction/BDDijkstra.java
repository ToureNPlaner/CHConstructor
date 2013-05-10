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

import java.util.PriorityQueue;

// functionality:
// * compute shortest path between given source and target
// * run bidirectional Dijkstra from source and target (on reverse edges)
// * termination:
// 		let D be the min-distance from the PQ when we first settle a node
//		that has been settled from the opposite direction => tentative dist. U
//		we go over all settled nodes of Fwd and take the minimum of
//		distFwd+distBwd

public class BDDijkstra {
    SGraph myGraph;
    int[] distFwd;
    int[] distBwd;
    boolean[] settledFwd;
    boolean[] settledBwd;
    PriorityQueue<BDPQElement> myQueue;
    int[] touchedNodes;
    int nofTouchedNodes;

    BDDijkstra(SGraph _myGraph) {
        myGraph = _myGraph;
        distFwd = new int[myGraph.nofNodes()];
        distBwd = new int[myGraph.nofNodes()];
        settledFwd = new boolean[myGraph.nofNodes()];
        settledBwd = new boolean[myGraph.nofNodes()];
        touchedNodes = new int[myGraph.nofNodes()];
        for (int i = 0; i < myGraph.nofNodes(); i++) {
            distFwd[i] = distBwd[i] = Integer.MAX_VALUE;
            settledFwd[i] = settledBwd[i] = false;
        }
        myQueue = new PriorityQueue<BDPQElement>();
        nofTouchedNodes = 0;
    }

    void labelFwd(int v, int d) {
        if ((distFwd[v] == Integer.MAX_VALUE) && (distBwd[v] == Integer.MAX_VALUE)) {
            touchedNodes[nofTouchedNodes++] = v;
        }
        distFwd[v] = d;
        myQueue.add(new BDPQElement(d, v, 0));
    }

    void labelBwd(int v, int d) {
        if ((distFwd[v] == Integer.MAX_VALUE) && (distBwd[v] == Integer.MAX_VALUE)) {
            touchedNodes[nofTouchedNodes++] = v;
        }
        distBwd[v] = d;
        myQueue.add(new BDPQElement(d, v, 1));
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
                    if (settledBwd[cur_node])        // ONLY WITHOUT CH !!!
                        phaseFinished = true;
                    for (int i = 0; i < myGraph.nofOutEdges(cur_node); i++) {
                        int cur_edge = myGraph.outEdgeID(cur_node, i);
                        int cur_trg = myGraph.getTarget(cur_edge);
                        int cur_weight = myGraph.getWeight(cur_edge);
                        if ((myGraph.getLevel(cur_trg) >= myGraph.getLevel(cur_node)) && (distFwd[cur_trg] > cur_dist + cur_weight)) {
                            labelFwd(cur_trg, cur_dist + cur_weight);
                            if (settledBwd[cur_trg] && (distFwd[cur_trg] + distBwd[cur_trg] < bestDist))
                                bestDist = distFwd[cur_trg] + distBwd[cur_trg];
                        }
                    }
                }
            } else    // we are in backward search
            {
                if (cur_dist == distBwd[cur_node]) {
                    settledBwd[cur_node] = true;
                    if (settledFwd[cur_node])
                        phaseFinished = true;
                    for (int i = 0; i < myGraph.nofInEdges(cur_node); i++) {
                        int cur_edge = myGraph.inEdgeID(cur_node, i);
                        int cur_trg = myGraph.getSource(cur_edge);
                        int cur_weight = myGraph.getWeight(cur_edge);
                        if ((myGraph.getLevel(cur_trg) >= myGraph.getLevel(cur_node)) && (distBwd[cur_trg] > cur_dist + cur_weight)) {
                            labelBwd(cur_trg, cur_dist + cur_weight);
                            if (settledFwd[cur_trg] && (distFwd[cur_trg] + distBwd[cur_trg] < bestDist))
                                bestDist = distFwd[cur_trg] + distBwd[cur_trg];
                        }
                    }
                }
            }
            //System.err.println(myQueue.size());
        }
        // now we iterate over all

        //System.err.println("BD has touched "+nofTouchedNodes);
        /*
          bestDist=Integer.MAX_VALUE;
          for(int i=0; i<nofTouchedNodes; i++)
          {
              int cur_node=touchedNodes[i];
              if (settledFwd[cur_node] && settledBwd[cur_node] &&
                      (distFwd[cur_node]+distBwd[cur_node]<bestDist))
                  bestDist=distFwd[cur_node]+distBwd[cur_node];
          }
          */
        return bestDist;

    }
}

class BDPQElement implements Comparable<BDPQElement> {

    public int key;
    public int value;
    public int queue;

    BDPQElement(int a, int b, int c) {
        key = a;
        value = b;
        queue = c;
    }

    public int compareTo(BDPQElement o) {
        if (key > o.key) return 1;
        else if (key == o.key) return 0;
        else return -1;
    }
}
