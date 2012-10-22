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


// functionality:
// * compute shortest path between given source and all nodes closer than target
// * no special handling of memory or CH stuff
// * implements ToClean list and source-history, though
// * baseline algorithm, probably not used in practice later on

public class Dijkstra {
    SGraph myGraph;
    int[] dist;
    int[] pred;
    boolean[] settled;
    PriorityQueue<PQElement> myQueue;
    int[] touchedNodes;
    int nofTouchedNodes;
    int lastSource;

    Dijkstra(SGraph _myGraph) {
        myGraph = _myGraph;
        dist = new int[myGraph.nofNodes()];
        pred = new int[myGraph.nofNodes()];
        settled = new boolean[myGraph.nofNodes()];
        touchedNodes = new int[myGraph.nofNodes()];
        for (int i = 0; i < myGraph.nofNodes(); i++) {
            dist[i] = Integer.MAX_VALUE;
            pred[i] = -1;
            settled[i] = false;
        }
        myQueue = new PriorityQueue<PQElement>();
        lastSource = -1;
        nofTouchedNodes = 0;
    }

    void label(int v, int d, int p) {
        if (dist[v] == Integer.MAX_VALUE) {
            touchedNodes[nofTouchedNodes++] = v;
        }
        dist[v] = d;
        pred[v] = p;
        myQueue.add(new PQElement(d, v));
    }

    int pred(int v)    // return predecessor of v in current SP-tree (-1 if none exists)
    {
        return pred[v];
    }

    int runDijkstra(int src, int trg)    // returns distance from src to trg
    {
        if ((lastSource != src)) {
            // clean up previously touched nodes
            for (int i = 0; i < nofTouchedNodes; i++) {
                dist[touchedNodes[i]] = Integer.MAX_VALUE;
                settled[touchedNodes[i]] = false;
            }
            nofTouchedNodes = 0;
            myQueue.clear();
            // start with src
            label(src, 0, -1);
            lastSource = src;
        } else if (settled[trg]) {
            return dist[trg];
        }
        // otherwise we have to process pq until settling trg
        boolean targetFound = false;
        while ((!myQueue.isEmpty()) && !targetFound) {
            PQElement cur = myQueue.remove();
            int cur_dist = cur.key;
            int cur_node = cur.value;

            if (cur_dist == dist[cur_node]) {
                settled[cur_node] = true;
                if (cur_node == trg)
                    targetFound = true;
                for (int i = 0; i < myGraph.nofOutEdges(cur_node); i++) {
                    int cur_edge = myGraph.outEdgeID(cur_node, i);
                    int cur_trg = myGraph.edgeTarget(cur_edge);
                    int cur_weight = myGraph.edgeWeight(cur_edge);
                    if (dist[cur_trg] > cur_dist + cur_weight) {
                        label(cur_trg, cur_dist + cur_weight, cur_node);
                    }
                }
            }
        }
        // System.err.println("Dijkstra has touched "+nofTouchedNodes);
        return dist[trg];

    }

    void printPath(int trg) {
        if (!settled[trg])
            return;
        int cur_node = trg;
        do {
            System.err.print(cur_node + "(" + myGraph.level(cur_node) + ")-");
            cur_node = pred[cur_node];
        } while (cur_node != lastSource);
        System.err.println(cur_node + "(" + myGraph.level(cur_node) + ")");

    }

    void printGeoPath(int trg) {
        if (!settled[trg])
            return;
        int cur_node = trg;
        System.err.println("***********************************************");
        do {
            System.err.println(myGraph.yCoord(cur_node) + "," + myGraph.xCoord(cur_node));
            cur_node = pred[cur_node];
        } while (cur_node != lastSource);
        System.err.println(myGraph.yCoord(cur_node) + "," + myGraph.xCoord(cur_node));
        System.err.println("***********************************************");
    }

}


