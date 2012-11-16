/*
 * (C) Copyright 2012 Peter Vollmer
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

public class CSPDijkstra {
    SGraph myGraph;
    int[] dist;
    int[] pred;
    boolean[] settled;
    PriorityQueue<PQElement> myQueue;
    int[] touchedNodes;
    int nofTouchedNodes;
    int lastSource;

    CSPDijkstra(SGraph _myGraph) {
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
        if (v < 0){
            System.err.println("CSPDijkstra v: "+v);
        }
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
    // writes the result of distance and altitudeSum from src to trg into returnLength and returnAltitude
    int runDijkstra(int src, int trg, int lambda, int maxLambda)
    {
         return internRunDijkstra(src,trg,lambda,maxLambda,false);
    }

    int runDijkstraWithoutSC(int src, int trg, int lambda, int maxLambda)
    {
        return internRunDijkstra(src,trg,lambda, maxLambda,true);
    }

    private int internRunDijkstra (int src, int trg, int lambda, int maxLambda, boolean withoutSC){
        //if ((lastSource != src)) {
        // clean up previously touched nodes
        //System.err.println("CSPDijkstra: Src= "+src+", Trg= "+trg);
        for (int i = 0; i < nofTouchedNodes; i++) {
            dist[touchedNodes[i]] = Integer.MAX_VALUE;
            settled[touchedNodes[i]] = false;
        }
        nofTouchedNodes = 0;
        myQueue.clear();
        // start with src
        label(src, 0,  -1);
        lastSource = src;
        //} else if (settled[trg]) {
        //    return new CSPPathAttributes(dist[trg], weightSum[trg], altitudeSum[trg]);
        //}
        // otherwise we have to process pq until settling trg                     CSPPathAttributes
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
                    if (withoutSC){
                        if (myGraph.edgeSkippedA(cur_edge) >= 0 ){
                            continue;
                        }
                    }
                    int cur_trg = myGraph.edgeTarget(cur_edge);
                    int cur_weight = myGraph.edgeWeight(cur_edge);
                    int cur_altitude = myGraph.edgeAltitudeDifference(cur_edge);
                    int newDist = cur_dist + (cur_weight-cur_altitude)*lambda + maxLambda*cur_altitude;
                    if (dist[cur_trg] > newDist) {
                        label(cur_trg, newDist , cur_edge);
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


