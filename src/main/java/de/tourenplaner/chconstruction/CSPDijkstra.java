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
    int internRunDijkstra(int src,int trg, int lambda, int maxLambda, boolean withoutSC){
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
        //printGeoPath(trg);
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
        if (!settled[trg]||lastSource == trg)
            return;
        int cur_node = trg;
        int cur_edge;
        System.out.println("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>");
        System.out.println("<gpx xmlns=\"http://www.topografix.com/GPX/1/1\" xmlns:gpxx=\"http://www.garmin.com/xmlschemas/GpxExtensions/v3\" xmlns:gpxtpx=\"http://www.garmin.com/xmlschemas/TrackPointExtension/v1\" creator=\"Oregon 400t\" version=\"1.1\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.topografix.com/GPX/1/1 http://www.topografix.com/GPX/1/1/gpx.xsd http://www.garmin.com/xmlschemas/GpxExtensions/v3 http://www.garmin.com/xmlschemas/GpxExtensionsv3.xsd http://www.garmin.com/xmlschemas/TrackPointExtension/v1 http://www.garmin.com/xmlschemas/TrackPointExtensionv1.xsd\">");
        System.out.println("  <trk>\n" + "    <name>Example GPX Document</name>");

        System.out.println("<trkseg>");
        do {
            cur_edge = pred(cur_node);
            System.out.println("<trkpt lat=\"" + myGraph.xCoord(myGraph.edgeTarget(cur_edge)) + "\" lon=\"" + myGraph.yCoord(myGraph.edgeTarget(cur_edge)) + "\"></trkpt>");
            cur_node = myGraph.edgeSource(cur_edge);
        } while (cur_node != lastSource);
        System.out.println("<trkpt lat=\"" + myGraph.xCoord(myGraph.edgeSource(cur_edge)) + "\" lon=\"" + myGraph.yCoord(myGraph.edgeSource(cur_edge)) + "\"></trkpt>");
        System.out.println("</trkseg>\n");
        System.out.println("</trk>\n</gpx>");
    }



}


