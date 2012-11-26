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

public class CSPCHDijkstra extends BDDijkstra {

    int[] predFwd;
    int[] predBwd;
    int[] pred;
    int bestNode;
    int src;

    CSPCHDijkstra(SGraph _myGraph) {
        super(_myGraph);
        predFwd = new int[myGraph.nofNodes()];
        predBwd = new int[myGraph.nofNodes()];
        pred = new int[myGraph.nofNodes()];
        for (int i = 0; i < myGraph.nofNodes(); i++) {
            predFwd[i]=predBwd[i]=pred[i]=-1;
        }
    }

    int runDijkstra(int src, int trg, int lambda, int maxLambda)    // returns distance from src to trg
    {
        this.src = src;
        // clean up previously touched nodes
        for (int i = 0; i < nofTouchedNodes; i++) {
            distFwd[touchedNodes[i]] = distBwd[touchedNodes[i]] = Integer.MAX_VALUE;
            settledFwd[touchedNodes[i]] = settledBwd[touchedNodes[i]] = false;
        }
        nofTouchedNodes = 0;
        myQueue.clear();
        // start with src
        labelFwd(src, 0,-1);
        labelBwd(trg, 0,-1);

        int bestDist = Integer.MAX_VALUE;
        bestNode = -1;

        int edgeCount = 0;
        // otherwise we have to process pq until settling trg
        boolean phaseFinished = false;
        while ((!myQueue.isEmpty()) && !phaseFinished) {
            //System.err.print(".");
            BDPQElement cur = myQueue.remove();
            int cur_dist = cur.key;
            int cur_node = cur.value;
            int cur_side = cur.queue;

            //if (cur_dist>bestDist)
            //	phaseFinished=true;

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
                        int tmp_resource = myGraph.edgeAltitudeDifference(tmp_edge);
                        int tmp_wgt = (myGraph.edgeWeight(tmp_edge)-tmp_resource)*lambda + maxLambda*tmp_resource;


                        int tmp_node = myGraph.edgeSource(tmp_edge);
                        if (distFwd[cur_node] - tmp_wgt > distFwd[tmp_node]) {
                            stalled = true;
                            break;
                        }
                        if (myGraph.level(tmp_node) < myGraph.level(cur_node))
                            break;
                    }

                    if ((settledBwd[cur_node]) && (distFwd[cur_node] + distBwd[cur_node] < bestDist)){
                        bestDist = distFwd[cur_node] + distBwd[cur_node];
                        bestNode = cur_node;
                    }

                    if (!stalled)
                        for (int i = myGraph.nofOutEdges(cur_node) - 1; i >= 0; i--) {
                            int cur_edge = myGraph.outEdgeID(cur_node, i);
                            int cur_trg = myGraph.edgeTarget(cur_edge);
                            int cur_weight = myGraph.edgeWeight(cur_edge);
                            int cur_resource = myGraph.edgeAltitudeDifference(cur_edge);
                            int new_Dist = cur_dist+ (cur_weight-cur_resource)*lambda + maxLambda*cur_resource;
                            if (myGraph.level(cur_trg) >= myGraph.level(cur_node))
                                edgeCount++;
                            else
                                break;
                            if ((myGraph.level(cur_trg) >= myGraph.level(cur_node)) && (distFwd[cur_trg] > new_Dist)) {
                                labelFwd(cur_trg, new_Dist, cur_edge);
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
                        int tmp_resource = myGraph.edgeAltitudeDifference(tmp_edge);
                        int tmp_wgt = (myGraph.edgeWeight(tmp_edge)-tmp_resource)*lambda + maxLambda*tmp_resource;
                        int tmp_node = myGraph.edgeTarget(tmp_edge);
                        if (distBwd[cur_node] - tmp_wgt > distBwd[tmp_node]) {
                            stalled = true;
                            break;
                        }
                        if (myGraph.level(cur_node) > myGraph.level(tmp_node))
                            break;
                    }


                    if ((settledFwd[cur_node]) && (distFwd[cur_node] + distBwd[cur_node] < bestDist)){
                        bestDist = distFwd[cur_node] + distBwd[cur_node];
                        bestNode = cur_node;
                    }
                    if (!stalled)
                        for (int i = myGraph.nofInEdges(cur_node) - 1; i >= 0; i--) {
                            int cur_edge = myGraph.inEdgeID(cur_node, i);
                            int cur_trg = myGraph.edgeSource(cur_edge);
                            int cur_weight = myGraph.edgeWeight(cur_edge);
                            int cur_resource = myGraph.edgeAltitudeDifference(cur_edge);
                            int new_Dist = cur_dist + (cur_weight-cur_resource)*lambda + maxLambda*cur_resource;
                            if (myGraph.level(cur_trg) >= myGraph.level(cur_node))
                                edgeCount++;
                            else
                                break;
                            if ((myGraph.level(cur_trg) >= myGraph.level(cur_node)) && (distBwd[cur_trg] > new_Dist)) {
                                labelBwd(cur_trg, new_Dist, cur_edge);
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
                        if (myGraph.level(curNode) > 80)
                            System.err.print(curNode + "(" + myGraph.level(curNode) + ") ");
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
        if (bestNode != -1){
            meldPredFwdBwd(predFwd,predBwd,bestNode);
        }
        //printGeoPath(trg);
        return bestDist;

    }

    private void meldPredFwdBwd (int[] predFwd, int[] predBwd, int bestNode){
        int curNode = bestNode;
        while (predFwd[curNode] != -1){
            pred[curNode] = predFwd[curNode];
            curNode = myGraph.edgeSource(predFwd[curNode]);
        }
        curNode = bestNode;
        while (predBwd[curNode] != -1){
            pred[myGraph.edgeTarget(predBwd[curNode])] = predBwd[curNode];
            curNode = myGraph.edgeTarget(predBwd[curNode]);
        }
    }


    void labelFwd(int v, int d, int p) {
        if ((distFwd[v] == Integer.MAX_VALUE) && (distBwd[v] == Integer.MAX_VALUE)) {
            touchedNodes[nofTouchedNodes++] = v;
        }
        distFwd[v] = d;
        predFwd[v] = p;
        myQueue.add(new BDPQElement(d, v, 0));
    }

    void labelBwd(int v, int d, int p) {
        if ((distFwd[v] == Integer.MAX_VALUE) && (distBwd[v] == Integer.MAX_VALUE)) {
            touchedNodes[nofTouchedNodes++] = v;
        }
        distBwd[v] = d;
        predBwd[v] = p;
        myQueue.add(new BDPQElement(d, v, 1));
    }

    int pred(int v)    // return predecessor of v in current SP-tree (-1 if none exists)
    {
        return pred[v];
    }


    void checkCHreqs() {
        for (int i = 0; i < myGraph.nofNodes(); i++) {
            for (int j = 0; j < myGraph.nofOutEdges(i) - 1; j++) {
                int cur_edge = myGraph.outEdgeID(i, j);
                int next_edge = myGraph.outEdgeID(i, j + 1);
                int trg_node = myGraph.edgeTarget(cur_edge);
                int next_trg_node = myGraph.edgeTarget(next_edge);
                assert (myGraph.level(trg_node) <= myGraph.level(next_trg_node));
            }
            for (int j = 0; j < myGraph.nofInEdges(i) - 1; j++) {
                int cur_edge = myGraph.inEdgeID(i, j);
                int next_edge = myGraph.inEdgeID(i, j + 1);
                int src_node = myGraph.edgeSource(cur_edge);
                int next_src_node = myGraph.edgeSource(next_edge);
                assert (myGraph.level(src_node) <= myGraph.level(next_src_node));
            }
        }
        System.err.println("CH Reqs ok!");
    }

    void printGeoPath(int trg) {
        if (bestNode == -1||src == trg)
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
        } while (cur_node != src);
        System.out.println("<trkpt lat=\"" + myGraph.xCoord(myGraph.edgeSource(cur_edge)) + "\" lon=\"" + myGraph.yCoord(myGraph.edgeSource(cur_edge)) + "\"></trkpt>");
        System.out.println("</trkseg>\n");
        System.out.println("</trk>\n</gpx>");
    }
}
