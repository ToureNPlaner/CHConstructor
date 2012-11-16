/*
 * (C) Copyright 2012 Peter Vollmer
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */

package de.tourenplaner.chconstruction;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.Random;
import java.util.concurrent.DelayQueue;

/**
 * User: Peter Vollmer
 * Date: 11/9/12
 * Time: 1:25 PM
 */
public class CSPGraphInspector {
    RAMGraph graph;
    CSPDijkstra cspDijkstra;
    CSPCHDijkstra cspCHDijkstra;

    CSPGraphInspector (RAMGraph graph){
        this.graph = graph;
        cspDijkstra = new CSPDijkstra(graph);
        cspCHDijkstra = new CSPCHDijkstra(graph);
    }

    void inspectGraph (){
        Random random = new Random();
        int src;
        int trg;
        int lambda;
        int maxLambda = 16384;
        int counter = 0;

        int cspDist;
        int[] cspPath;

        int cspCHDist;
        int[] cspCHPath;

        while (true){

            counter++;
            System.err.println("Counter: "+counter);
            src = random.nextInt(graph.nofNodes());
            trg = random.nextInt(graph.nofNodes());
            lambda = random.nextInt(maxLambda);
            System.err.println("CSPGraphInspector: src="+src+", trg="+trg+", lambda="+lambda);
            cspDist = cspDijkstra.runDijkstraWithoutSC(src, trg, lambda, maxLambda);
            cspCHDist = cspCHDijkstra.runDijkstra(src,trg,lambda,maxLambda);

            if (cspDist == Integer.MAX_VALUE){
                if (cspCHDist == Integer.MAX_VALUE){
                    continue;
                } else {
                    System.err.println("csp finds no path from src("+src+") to trg("+trg+") and equals not to cspCH which found a path with dist("+cspCHDist+").");
                    return;
                }
            }
            if (cspCHDist == Integer.MAX_VALUE){
                System.err.println("cspCH finds no path from src("+src+") to trg("+trg+") and equals not to csp which found a path with dist("+cspDist+").");
                return;
            }
            if (cspDist != cspCHDist){
                System.err.println("cspDist("+cspDist+") equals not to cspCHDist("+cspCHDist+").");
                cspPath = cspBacktrack(src,trg);
                cspCHPath = cspCHBacktrack(src,trg);
                for(int i = 0 ; i < cspCHPath.length ; ++i) {
                    if (!(cspPath[i] == cspCHPath[i])){
                        System.err.println("On path element("+i+") cspPath("+cspPath[i]+")  equals not to cspCHPath("+cspCHPath[i]+")");
                        System.err.println("skippedAcsp "+graph.edgeSkippedA(cspPath[i])+" src "+graph.edgeSource(cspPath[i])+" weigth "+graph.edgeWeight(cspPath[i])+" altitude "+graph.edgeAltitudeDifference(cspPath[i]));
                        System.err.println("skippedAcspCH "+graph.edgeSkippedA(cspCHPath[i])+" src "+graph.edgeSource(cspPath[i])+" weigth "+graph.edgeWeight(cspCHPath[i])+" altitude "+graph.edgeAltitudeDifference(cspCHPath[i]));
                    }
                }
                return;
            }
            cspPath = cspBacktrack(src,trg);
            cspCHPath = cspCHBacktrack(src,trg);
            if (cspPath.length != cspCHPath.length){
                System.err.println("csp path elements("+cspPath.length+") equals not to csp path elements("+cspCHPath.length+").");
                return;
            }
            for(int i = 0 ; i < cspCHPath.length ; ++i) {
                if (!(cspPath[i] == cspCHPath[i])){
                    System.err.println("On path element("+i+") cspPath("+cspPath[i]+")  equals not to cspCHPath("+cspCHPath[i]+")");
                    return;
                }
            }
        }


    }

    private int[] cspBacktrack(int src, int trg) {
        int curNode = trg;
        int nofRouteEdges = 0;
        while (curNode != src) {
            nofRouteEdges++;
            curNode = (graph.edgeSource(cspDijkstra.pred(curNode)));
        }

        // Add them without values we set the values in the next step
        int [] path = new int[nofRouteEdges];

        // backtracking here
        int curEdge;
        curNode = trg;
        for (int i = nofRouteEdges-1 ; i >= 0; --i){
            curEdge = cspDijkstra.pred(curNode);
            path[i] = curEdge;
            curNode = graph.edgeSource(curEdge);
        }
        // add source node to the result.
        return path;

    }

    private int[] cspCHBacktrack(int src, int trg) {
        ArrayDeque<Integer> deque = new ArrayDeque<Integer>();
        ArrayDeque<Integer> pathQueue = new ArrayDeque<Integer>();
        int currNode = trg;
        int curEdge;
        int edgeSkippedA;
        int edgeSkippedB;

        while (currNode != src) {
            curEdge = cspCHDijkstra.pred(currNode);
            deque.addFirst(curEdge);
            currNode = graph.edgeSource(curEdge);
        }
        while (!deque.isEmpty()) {

            curEdge = deque.removeFirst();
            edgeSkippedA = graph.edgeSkippedA(curEdge);
            if (edgeSkippedA >= 0) {
                // We have a shortcut unpack it
                edgeSkippedB = graph.edgeSkippedB(curEdge);
                deque.addFirst(edgeSkippedB);
                deque.addFirst(edgeSkippedA);
            } else {
                // No shortcut remember it
                pathQueue.addLast(curEdge);
            }
        }

        int[] path = new int[pathQueue.size()];

        for (int i = 0 ; i < path.length; ++i){
           path[i] = pathQueue.removeFirst();
        }

        return path;

    }

}
