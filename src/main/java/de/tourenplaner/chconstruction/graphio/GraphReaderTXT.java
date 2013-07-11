/*
 * (C) Copyright 2012 Peter Vollmer
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */

package de.tourenplaner.chconstruction.graphio;

import de.tourenplaner.chconstruction.graph.RAMGraph;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.regex.Pattern;

/**
 * User: Peter Vollmer
 * Date: 10/9/12
 * Time: 3:38 PM
 */
public class GraphReaderTXT implements GraphReader {

    private static final Pattern COMPILE = Pattern.compile(" ");

    private int calcWeight(int length, int roadType) {
        double weigth;
        switch (roadType) {
            //motorway
            case 1:
                weigth = (length * 1.3) / 1.3;
                break;
            //motorway link
            case 2:
                weigth = (length * 1.3) / 1.0;
                break;
            //primary
            case 3:
                weigth = (length * 1.3) / 0.7;
                break;
            //primary link
            case 4:
                weigth = (length * 1.3) / 0.7;
                break;
            //secondary
            case 5:
                weigth = (length * 1.3) / 0.65;
                break;
            //secondary link
            case 6:
                weigth = (length * 1.3) / 0.65;
                break;
            //tertiary
            case 7:
                weigth = (length * 1.3) / 0.6;
                break;
            //tertiary link
            case 8:
                weigth = (length * 1.3) / 0.6;
                break;
            //trunk
            case 9:
                weigth = (length * 1.3) / 0.8;
                break;
            //trunk link
            case 10:
                weigth = (length * 1.3) / 0.8;
                break;
            //unclassified
            case 11:
                weigth = (length * 1.3) / 0.25;
                break;
            //residential
            case 12:
                weigth = (length * 1.3) / 0.45;
                break;
            //living street
            case 13:
                weigth = (length * 1.3) / 0.3;
                break;
            //road
            case 14:
                weigth = (length * 1.3) / 0.25;
                break;
            //service
            case 15:
                weigth = (length * 1.3) / 0.3;
                break;
            //turning circle
            case 16:
                weigth = (length * 1.3) / 0.3;
                break;
            default:
                weigth = (length * 1.3) / 0.5;
        }


        return (int) weigth;
    }

    @Override
    public RAMGraph createRAMGraph(InputStream in) throws IOException {

        int nofNodes;
        int nofEdges;


        long curTime = System.currentTimeMillis();

        BufferedReader inb = new BufferedReader(new InputStreamReader(in));
        String line = inb.readLine();
        while (line != null && line.trim().startsWith("#")) {
            line = inb.readLine();
        }
        nofNodes = line != null ? Integer.parseInt(line) : 0;

        line = inb.readLine();
        nofEdges = line != null ? Integer.parseInt(line) : 0;

        RAMGraph graph = new RAMGraph(nofNodes, nofEdges, null);
        float x, y;
        int altID, OSMID, height;
        String[] splittedLine;
        for (int i = 0; i < nofNodes; i++) {
            splittedLine = COMPILE.split(inb.readLine());
            //altID[i]=i;
            if (splittedLine.length == 5) {
                altID = Integer.parseInt(splittedLine[0]);
                OSMID = Integer.parseInt(splittedLine[1]);
                x = Float.parseFloat(splittedLine[2]);
                y = Float.parseFloat(splittedLine[3]);
                height = Integer.parseInt(splittedLine[4]);
                graph.addNode(x, y, altID, height, OSMID);

                if ((i % (nofNodes / 10)) == 0) {
                    System.err.print((10 * i / (nofNodes / 10) + "% "));
                }
            } else {
                nofNodes--;
            }
        }

        int edgeSource, edgeTarget, edgeWeight, edgeLength, edgeHeight;
        for (int i = 0; i < nofEdges; i++) {
            splittedLine = COMPILE.split(inb.readLine());
            if (splittedLine.length == 4) {
                edgeSource = Integer.parseInt(splittedLine[0]);
                edgeTarget = Integer.parseInt(splittedLine[1]);
                edgeLength = Integer.parseInt(splittedLine[2]);
                edgeWeight = calcWeight(edgeLength, Integer.parseInt(splittedLine[3]));
                edgeHeight = graph.getHeight(edgeTarget)-graph.getHeight(edgeSource);
                if (edgeHeight < 0) {
                    edgeHeight = 0;
                }
                graph.addEdge(edgeSource, edgeTarget, edgeWeight, edgeLength, edgeHeight);

                if ((i % (nofEdges / 10)) == 0) {
                    System.err.print((10 * i / (nofEdges / 10) + "% "));
                }
            } else {
                nofEdges--;
            }
        }

        if (nofNodes != graph.nofNodes() || nofEdges != graph.nofEdges()) {
            graph = new RAMGraph(graph, nofNodes, nofEdges);
        }

        System.err.println("Parsing took " + (System.currentTimeMillis() - curTime));
        graph.setupOffsets();

        System.err.println("Read graph with " + nofNodes +
                " vertices and " + nofEdges + " edges in time " + (System.currentTimeMillis() - curTime) + "ms");
        return graph;
    }


}
