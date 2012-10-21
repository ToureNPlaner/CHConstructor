/*
 * (C) Copyright 2012 Peter Vollmer
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */

package de.tourenplaner.chconstruction;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.regex.Pattern;

/**
 * User: Peter Vollmer
 * Date: 10/8/12
 * Time: 12:43 PM
 */
public class GraphReaderTXTFunke implements GraphReader {
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

        RAMGraph graph = new RAMGraph(nofNodes, nofEdges);
        float x, y;
        int altID;
        String[] splittedLine;
        for (int i = 0; i < nofNodes; i++) {
            splittedLine = COMPILE.split(inb.readLine());
            x = Float.parseFloat(splittedLine[1]);
            y = Float.parseFloat(splittedLine[2]);
            //altID[i]=i;
            altID = Integer.parseInt(splittedLine[3]);
            graph.addNode(x, y, altID, -1, -1);

            if ((i % (nofNodes / 10)) == 0) {
                System.err.print((10 * i / (nofNodes / 10) + "% "));
            }
        }

        int edgeSource, edgeTarget, edgeWeight, edgeHeight;
        for (int i = 0; i < nofEdges; i++) {
            splittedLine = COMPILE.split(inb.readLine());
            edgeSource = Integer.parseInt(splittedLine[0]);
            edgeTarget = Integer.parseInt(splittedLine[1]);
            edgeWeight = Integer.parseInt(splittedLine[2]);
            edgeHeight = graph.height(edgeTarget)-graph.height(edgeSource);
            if (edgeHeight < 0) {
                edgeHeight = 0;
            }
            graph.addEdge(edgeSource, edgeTarget, edgeWeight, edgeWeight, edgeHeight);

            if ((i % (nofEdges / 10)) == 0) {
                System.err.print((10 * i / (nofEdges / 10) + "% "));
            }
        }

        System.err.println("Parsing took " + (System.currentTimeMillis() - curTime));
        graph.setupOffsets();

        System.err.println("Read graph with " + nofNodes +
                " vertices and " + nofEdges + " edges in time " + (System.currentTimeMillis() - curTime) + "ms");
        return graph;
    }

    private static final Pattern COMPILE = Pattern.compile(" ");
}
