/*
 * (C) Copyright 2012 Stefan Funke
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */

package de.tourenplaner.chconstruction;

import de.tourenplaner.chconstruction.graph.SGraph;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;

/**
 * @author Stefan Funke, Niklas Schnelle
 */
public class CoreExporter {
    public void exportCore(SGraph g, OutputStream os, int lvlBound) {
        long curTime = System.currentTimeMillis();
        try {
            BufferedWriter data_out = new BufferedWriter(new OutputStreamWriter(os));

            int highNodes = 0;
            int highEdges = 0;
            int coreEdges = 0;

            for (int i = 0; i < g.nofNodes(); i++)
                if (g.getLevel(i) >= lvlBound) highNodes++;

            for (int j = 0; j < g.nofEdges(); j++) {
                int src = g.getSource(j);
                int trg = g.getTarget(j);
                if ((g.getLevel(src) >= lvlBound) && (g.getLevel(trg) >= lvlBound)) {
                    // but we dont want to write out edges that were created as shortcuts
                    // very late (contracting a node that had very high  level !)
                    boolean outputEdge = true;
                    int skipEdgeA = g.getSkippedA(j);
                    if (skipEdgeA >= 0) {
                        int contrNode = g.getTarget(skipEdgeA);
                        if (g.getLevel(contrNode) >= lvlBound) outputEdge = false;
                    }
                    highEdges++;
                    if (outputEdge) coreEdges++;

                }
            }
            // write-out nof vertices and edges
            System.err.println("We have a CORE of size  " + highNodes + " / " + coreEdges + " vs. " + highEdges + " highEdges");
            data_out.write(coreEdges + "\n");
            // write-out coordinates and levels

            System.err.print("\n Edges: ");
            // write-out edges
            for (int j = 0; j < g.nofEdges(); j++) {
                int src = g.getSource(j);
                int trg = g.getTarget(j);
                if ((g.getLevel(src) >= lvlBound) && (g.getLevel(trg) >= lvlBound)) {
                    boolean outputEdge = true;
                    int skipEdgeA = g.getSkippedA(j);
                    if (skipEdgeA >= 0) {
                        int contrNode = g.getTarget(skipEdgeA);
                        if (g.getLevel(contrNode) >= lvlBound) outputEdge = false;
                    }

                    if (outputEdge) {
                        data_out.write(g.getSource(j) + " ");
                        data_out.write(g.getTarget(j) + " ");
                        data_out.write(g.getWeight(j) + " ");
                        data_out.write(g.getEuclidianLength(j) + " ");
                        data_out.write("\n");
                    }
                }
                if ((j % (g.nofEdges() / 10)) == 0) {
                    System.err.print((10 * j / (g.nofEdges() / 10) + "% "));
                }
            }
            data_out.close();
        } catch (IOException e) {
            System.err.println("IO exception = " + e);
        }
        System.err.println("Writing GTXT took  " + (System.currentTimeMillis() - curTime) + "ms");
    }
}
