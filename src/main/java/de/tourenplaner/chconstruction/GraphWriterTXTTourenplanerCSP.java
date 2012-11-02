/*
 * (C) Copyright 2012 Peter Vollmer
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */

/*
 * (C) Copyright 2012 Peter Vollmer
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */

package de.tourenplaner.chconstruction;

import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;

/**
 * User: Peter Vollmer
 * Date: 10/10/12
 * Time: 5:38 PM
 */
public class GraphWriterTXTTourenplanerCSP implements GraphWriter {
    @Override
    public void writeRAMGraph(OutputStream out, RAMGraph ramGraph) throws IOException {
        long curTime = System.currentTimeMillis();

        OutputStreamWriter data_out = new OutputStreamWriter(out);

        // write-out nof vertices and edges
        data_out.write(ramGraph.nofNodes() + "\n");
        data_out.write(ramGraph.nofEdges() + "\n");
        // write-out coordinates and levels
        System.err.print("\n Nodes:");
        for (int i = 0; i < ramGraph.nofNodes(); i++) {
            data_out.write((int) (ramGraph.xCoord[i] * 10000000.0) + " ");
            data_out.write((int) (ramGraph.yCoord[i] * 10000000.0) + " ");
            data_out.write(ramGraph.height[i] + " ");
            data_out.write(ramGraph.level[i] + "\n");

            if ((i % (ramGraph.nofNodes / 10)) == 0) {
                System.err.print((10 * i / (ramGraph.nofNodes / 10) + "% "));
            }
        }
        System.err.print("\n Edges: ");
        // write-out edges
        for (int j = 0; j < ramGraph.nofEdges; j++) {
            data_out.write(ramGraph.edgeSource[j] + " ");
            data_out.write(ramGraph.edgeTarget[j] + " ");
            data_out.write(ramGraph.edgeWeight[j] + " ");
            data_out.write(ramGraph.edgeLength[j] + " ");
            data_out.write(ramGraph.edgeAltitudeDifference[j] + " ");
            data_out.write(ramGraph.edgeSkippedA[j] + " ");
            data_out.write(ramGraph.edgeSkippedB[j] + "\n");

            if ((j % (ramGraph.nofEdges / 10)) == 0) {
                System.err.print((10 * j / (ramGraph.nofEdges / 10) + "% "));
            }
        }
        data_out.close();
        System.err.println("Writing GTXT took " + (System.currentTimeMillis() - curTime) + "ms");
    }
}
