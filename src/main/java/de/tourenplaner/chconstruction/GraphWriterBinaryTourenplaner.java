/*
 * (C) Copyright 2012 Peter Vollmer
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */

package de.tourenplaner.chconstruction;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.Charset;

/**
 * User: Peter Vollmer
 * Date: 10/8/12
 * Time: 4:20 PM
 */
public class GraphWriterBinaryTourenplaner implements GraphWriter{

    private static final int version = 2;

    @Override
    public void writeRAMGraph(OutputStream out, RAMGraph ramGraph) throws IOException {

        //TODO

        // 8388608 Bytes = 8 MB
        BufferedOutputStream bout = new BufferedOutputStream(out, 8388608);
        DataOutputStream dout = new DataOutputStream(bout);
        // Write magic bytes so we can identify the ToureNPlaner Graphfileformat
        dout.write("TPG\n".getBytes(Charset.forName("UTF-8")));
        // Write format version
        dout.writeInt(version);

        // Write number of nodes and edges
        int numNodes = ramGraph.nofNodes();
        int numEdges = ramGraph.nofEdges();
        dout.writeInt(numNodes);
        dout.writeInt(numEdges);

        // Write nodes
        for (int i = 0; i < numNodes; i++) {
            dout.writeInt((int)(ramGraph.xCoord(i)*10000000.0));
            dout.writeInt((int)(ramGraph.yCoord(i)*10000000.0));
            dout.writeInt(ramGraph.height(i));
            dout.writeInt(ramGraph.level(i));
        }

        // Write edges
        for (int i = 0; i < numEdges; i++) {
            dout.writeInt(ramGraph.edgeSource(i));
            dout.writeInt(ramGraph.edgeTarget(i));
            dout.writeInt(ramGraph.edgeWeight(i));
            dout.writeInt(ramGraph.edgeWeight(i));
            dout.writeInt(ramGraph.edgeSkippedA(i));
            dout.writeInt(ramGraph.edgeSkippedB(i));
        }

        bout.flush();
        dout.flush();
        dout.close();
        //log.info("Successfully wrote graph");
    }
}
