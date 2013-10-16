/*
 * (C) Copyright 2012 FMI Universität Stuttgart
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */
package de.tourenplaner.chconstruction.graphio;

import de.tourenplaner.chconstruction.graph.RAMGraph;
import fmi.graph.chgraph.Edge;
import fmi.graph.chgraph.Node;
import fmi.graph.definition.GraphException;
import fmi.graph.exceptions.InvalidFunctionException;
import fmi.graph.metaio.MetaData;
import fmi.graph.chgraph.Writer;

import java.io.IOException;
import java.io.OutputStream;

/**
 * @author Niklas Schnelle
 */
public class GraphWriterNewStandard implements GraphWriter {
    private boolean binary;

    public GraphWriterNewStandard(boolean binary) {
        this.binary = binary;
    }

    @Override
    public void writeRAMGraph(OutputStream out, RAMGraph ramGraph) throws IOException {
        long curTime = System.currentTimeMillis();
        // Write number of nodes and edges
        int numNodes = ramGraph.nofNodes();
        int numEdges = ramGraph.nofEdges();

        Writer w = new Writer();

        try {
            if (binary) {
                w.writeBin(out);
            } else {
                w.write(out);
            }
            w.setNodeCount(numNodes);
            w.setEdgeCount(numEdges);

            MetaData data = w.prepareMetaData(ramGraph.getMetaData());
            data.add("Origin", "CHConstructor");
            w.writeMetaData(data);

            for (int i = 0; i < numNodes; i++) {
                w.writeNode(new Node(i, ramGraph.getOSMID(i), ramGraph.getLat(i), ramGraph.getLon(i), ramGraph.getHeight(i), ramGraph.getLevel(i)));
            }

            for (int i = 0; i < numEdges; i++) {
                w.writeEdge(new Edge(ramGraph.getSource(i), ramGraph.getTarget(i), 0, ramGraph.getWeight(i),
                        ramGraph.getEuclidianLength(i), ramGraph.getSkippedA(i), ramGraph.getSkippedB(i)));
            }
            w.close();
        } catch (InvalidFunctionException e) {
            e.printStackTrace();
        } catch (GraphException ex) {
            throw new IOException(ex);
        }
        System.err.println("Writing graph took " + (System.currentTimeMillis() - curTime) + "ms");
    }
}
