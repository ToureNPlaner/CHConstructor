/*
 * (C) Copyright 2012 Dr. Stefan Funke
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */

package de.tourenplaner.chconstruction;
/* Graph classes 
 * Desired Properties: 	- static
 * 						- explicit in- and outedge lists
 * 						- small memory footprints
 * Implementation:		- array based
 * 						- edges are stored once in a large list of edges
 * 						- two arrays inEdges/outEdges refer to this large list of edges
 *
 * RAM and ExtMEMORY
 */

public abstract class SGraph {

    int nofEdges;
    int nofNodes;

    SGraph() {
        nofEdges = nofNodes = 0;
    }


    // interface for accessing Graph data
    int nofEdges() {
        return nofEdges;
    }

    int nofNodes() {
        return nofNodes;
    }

    abstract float xCoord(int nodeID);

    abstract float yCoord(int nodeID);

    abstract int altNodeID(int nodeID);

    abstract int height(int nodeID);

    abstract int OSMID(int nodeID);

    abstract int level(int nodeID);


    abstract int nofOutEdges(int nodeID);

    abstract int nofInEdges(int nodeID);

    abstract int outEdgeID(int nodeID, int edgePos);    // returns edge ID of edgePos-th outEdge of nodeID

    abstract int inEdgeID(int nodeID, int edgePos);    // returns edge ID of edgePos-th inEdge of nodeID


    abstract int edgeSource(int edgeID);

    abstract int edgeTarget(int edgeID);

    abstract int edgeWeight(int edgeID);

    abstract int edgeLength(int edgeID);

    abstract int edgeHeight(int edgeID);

    abstract int edgeSkippedA(int edgeID);

    abstract int edgeSkippedB(int edgeID);
}
