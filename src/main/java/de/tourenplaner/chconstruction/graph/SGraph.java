/*
 * (C) Copyright 2012 Dr. Stefan Funke
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */

package de.tourenplaner.chconstruction.graph;
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

import fmi.graph.metaio.MetaData;

public abstract class SGraph {
    /**
     * access the MetaData object associated with this graph
     * @return
     */
    public abstract MetaData getMetaData();

    // interface for accessing Graph data
    public abstract int nofEdges();

    public abstract int nofNodes();


    public abstract float getLat(int nodeID);

    public abstract float getLon(int nodeID);

    public abstract int getAltNodeID(int nodeID);

    public abstract int getHeight(int nodeID);

    public abstract int getOSMID(int nodeID);

    public abstract int getLevel(int nodeID);


    public abstract int nofOutEdges(int nodeID);

    public abstract int nofInEdges(int nodeID);

    public abstract int outEdgeID(int nodeID, int edgePos);    // returns edge ID of edgePos-th outEdge of nodeID

    public abstract int inEdgeID(int nodeID, int edgePos);    // returns edge ID of edgePos-th inEdge of nodeID


    public abstract int getSource(int edgeID);

    public abstract int getTarget(int edgeID);

    public abstract int getWeight(int edgeID);

    public abstract int getEuclidianLength(int edgeID);

    public abstract int getAltitudeDifference(int edgeID);

    public abstract int getSkippedA(int edgeID);

    public abstract int getSkippedB(int edgeID);
}
