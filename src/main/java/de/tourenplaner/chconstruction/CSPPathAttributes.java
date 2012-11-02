/*
 * (C) Copyright 2012 Peter Vollmer
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; version 2
 * of the License.
 */

package de.tourenplaner.chconstruction;

/**
 * User: Peter Vollmer
 * Date: 11/2/12
 * Time: 10:11 AM
 */
public class CSPPathAttributes{
    private int dist;
    private int length;
    private int altitude;

    CSPPathAttributes(int dist, int length, int altitude){
        this.dist = dist;
        this.length = length;
        this.altitude = altitude;
    }

    int getDist(){
        return dist;
    }

    int getLength(){
        return length;
    }

    int getAltitude(){
        return altitude;
    }

}
