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
 * Date: 11/9/12
 * Time: 1:08 PM
 */
public abstract class Constructor {
    RAMGraph myCHGraph; // store the augmented result graph
    abstract int contractLevel(int newLevel);

}
